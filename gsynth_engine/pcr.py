"""Conventional and restriction-cloning PCR primer design."""
from __future__ import annotations

from dataclasses import dataclass, field

from gsynth_engine.cloning import Digest, digest_linear, find_sites, translate
from gsynth_engine.constants import (
    ALL_ENZYMES,
    left_remainders,
    supplies_start_codon,
)
from gsynth_engine.sequence import (
    SequenceError,
    clean_dna,
    gc_content,
    longest_homopolymer,
    reverse_complement,
)
from gsynth_engine.thermo import PCR, melting_temperature

START_CODON_MODES = ("use_site", "keep_both")

#: Annealing-region length bounds.
MIN_ANNEAL: int = 18
MAX_ANNEAL: int = 30

#: Target melting temperature and accepted deviation.
TARGET_TM: float = 60.0
TM_TOLERANCE: float = 3.0

#: Maximum accepted Tm difference within a primer pair.
MAX_TM_DIFFERENCE: float = 5.0

#: Terminal clamp length used for restriction-site cleavage.
DEFAULT_CLAMP: int = 6

#: Validated clamp sequence without supported restriction sites.
_CLAMP_BASES: str = "GGAGGTGAAGACAGTAACTG"

#: Annealing-temperature offset from the lower primer Tm.
TA_OFFSET: float = 5.0

#: Maximum proposed annealing temperature.
MAX_TA: float = 72.0


def _clamp_of(length: int) -> str:
    """Clamp bases for a tail, taken from the validated terminal run."""
    if length <= 0:
        return ""
    return _CLAMP_BASES[:length]


@dataclass(frozen=True)
class PcrPrimer:
    """PCR primer with separate 5′ tail and template-annealing region."""

    name: str
    sequence: str
    #: The 5' addition. Empty for conventional PCR.
    tail: str
    #: The 3' portion that binds the template.
    anneals: str
    direction: int          #: 1 forward, -1 reverse
    #: Footprint of the annealing portion on the template, top-strand
    #: coordinates, half-open and 0-based. A reverse primer's footprint is
    #: still written left-to-right on the top strand.
    start: int
    end: int
    #: Tm of the annealing portion under PCR conditions. This sets Ta.
    tm: float
    #: Tm of the whole oligo, which applies from the third cycle onwards.
    tm_full: float
    gc: float
    enzyme: str | None = None
    restriction_site: str = ""
    warnings: tuple[str, ...] = ()

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def anneal_length(self) -> int:
        return len(self.anneals)

    @property
    def has_gc_clamp(self) -> bool:
        """A G or C among the last two bases.

        The 3' end is the end the polymerase extends from; a run of A/T there
        breathes open and primes poorly. Three or more G/C is the opposite
        problem and is warned about separately.
        """
        return any(base in "GC" for base in self.anneals[-2:])

    @property
    def as_row(self) -> dict[str, object]:
        return {
            "Name": self.name,
            "Sequence (5'->3')": self.sequence,
            "Length (nt)": self.length,
            "Tm anneal (°C)": self.tm,
            "Tm full (°C)": self.tm_full,
            "GC (%)": self.gc,
            "Tail": self.tail or "—",
            "Site": self.enzyme or "—",
        }


@dataclass
class PcrResult:
    """A designed reaction: the two primers, the product, and the verdict."""

    forward: PcrPrimer
    reverse: PcrPrimer
    #: The amplicon's top strand, tails included.
    product: str
    #: The stretch of template between the primers' outer edges, without
    #: tails — what conventional PCR would have produced.
    amplified_region: str
    template_start: int
    template_end: int
    #: Suggested annealing temperature, from the lower annealing Tm.
    annealing_temperature: float
    #: Present only for a cloning design.
    left_enzyme: str | None = None
    right_enzyme: str | None = None
    digest: Digest | None = None
    #: Where the reading frame starts in the digested insert, when asked for.
    insert_orf_start: int | None = None
    problems: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def product_length(self) -> int:
        return len(self.product)

    @property
    def is_clean(self) -> bool:
        """No problems — this reaction can be run and its product cut.

        Warnings are not consulted: a primer with a middling GC clamp is
        worth knowing about and is still worth ordering, whereas an enzyme
        that cuts inside the gene is not a reaction to optimise.
        """
        return not self.problems

    @property
    def insert(self) -> str:
        """The insert as it goes into the ligation."""
        return self.digest.top if self.digest else self.amplified_region


def _tm(sequence: str) -> float:
    """Tm under the PCR reaction, rounded the way the interface shows it."""
    return round(melting_temperature(sequence, conditions=PCR), 1)


def _self_complementarity(sequence: str) -> int:
    """Longest run where the oligo's 3' end pairs with itself elsewhere.

    A primer whose 3' end is complementary to its own body forms a hairpin
    or a primer-dimer, and both consume primer that should be extending the
    template. Only the 3' end is checked, because that is the end that gets
    extended — an internal match stalls but does not produce a product.
    """
    seq = sequence.upper()
    tail_3 = seq[-6:]
    if len(tail_3) < 4:
        return 0
    partner = reverse_complement(tail_3)
    best = 0
    for size in range(4, len(partner) + 1):
        probe = partner[-size:]
        # Exclude the 3' end matching itself, which is not a real duplex.
        if probe in seq[: len(seq) - size]:
            best = size
    return best


def _cross_dimer(forward: str, reverse: str) -> int:
    """Longest 3'-end complementarity between the two primers.

    When one primer's 3' end pairs with the other's, the pair extends each
    other into a short product that then amplifies far better than the
    template does — the reaction fills with primer-dimer and the band of
    interest never appears.
    """
    best = 0
    for size in range(4, 9):
        if len(forward) < size or len(reverse) < size:
            break
        if reverse_complement(forward[-size:]) in reverse:
            best = max(best, size)
        if reverse_complement(reverse[-size:]) in forward:
            best = max(best, size)
    return best


def _primer_warnings(anneals: str, *, label: str) -> list[str]:
    """What is wrong with one annealing region, said in reaction terms."""
    notes: list[str] = []
    gc = gc_content(anneals)

    if gc < 40:
        notes.append(
            f"The {label} primer's annealing region is {gc:.0f}% GC. Below 40% "
            f"it binds weakly; extend it or shift it into a richer stretch."
        )
    elif gc > 65:
        notes.append(
            f"The {label} primer's annealing region is {gc:.0f}% GC. Above 65% "
            f"it is prone to mispriming and to secondary structure."
        )

    if not any(base in "GC" for base in anneals[-2:]):
        notes.append(
            f"The {label} primer ends in A or T. A G or C at the 3' end holds "
            f"the extending end down; without one, priming is less efficient."
        )

    if sum(base in "GC" for base in anneals[-5:]) >= 4:
        notes.append(
            f"The {label} primer has four or more G/C in its last five bases. "
            f"That 3' end binds partial matches elsewhere as readily as the "
            f"intended site."
        )

    run = longest_homopolymer(anneals)
    if run >= 5:
        notes.append(
            f"The {label} primer contains a run of {run} identical bases, "
            f"which slips during synthesis and during extension."
        )

    if _self_complementarity(anneals) >= 5:
        notes.append(
            f"The {label} primer's 3' end is complementary to its own "
            f"sequence, so it will fold back on itself instead of priming."
        )

    return notes


def _pick_anneal(
    template: str, *, start: int, direction: int,
) -> tuple[str, int, int]:
    """Choose the annealing region from one end of the target.

    Walks lengths from MIN_ANNEAL upward and keeps the one whose Tm is
    nearest TARGET_TM, preferring a G/C at the 3' end where the choice is
    otherwise close. Returns the region as it binds, plus its footprint on
    the top strand.
    """
    best: tuple[float, str, int, int] | None = None

    for length in range(MIN_ANNEAL, MAX_ANNEAL + 1):
        if direction == 1:
            lo, hi = start, start + length
            if hi > len(template):
                break
            region = template[lo:hi]
        else:
            lo, hi = start - length, start
            if lo < 0:
                break
            # A reverse primer is written 5'→3' in its own direction, so it
            # is the reverse complement of the top strand it covers.
            region = reverse_complement(template[lo:hi])

        penalty = abs(_tm(region) - TARGET_TM)
        # A 3' G/C is worth about a degree of Tm accuracy: it buys priming
        # efficiency that Tm alone does not express.
        if not any(base in "GC" for base in region[-2:]):
            penalty += 1.0

        if best is None or penalty < best[0]:
            best = (penalty, region, lo, hi)

    if best is None:
        raise SequenceError(
            f"There is not enough template to place a primer: at least "
            f"{MIN_ANNEAL} bases are needed at each end of the region you "
            f"want to amplify."
        )

    _, region, lo, hi = best
    return region, lo, hi


def _site_of(enzyme: str) -> str:
    """The recognition sequence, or a message naming what is available."""
    entry = ALL_ENZYMES.get(enzyme)
    if entry is None:
        raise SequenceError(
            f"{enzyme} is not an enzyme this tool knows. Choose one from the "
            f"restriction enzyme list."
        )
    return str(entry["recognition"])


def _resolve_custom_primer(
    sequence: str,
    *,
    template: str,
    boundary: int,
    direction: int,
    allow_tail: bool,
) -> tuple[str, str, int, int]:
    """Split an entered oligo into its 5′ addition and 3′ template match."""
    primer = clean_dna(sequence)
    if not primer:
        raise SequenceError("Enter both primer sequences before revalidating them.")
    if len(primer) > 300:
        raise SequenceError("A primer cannot exceed 300 nucleotides.")

    available = len(template) - boundary if direction == 1 else boundary
    longest = min(len(primer), available)
    match: tuple[str, str, int, int] | None = None
    for length in range(longest, MIN_ANNEAL - 1, -1):
        if direction == 1:
            start, end = boundary, boundary + length
            expected = template[start:end]
        else:
            start, end = boundary - length, boundary
            expected = reverse_complement(template[start:end])
        if primer.endswith(expected):
            match = primer[:-length], primer[-length:], start, end
            break

    label = "forward" if direction == 1 else "reverse"
    if match is None:
        raise SequenceError(
            f"The {label} primer needs an exact 3′ match of at least "
            f"{MIN_ANNEAL} bases at the requested target boundary."
        )
    tail, anneals, start, end = match
    if not allow_tail and tail:
        raise SequenceError(
            f"The conventional {label} primer contains {len(tail)} unpaired "
            "5′ base(s). Remove the addition or select cloning PCR."
        )
    return tail, anneals, start, end


def design_pcr(
    template: str,
    *,
    target_start: int = 0,
    target_end: int | None = None,
    left_enzyme: str | None = None,
    right_enzyme: str | None = None,
    clamp: int = DEFAULT_CLAMP,
    keep_frame: bool = False,
    start_codon_mode: str = "use_site",
    name: str = "product",
    forward_primer: str | None = None,
    reverse_primer: str | None = None,
) -> PcrResult:
    """Design a PCR, with or without cloning tails, and simulate it.

    With neither enzyme named this is conventional PCR: two primers that
    copy `template[target_start:target_end]` and give back exactly that.
    Name both enzymes and each primer gains a tail, so the product carries
    the sites and can be cut into an insert.

    Args:
        target_start / target_end: the region to amplify, half-open and
            0-based. `target_end` defaults to the end of the template.
        left_enzyme / right_enzyme: as they will sit on the *product*, left
            first. Both or neither.
        clamp: bases placed outside each site so the enzyme can cut. Fewer
            than the default is accepted and warned about, because the digest
            is what turns the product into an insert.
        keep_frame: pad the left tail so the amplified region stays in frame
            with the vector's reading frame downstream of the site. Only
            meaningful for an N-terminal fusion.
        start_codon_mode: when the left restriction site supplies an ``ATG``
            and the selected region also begins with ``ATG``, ``use_site``
            omits the template codon so the protein starts with one
            methionine. ``keep_both`` preserves both codons deliberately.

    Returns:
        A :class:`PcrResult`. Check `is_clean` before ordering: an enzyme
        that cuts inside the region is reported as a problem rather than
        raised, so the interface can show what is wrong with the choice.
    """
    working = clean_dna(template)
    if not working:
        raise SequenceError("Enter a template sequence to amplify.")

    stop = len(working) if target_end is None else target_end
    if not 0 <= target_start < stop <= len(working):
        raise SequenceError(
            f"The region to amplify must lie inside the template: "
            f"{target_start + 1}–{stop} was asked for on {len(working)} bases."
        )

    if (left_enzyme is None) != (right_enzyme is None):
        raise SequenceError(
            "Name an enzyme for both ends or for neither. One site alone "
            "cannot open both ends of the product."
        )

    if not 0 <= clamp <= len(_CLAMP_BASES):
        raise SequenceError(
            f"Clamp length must be between 0 and {len(_CLAMP_BASES)} bases."
        )

    if left_enzyme is not None and left_enzyme == right_enzyme:
        raise SequenceError(
            "The two enzymes must differ — using one enzyme at both ends "
            "does not force the insert's orientation."
        )

    if start_codon_mode not in START_CODON_MODES:
        raise SequenceError(
            "Start-codon mode must be 'use_site' or 'keep_both'."
        )

    if (forward_primer is None) != (reverse_primer is None):
        raise SequenceError("Enter both custom primer sequences, or neither.")

    problems: list[str] = []
    warnings: list[str] = []

    # NdeI's CATATG is the common case: the recognition site already contains
    # the translation start. The recommended primer therefore anneals after
    # a template-leading ATG. Keeping that codon as well changes the
    # expressed protein to Met-Met, so it is an explicit opt-in rather than
    # the default.
    effective_start = target_start
    left_site = _site_of(left_enzyme) if left_enzyme is not None else ""
    supplies_atg = (
        supplies_start_codon(left_enzyme) if left_enzyme is not None else False
    )
    template_supplies_atg = working[target_start:stop].startswith("ATG")
    if supplies_atg and template_supplies_atg and start_codon_mode == "use_site":
        effective_start += 3
        warnings.append(
            f"The template's first ATG was omitted because {left_enzyme}'s "
            "restriction site supplies the start codon. The expressed "
            "protein therefore begins with one methionine."
        )

    custom_primers = forward_primer is not None and reverse_primer is not None
    if custom_primers:
        fwd_tail, fwd_anneal, fwd_lo, fwd_hi = _resolve_custom_primer(
            forward_primer,
            template=working,
            boundary=effective_start,
            direction=1,
            allow_tail=left_enzyme is not None,
        )
        rev_tail, rev_anneal, rev_lo, rev_hi = _resolve_custom_primer(
            reverse_primer,
            template=working,
            boundary=stop,
            direction=-1,
            allow_tail=right_enzyme is not None,
        )
    else:
        fwd_anneal, fwd_lo, fwd_hi = _pick_anneal(
            working, start=effective_start, direction=1,
        )
        rev_anneal, rev_lo, rev_hi = _pick_anneal(
            working, start=stop, direction=-1,
        )
        fwd_tail = rev_tail = ""

    if fwd_hi > rev_lo:
        raise SequenceError(
            "The region to amplify is too short for two primers that do not "
            f"overlap: {stop - target_start} bases were asked for, and "
            f"{2 * MIN_ANNEAL} is the minimum."
        )

    region = working[effective_start:stop]

    if left_enzyme is not None and right_enzyme is not None:
        right_site = _site_of(right_enzyme)
        if custom_primers:
            reverse_site = reverse_complement(right_site)
            if left_site not in fwd_tail:
                raise SequenceError(
                    f"The forward primer's 5′ addition does not contain the "
                    f"selected {left_enzyme} site ({left_site})."
                )
            if reverse_site not in rev_tail:
                raise SequenceError(
                    f"The reverse primer's 5′ addition does not contain the "
                    f"selected {right_enzyme} site in primer orientation "
                    f"({reverse_site})."
                )
            left_site_at = fwd_tail.index(left_site)
            reverse_site_at = rev_tail.index(reverse_site)
            if left_site_at < DEFAULT_CLAMP or reverse_site_at < DEFAULT_CLAMP:
                warnings.append(
                    "One or both edited primers provide fewer than six bases "
                    "outside the restriction site; confirm the enzyme supplier's "
                    "end-cleavage requirement before ordering."
                )
        else:
            clamp_bases = _clamp_of(clamp)
            fwd_tail = clamp_bases + left_site
            rev_tail = clamp_bases + reverse_complement(right_site)
            left_site_at = clamp
            reverse_site_at = clamp

        # Where the digest will cut, so frame can be measured on the insert
        # rather than on the product — they differ by the stub that falls off.
        left_cut_offset = clamp + int(ALL_ENZYMES[left_enzyme]["cut_top"])  # type: ignore[arg-type]
        if supplies_atg and template_supplies_atg and start_codon_mode == "keep_both":
            warnings.append(
                f"{left_enzyme}'s site carries its own ATG and the region "
                f"already begins with one, so the protein starts Met-Met. "
                f"Drop the region's first codon if the extra residue is "
                f"unwanted — with {left_enzyme} the site supplies the start."
            )

        if keep_frame and not supplies_atg and not custom_primers:
            # The vector supplies the ATG upstream, so the region must sit a
            # whole number of codons from the insert's own 5' end. Anchoring
            # on the product's end instead would be out by the stub the
            # digest removes.
            padding = (len(fwd_tail) - left_cut_offset) % 3
            padding = (3 - padding) % 3
            fwd_tail += "A" * padding
            if padding:
                warnings.append(
                    f"{padding} base(s) were added between the {left_enzyme} "
                    f"site and the insert to keep the vector's reading frame. "
                    f"They are translated as part of the junction."
                )

        if not custom_primers and clamp < DEFAULT_CLAMP:
            warnings.append(
                f"Only {clamp} clamp base(s) sit outside each site. Enzymes cut "
                f"poorly at a fragment's end; {DEFAULT_CLAMP} is the usual "
                f"minimum, and a digest that under-cuts looks like a ligation "
                f"that failed."
            )

    forward_seq = fwd_tail + fwd_anneal
    reverse_seq = rev_tail + rev_anneal

    # ── The product ─────────────────────────────────────────────────────
    # Top strand: the forward primer as written, the template between the
    # primers' inner edges, then the reverse primer's complement. Built from
    # the primers rather than assembled from the design, so that what is
    # reported is what those two oligos would actually produce.
    product = fwd_tail + region + reverse_complement(rev_tail)

    # Validate the molecule that will actually be digested, not just the
    # template region. Clamp/site and site/insert boundaries can create an
    # overlapping recognition sequence even when the original gene is clean.
    if left_enzyme is not None and right_enzyme is not None:
        right_site = _site_of(right_enzyme)
        expected = {
            left_enzyme: [left_site_at],
            right_enzyme: [
                len(fwd_tail) + len(region) + len(rev_tail)
                - reverse_site_at - len(right_site)
            ],
        }
        for enzyme in (left_enzyme, right_enzyme):
            found = find_sites(product, enzyme, circular=False)
            if found != expected[enzyme]:
                wanted = expected[enzyme][0] + 1
                observed = ", ".join(str(position + 1) for position in found) or "none"
                problems.append(
                    f"The complete PCR product has {enzyme} sites at {observed}; "
                    f"only the intended tail site at position {wanted} is "
                    "allowed. The enzyme cuts inside the amplified region, "
                    "or a clamp/site junction created an overlapping site. "
                    "Choose another enzyme or optimise the "
                    "insert before ordering the primers."
                )

    fwd_warnings = _primer_warnings(fwd_anneal, label="forward")
    rev_warnings = _primer_warnings(rev_anneal, label="reverse")

    forward = PcrPrimer(
        name=f"{name}_F", sequence=forward_seq, tail=fwd_tail, anneals=fwd_anneal,
        direction=1, start=fwd_lo, end=fwd_hi,
        tm=_tm(fwd_anneal), tm_full=_tm(forward_seq),
        gc=round(gc_content(forward_seq), 1), enzyme=left_enzyme,
        restriction_site=left_site,
        warnings=tuple(fwd_warnings),
    )
    reverse = PcrPrimer(
        name=f"{name}_R", sequence=reverse_seq, tail=rev_tail, anneals=rev_anneal,
        direction=-1, start=rev_lo, end=rev_hi,
        tm=_tm(rev_anneal), tm_full=_tm(reverse_seq),
        gc=round(gc_content(reverse_seq), 1), enzyme=right_enzyme,
        restriction_site=(reverse_complement(_site_of(right_enzyme)) if right_enzyme else ""),
        warnings=tuple(rev_warnings),
    )

    warnings.extend(fwd_warnings)
    warnings.extend(rev_warnings)

    # Ta from the annealing portions, which is what binds in cycle one.
    annealing_temperature = round(min(min(forward.tm, reverse.tm) - TA_OFFSET, MAX_TA), 1)

    if abs(forward.tm - reverse.tm) > MAX_TM_DIFFERENCE:
        warnings.append(
            f"The two annealing regions differ by "
            f"{abs(forward.tm - reverse.tm):.1f} °C ({forward.tm} and "
            f"{reverse.tm}). One primer will be near its limit at any single "
            f"annealing temperature; a gradient will find the workable band."
        )

    dimer = _cross_dimer(forward.sequence, reverse.sequence)
    if dimer >= 5:
        warnings.append(
            f"The two primers are complementary over {dimer} bases at a 3' "
            f"end. They will extend one another into primer-dimer, which "
            f"amplifies faster than the template does."
        )

    result = PcrResult(
        forward=forward, reverse=reverse,
        product=product, amplified_region=region,
        template_start=effective_start, template_end=stop,
        annealing_temperature=annealing_temperature,
        left_enzyme=left_enzyme, right_enzyme=right_enzyme,
        problems=problems, warnings=warnings,
    )

    # ── The digest ──────────────────────────────────────────────────────
    if left_enzyme is not None and right_enzyme is not None and not problems:
        result.digest = digest_linear(
            product, left_enzyme=left_enzyme, right_enzyme=right_enzyme,
        )
        if keep_frame:
            # Anchored on whichever ATG actually starts translation. With an
            # enzyme like NdeI that is the one inside the site, and the
            # region follows it in frame; with any other enzyme the vector
            # supplies it upstream and the frame begins where the region
            # does. `clone()` is what checks this against the real vector —
            # here there is no vector to check against.
            region_at = len(fwd_tail) - result.digest.trimmed_left
            if supplies_atg:
                retained, _ = left_remainders(left_enzyme)
                result.insert_orf_start = len(retained) - 3
            else:
                result.insert_orf_start = region_at
                if custom_primers and region_at % 3:
                    result.warnings.append(
                        "The edited forward primer places the target "
                        f"{region_at % 3} nucleotide(s) out of the requested "
                        "reading frame after digestion."
                    )
            _check_frame(result, region)

    return result


def _check_frame(result: PcrResult, region: str) -> None:
    """Warn when an in-frame request cannot be honoured by the region itself.

    Measured on the amplified region, not on the insert. The insert carries a
    single-stranded overhang at each end that is not part of any codon, so its
    length is almost never a multiple of three — testing that instead would
    fire on every well-formed design, and a warning that always fires is one
    the reader learns to skip past.
    """
    if len(region) % 3:
        result.warnings.append(
            f"The amplified region is {len(region)} bases, which is not a "
            f"whole number of codons. Anything fused after it — a C-terminal "
            f"tag on the vector, for instance — will be read out of frame."
        )

    protein = translate(region[: len(region) - len(region) % 3])
    if "*" in protein[:-1]:
        result.warnings.append(
            "The amplified region contains an in-frame stop codon before its "
            "end, so a C-terminal tag on the vector will not be translated."
        )
