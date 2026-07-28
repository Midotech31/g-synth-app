"""Codon optimisation — rewriting a gene for the host that will express it.

A bacteriocin gene taken from *Enterococcus* and put into *E. coli* is the
same protein and a different translation problem: the codons the donor
prefers are often the ones the host reads slowly, and a run of them stalls
the ribosome. This module rewrites the coding sequence for the host while
leaving the protein untouched.

**The protein is the invariant.** Every path through this module ends with
the same amino-acid sequence it started with, and a test asserts it for
every optimisation the suite performs. Everything else — codon choice, GC
content, which sites are avoided — is negotiable; that is not.

**Constraints, not just preferences.** Picking the host's favourite codon
everywhere produces a sequence that is easy to translate and often
impossible to clone: it will contain the very restriction sites the
construct is going to be cut with. So optimisation runs as a two-stage
process — choose codons by adaptiveness, then repair the sequence against
hard constraints by swapping synonymous codons, iterating until nothing
changes. Anything that cannot be repaired is reported rather than left for
the user to discover at the bench.

**About the usage table.** The shipped table is the Sharp & Li (1987)
relative-adaptiveness index for *E. coli*, the same index CAI was defined
against. It drives two different things, and they are not equally
sensitive: codon *choice* depends only on the ranking within each
amino-acid family, which is not controversial, while the CAI *number*
depends on the exact weights. If a CAI is going into a manuscript, build a
table from your own reference genes with `build_table` and quote that —
that is what the metric was designed for, and it removes any dependence on
a table you did not measure.

References
    Sharp P.M. & Li W.-H. (1987) Nucleic Acids Res 15:1281–1295.
    Welch M. et al. (2009) PLoS ONE 4:e7002 — expression vs codon choice.
"""
from __future__ import annotations

import math
import random
from collections import Counter
from dataclasses import dataclass, field

from gsynth_engine.cloning import find_sites, translate
from gsynth_engine.constants import RESTRICTION_ENZYMES, STOP_CODONS
from gsynth_engine.sequence import (
    SequenceError,
    clean_dna,
    gc_content,
    longest_homopolymer,
    validate_dna,
)

#: Codons grouped by the amino acid they encode. Built once from the
#: translation table so the two can never disagree.
SYNONYMS: dict[str, tuple[str, ...]] = {}
_CODON_TO_AA: dict[str, str] = {}

for _first in "TCAG":
    for _second in "TCAG":
        for _third in "TCAG":
            _codon = _first + _second + _third
            _aa = translate(_codon)
            _CODON_TO_AA[_codon] = _aa
            SYNONYMS.setdefault(_aa, ())
            SYNONYMS[_aa] += (_codon,)


@dataclass(frozen=True)
class CodonTable:
    """Relative adaptiveness per codon, plus where the numbers came from.

    `weights` are w values in the CAI sense: within each amino acid the most
    used codon is 1.0 and the others are its fraction. Provenance is carried
    with the numbers because a CAI is only meaningful against a stated
    reference set.
    """

    name: str
    source: str
    weights: dict[str, float]

    def weight(self, codon: str) -> float:
        return self.weights.get(codon.upper(), 0.0)

    def best(self, amino_acid: str) -> str:
        """The host's most-used codon for this amino acid."""
        return max(SYNONYMS[amino_acid], key=self.weight)

    def ranked(self, amino_acid: str) -> list[str]:
        """Synonymous codons, most used first."""
        return sorted(SYNONYMS[amino_acid], key=self.weight, reverse=True)

    def rare(self, threshold: float = 0.1) -> frozenset[str]:
        """Codons the host reads slowly enough to stall on a run of them."""
        return frozenset(
            codon for codon, w in self.weights.items()
            if w < threshold and len(SYNONYMS[_CODON_TO_AA[codon]]) > 1
        )


#: Sharp & Li (1987) relative adaptiveness for *E. coli*. See the module
#: docstring on when to replace this with a table of your own.
ECOLI = CodonTable(
    name="E. coli",
    source="Sharp & Li (1987), relative adaptiveness index",
    weights={
        "TTT": 0.296, "TTC": 1.000, "TTA": 0.020, "TTG": 0.020,
        "CTT": 0.042, "CTC": 0.037, "CTA": 0.007, "CTG": 1.000,
        "ATT": 0.185, "ATC": 1.000, "ATA": 0.003, "ATG": 1.000,
        "GTT": 1.000, "GTC": 0.066, "GTA": 0.495, "GTG": 0.221,
        "TCT": 1.000, "TCC": 0.744, "TCA": 0.077, "TCG": 0.017,
        "AGT": 0.085, "AGC": 0.410,
        "CCT": 0.070, "CCC": 0.012, "CCA": 0.135, "CCG": 1.000,
        "ACT": 0.965, "ACC": 1.000, "ACA": 0.076, "ACG": 0.099,
        "GCT": 1.000, "GCC": 0.122, "GCA": 0.586, "GCG": 0.424,
        "TAT": 0.239, "TAC": 1.000,
        "CAT": 0.291, "CAC": 1.000,
        "CAA": 0.124, "CAG": 1.000,
        "AAT": 0.051, "AAC": 1.000,
        "AAA": 1.000, "AAG": 0.253,
        "GAT": 0.434, "GAC": 1.000,
        "GAA": 1.000, "GAG": 0.259,
        "TGT": 0.500, "TGC": 1.000, "TGG": 1.000,
        "CGT": 1.000, "CGC": 0.356, "CGA": 0.004, "CGG": 0.004,
        "AGA": 0.004, "AGG": 0.002,
        "GGT": 1.000, "GGC": 0.724, "GGA": 0.010, "GGG": 0.019,
        # Stops are never chosen by the optimiser; TAA is the E. coli default.
        "TAA": 1.000, "TAG": 0.100, "TGA": 0.200,
    },
)

TABLES: dict[str, CodonTable] = {"E. coli": ECOLI}


def build_table(sequences: list[str], *, name: str, source: str = "") -> CodonTable:
    """Derive a usage table from reference genes.

    This is the honest way to get a CAI: measure it against genes you chose,
    usually the host's most highly expressed ones. A table built from twenty
    ribosomal protein genes says something specific; a table shipped with a
    program says only that somebody once measured something.

    Raises:
        SequenceError: if no complete codon could be read.
    """
    counts: Counter[str] = Counter()
    for entry in sequences:
        seq = clean_dna(entry)
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i : i + 3]
            if codon in _CODON_TO_AA:
                counts[codon] += 1
    if not counts:
        raise SequenceError("No complete codons were found in the reference set.")

    weights: dict[str, float] = {}
    for amino_acid, codons in SYNONYMS.items():
        top = max((counts[c] for c in codons), default=0)
        for codon in codons:
            weights[codon] = (counts[codon] / top) if top else 0.0
    return CodonTable(name=name, source=source or "user reference set", weights=weights)


def codon_adaptation_index(sequence: str, table: CodonTable = ECOLI) -> float:
    """CAI: the geometric mean of the codons' relative adaptiveness.

    Single-codon families (Met, Trp) are excluded — they carry no choice, so
    including them only pulls every gene towards 1.0. Codons with a weight of
    zero are excluded too rather than sending the mean to zero, since a
    reference set that never used a codon says nothing about how badly it
    reads.
    """
    seq = clean_dna(sequence)
    logs: list[float] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        amino_acid = _CODON_TO_AA.get(codon)
        if amino_acid is None or amino_acid == "*":
            continue
        if len(SYNONYMS[amino_acid]) < 2:
            continue
        weight = table.weight(codon)
        if weight > 0:
            logs.append(math.log(weight))
    if not logs:
        return 0.0
    return round(math.exp(sum(logs) / len(logs)), 3)


# ── Constraints ─────────────────────────────────────────────────────────────


@dataclass
class Constraints:
    """What the finished sequence must avoid.

    The defaults are what a synthesis supplier and a cloning strategy
    between them require: no site belonging to the enzymes the construct
    will be cut with, no homopolymer long enough to slip during synthesis,
    and GC in a range that both synthesises and amplifies.
    """

    #: Enzymes whose sites must not appear. Usually the cloning pair.
    avoid_enzymes: tuple[str, ...] = ()
    #: Any other motif to keep out — a supplier's blacklist, an internal site.
    avoid_motifs: tuple[str, ...] = ()
    max_homopolymer: int = 5
    gc_min: float = 30.0
    gc_max: float = 70.0
    #: Sliding window for local GC, which is what actually breaks synthesis.
    gc_window: int = 50
    #: Reject any repeated stretch at least this long.
    max_repeat: int = 15
    avoid_rare: bool = True
    rare_threshold: float = 0.1

    def motifs(self) -> tuple[str, ...]:
        """Every literal sequence to avoid, from the enzymes and the extras."""
        sites = tuple(
            str(RESTRICTION_ENZYMES[e]["recognition"])
            for e in self.avoid_enzymes
            if e in RESTRICTION_ENZYMES
        )
        return sites + tuple(m.upper() for m in self.avoid_motifs if m)


@dataclass
class OptimisationResult:
    """The rewritten gene, and what changed."""

    sequence: str
    protein: str
    table: str
    cai_before: float | None
    cai_after: float
    gc_before: float | None
    gc_after: float
    #: Motifs that were present before and are gone now.
    sites_removed: list[str] = field(default_factory=list)
    rare_codons_before: int = 0
    rare_codons_after: int = 0
    changed_codons: int = 0
    #: Breaches that stop the construct working: a forbidden restriction site
    #: makes it unclonable. Empty means the gene can be built and cut.
    problems: list[str] = field(default_factory=list)
    #: Breaches that cost quality rather than viability — a rare codon left
    #: in to satisfy the GC window, a repeat a supplier may charge more for.
    warnings: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def is_clean(self) -> bool:
        return not self.problems


def _violations(
    sequence: str, constraints: Constraints, table: CodonTable,
) -> list[tuple[int, int, str]]:
    """Every constraint breach, as (position, width, description).

    Positions and widths are in nucleotides, so the repair pass knows which
    codons to reach for — a forbidden site usually straddles two or three.
    """
    found: list[tuple[int, int, str]] = []

    for motif in constraints.motifs():
        start = sequence.find(motif)
        while start != -1:
            found.append((start, len(motif), f"contains {motif}"))
            start = sequence.find(motif, start + 1)

    run = 1
    for i in range(1, len(sequence)):
        run = run + 1 if sequence[i] == sequence[i - 1] else 1
        if run > constraints.max_homopolymer:
            found.append((i - run + 1, run, f"{run}x {sequence[i]} run"))

    window = constraints.gc_window
    if len(sequence) >= window:
        for start in range(0, len(sequence) - window + 1, 3):
            local = gc_content(sequence[start : start + window])
            if local < constraints.gc_min or local > constraints.gc_max:
                found.append((start, window, f"local GC {local:.0f}%"))

    if constraints.avoid_rare:
        rare = table.rare(constraints.rare_threshold)
        for i in range(0, len(sequence) - 2, 3):
            if sequence[i : i + 3] in rare:
                found.append((i, 3, f"rare codon {sequence[i:i + 3]}"))

    found.extend(_repeats(sequence, constraints.max_repeat))
    return found


def _repeats(sequence: str, minimum: int) -> list[tuple[int, int, str]]:
    """Stretches that occur more than once, which synthesis suppliers reject."""
    seen: set[str] = set()
    found: list[tuple[int, int, str]] = []
    for i in range(len(sequence) - minimum + 1):
        chunk = sequence[i : i + minimum]
        if chunk in seen:
            found.append((i, minimum, f"repeat of {chunk}"))
        else:
            seen.add(chunk)
    return found


def _cost(
    sequence: str,
    constraints: Constraints,
    table: CodonTable,
    *,
    window: tuple[int, int] | None = None,
    with_repeats: bool = True,
) -> float:
    """How badly a sequence breaks the constraints, as one number.

    A plain count of breaches makes a poor objective: widening a GC window
    fails every overlapping window at once, so no single codon swap ever
    reduces the count and the repair pass stalls with nothing to climb. A
    weighted cost that measures *how far* each window is out gives the swap
    something to improve, one codon at a time.

    The weights are a priority order, not a measurement. A forbidden site
    makes the construct unclonable, so it outranks everything; a rare codon
    only slows translation, so it yields to all of them.

    `window` restricts the sum to a stretch of the sequence. Swapping one
    codon changes three bases, so every term except repeats is unaffected
    outside a short neighbourhood — and scoring the whole gene for each
    candidate is what made a 3 kb optimisation take seven seconds of CPU,
    which is a denial of service anyone can trigger by pasting an operon.
    """
    if window is None:
        start, stop = 0, len(sequence)
    else:
        start, stop = max(0, window[0]), min(len(sequence), window[1])
    if stop <= start:
        return 0.0

    region = sequence[start:stop]
    total = 0.0

    for motif in constraints.motifs():
        total += 1000.0 * region.count(motif)

    run = 1
    for i in range(1, len(region)):
        run = run + 1 if region[i] == region[i - 1] else 1
        if run > constraints.max_homopolymer:
            total += 50.0

    # GC windows are enumerated on absolute positions that are multiples of
    # three, exactly as the global pass does. Stepping from the region's own
    # start instead enumerates a *different* set of windows, so a window that
    # is out of range can fall between the two and never be repaired.
    width = constraints.gc_window
    first_window = -(-start // 3) * 3
    for at in range(first_window, stop - width + 1, 3):
        local = gc_content(sequence[at : at + width])
        if local < constraints.gc_min:
            total += (constraints.gc_min - local)
        elif local > constraints.gc_max:
            total += (local - constraints.gc_max)

    if with_repeats:
        total += 20.0 * len(_repeats(sequence, constraints.max_repeat))

    if constraints.avoid_rare:
        rare = table.rare(constraints.rare_threshold)
        # Codon boundaries are absolute, so step from the frame rather than
        # from the window's own start.
        first_codon = -(-start // 3) * 3
        total += sum(
            1.0 for i in range(first_codon, stop - 2, 3)
            if sequence[i : i + 3] in rare
        )
    return total


def _fix_near(
    codons: list[str],
    position: int,
    width: int,
    protein: str,
    table: CodonTable,
    constraints: Constraints,
    *,
    local: bool = True,
) -> bool:
    """Swap one codon in or beside the breach. True when the cost dropped.

    Every codon overlapping the breach is tried, and within each the
    synonymous alternatives are tried best-first, so the sequence gives up as
    little adaptiveness as it has to. Two codons of slack either side catch
    sites that straddle the edge of the region.

    Candidates are compared on a *window* of the sequence rather than the
    whole of it. Every candidate differs from the others in the same three
    bases, so the terms outside that neighbourhood are identical and cancel —
    comparing them is the same decision at a fraction of the cost. Repeats
    are the exception, since a repeat's partner can be anywhere, so a breach
    of that kind falls back to scoring the whole sequence.
    """
    first = max(0, position // 3 - 2)
    last = min(len(codons) - 1, (position + width) // 3 + 2)

    # Wide enough that every GC window and every motif touching the edited
    # codons lies inside it.
    radius = constraints.gc_window + max(
        (len(m) for m in constraints.motifs()), default=0
    ) + constraints.max_homopolymer + 6
    span = (
        (first * 3 - radius, (last + 1) * 3 + radius) if local else None
    )
    kwargs = {"window": span, "with_repeats": not local}

    sequence = "".join(codons)
    baseline = _cost(sequence, constraints, table, **kwargs)
    best_cost, best_at, best_codon = baseline, -1, ""

    for index in range(first, last + 1):
        current = codons[index]
        for candidate in table.ranked(protein[index]):
            if candidate == current:
                continue
            codons[index] = candidate
            cost = _cost("".join(codons), constraints, table, **kwargs)
            if cost < best_cost:
                best_cost, best_at, best_codon = cost, index, candidate
        codons[index] = current

    if best_at < 0:
        return False
    codons[best_at] = best_codon
    return True


def _perturb(
    codons: list[str],
    position: int,
    width: int,
    protein: str,
    table: CodonTable,
    rng: random.Random,
) -> None:
    """Take a sideways step so a plateau does not end the repair pass."""
    first = max(0, position // 3 - 2)
    last = min(len(codons) - 1, (position + width) // 3 + 2)
    if last < first:
        return
    index = rng.randint(first, last)
    alternatives = [c for c in SYNONYMS[protein[index]] if c != codons[index]]
    if alternatives:
        codons[index] = rng.choice(alternatives)


def optimise(
    sequence: str,
    *,
    table: CodonTable = ECOLI,
    constraints: Constraints | None = None,
    is_protein: bool = False,
    keep_stop: bool = True,
    seed: int = 0,
    max_rounds: int = 40,
) -> OptimisationResult:
    """Rewrite a gene for the host, keeping the protein identical.

    Args:
        sequence: a coding sequence, or a protein when `is_protein` is set.
        table: the host's codon usage. Build your own with `build_table` when
            the CAI matters.
        constraints: what the result must avoid. Pass the cloning enzymes in
            `avoid_enzymes` — a gene carrying an internal NdeI site cannot be
            cloned NdeI/XhoI however well it translates.
        keep_stop: append the host's preferred stop codon. Turn it off for an
            insert destined for a C-terminal vector tag, where a stop would
            silently remove the tag.
        seed: choices are deterministic. The same input gives the same gene,
            which matters when a sequence has been ordered and someone wants
            to know it was this design that produced it.

    Returns:
        An :class:`OptimisationResult`. `problems` lists anything the repair
        pass could not fix — a stretch where no synonymous codon removes a
        site, for instance.

    Raises:
        SequenceError: for input that cannot be read as a gene or a protein.
    """
    constraints = constraints or Constraints()
    rng = random.Random(seed)

    if is_protein:
        protein = "".join(sequence.split()).upper()
        invalid = sorted(set(protein) - set(SYNONYMS) - {"*"})
        if invalid:
            raise SequenceError(
                "The protein contains characters that are not amino acids: "
                + ", ".join(invalid)
            )
        if not protein:
            raise SequenceError("The protein is empty.")
        original = None
    else:
        original = validate_dna(sequence, field="sequence")
        if len(original) % 3:
            raise SequenceError(
                f"The coding sequence is {len(original)} nt, which is not a "
                f"multiple of three. Trim it to whole codons first."
            )
        protein = translate(original)

    protein = protein.rstrip("*")
    if not protein:
        raise SequenceError("There is nothing to optimise: the protein is empty.")

    # ── Stage one: choose by adaptiveness ───────────────────────────────────
    codons = [table.best(aa) for aa in protein]

    # ── Stage two: repair against the constraints ───────────────────────────
    # A hill-climb on the cost above: fix the worst breach, re-measure, repeat.
    # Sideways steps break plateaus; giving up is reported, never silent.
    problems: list[str] = []
    stalled = 0
    for _ in range(max_rounds):
        breaches = _violations("".join(codons), constraints, table)
        if not breaches:
            break
        position, width, reason = breaches[0]
        # A repeat's partner can be anywhere, so that one breach is the one
        # kind a local comparison cannot judge.
        locally = not reason.startswith("repeat")
        if _fix_near(
            codons, position, width, protein, table, constraints, local=locally,
        ):
            stalled = 0
            continue
        stalled += 1
        if stalled > 4:
            problems.append(
                f"Could not remove: {reason} at position {position + 1}. No "
                f"synonymous codon in that stretch helps — the residues there "
                f"may have only one codon each."
            )
            break
        _perturb(codons, position, width, protein, table, rng)

    optimised = "".join(codons)
    if keep_stop:
        optimised += table.best("*")

    # Severity is about consequence, not about which constraint was set. A
    # site left in cannot be cloned around; a rare codon left in costs a
    # little translation speed, and calling both "problems" would train the
    # user to ignore the word.
    warnings: list[str] = []
    blocking = tuple(constraints.motifs())
    for position, _width, reason in _violations(optimised, constraints, table):
        message = f"{reason} at position {position + 1}."
        if any(motif in reason for motif in blocking):
            if message not in problems:
                problems.append(message)
        elif message not in warnings:
            warnings.append(message)

    # ── What changed ────────────────────────────────────────────────────────
    rare = table.rare(constraints.rare_threshold)
    sites_removed = [
        motif for motif in constraints.motifs()
        if original and motif in original and motif not in optimised
    ]
    if original and len(original) == len(optimised):
        changed = sum(
            1 for i in range(0, len(original), 3)
            if original[i : i + 3] != optimised[i : i + 3]
        )
    else:
        changed = len(codons)

    for enzyme in constraints.avoid_enzymes:
        if enzyme in RESTRICTION_ENZYMES and find_sites(
            optimised, enzyme, circular=False
        ):
            warnings.append(
                f"An {enzyme} site remains. The construct cannot be cut with "
                f"{enzyme} without cutting inside the gene."
            )

    return OptimisationResult(
        sequence=optimised,
        protein=protein,
        table=table.name,
        cai_before=codon_adaptation_index(original, table) if original else None,
        cai_after=codon_adaptation_index(optimised, table),
        gc_before=round(gc_content(original), 1) if original else None,
        gc_after=round(gc_content(optimised), 1),
        sites_removed=sites_removed,
        rare_codons_before=(
            sum(
                1 for i in range(0, len(original) - 2, 3)
                if original[i : i + 3] in rare
            ) if original else 0
        ),
        rare_codons_after=sum(
            1 for i in range(0, len(optimised) - 2, 3)
            if optimised[i : i + 3] in rare
        ),
        changed_codons=changed,
        problems=problems,
        warnings=warnings,
    )


def _swap(
    codons: list[str],
    index: int,
    protein: str,
    table: CodonTable,
    constraints: Constraints,
    rng: random.Random,
) -> bool:
    """Try a different synonymous codon at `index`. True when it helped.

    Candidates are tried best-first, so the sequence gives up as little
    adaptiveness as it has to. A swap is kept only if it reduces the number
    of breaches; otherwise the codon is put back, which stops the repair pass
    from wandering.
    """
    if not 0 <= index < len(codons):
        return False

    amino_acid = protein[index]
    alternatives = [c for c in table.ranked(amino_acid) if c != codons[index]]
    if not alternatives:
        return False

    before = len(_violations("".join(codons), constraints, table))
    original = codons[index]

    for candidate in alternatives:
        codons[index] = candidate
        after = len(_violations("".join(codons), constraints, table))
        if after < before:
            return True

    codons[index] = original
    # A last resort when no single swap improves the count: take a random
    # alternative anyway, so a stuck position can still move.
    if len(alternatives) > 1:
        codons[index] = rng.choice(alternatives)
        return True
    return False


def back_translate(protein: str, *, table: CodonTable = ECOLI) -> str:
    """The host's preferred codon for each residue, with no constraints.

    Kept separate from `optimise` because it is a different thing: this is
    what "reverse translate" means in a sequence editor, and it makes no
    claim about being clonable.
    """
    residues = "".join(protein.split()).upper()
    invalid = sorted(set(residues) - set(SYNONYMS))
    if invalid:
        raise SequenceError(
            "The protein contains characters that are not amino acids: "
            + ", ".join(invalid)
        )
    return "".join(table.best(aa) for aa in residues)
