"""Codon optimisation — rewriting a gene for the host that will express it.

A bacteriocin gene taken from *Enterococcus* and put into *E. coli* is the
same protein and a different translation context: synonymous-codon usage
differs across taxa, while translation also depends on tRNA supply, transcript
structure, growth state and the protein itself. This module rewrites the
coding sequence against a documented host profile while leaving the protein
untouched; it does not claim to predict expression yield.

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

**About the usage tables.** Bundled profiles are normalized from a committed
FDA HIVE-CUTs/CoCoPUTs snapshot. RefSeq genomic species aggregates are used
when present, with a disclosed GenBank fallback for organisms absent from that
RefSeq snapshot. Codon *choice* depends on the ranking within each amino-acid
family. A strict CAI, however, is defined against an explicit set of highly
expressed genes; a species-wide score is therefore labelled profile-relative,
not an expression-yield prediction. For a particular strain, tissue or cell
line, build a table from the relevant expressed genes with `build_table` and
report that reference set.

References
    Sharp P.M. & Li W.-H. (1987) Nucleic Acids Res 15:1281–1295.
    Athey J. et al. (2017) BMC Bioinformatics 18:391 — HIVE-CUTs.
    Alexaki A. et al. (2019) J Mol Biol 431:2434–2441 — CoCoPUTs.
    Ranaghan M.J. et al. (2021) BMC Biology 19:36 — algorithm inequality.
    Welch M. et al. (2009) PLoS ONE 4:e7002 — expression vs codon choice.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from gsynth_engine.cloning import translate
from gsynth_engine.constants import ALL_ENZYMES
from gsynth_engine.sequence import (
    SequenceError,
    clean_dna,
    gc_content,
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
    category: str = "Custom"
    taxon_id: int | None = None
    dataset: str = "User reference set"
    dataset_release: str = ""
    data_scope: str = "user-supplied reference genes"
    coding_sequences: int | None = None
    codon_count: int | None = None
    gc_percent: float | None = None
    source_url: str = ""

    def weight(self, codon: str) -> float:
        """Relative adaptiveness, 0–1: this codon's usage against the most
        used codon for the same amino acid, which scores 1.0.

        Not a frequency — the synonyms of one amino acid do not sum to 1.
        These are the w values CAI is the geometric mean of.

        An unknown codon returns 0.0 rather than raising, so a sequence
        containing an ambiguity code degrades the score instead of failing
        the whole optimisation.
        """
        return self.weights.get(codon.upper(), 0.0)

    def best(self, amino_acid: str) -> str:
        """The host's most-used codon for this amino acid."""
        return max(SYNONYMS[amino_acid], key=self.weight)

    def ranked(self, amino_acid: str) -> list[str]:
        """Synonymous codons, most used first."""
        return sorted(SYNONYMS[amino_acid], key=self.weight, reverse=True)

    def rare(self, threshold: float = 0.1) -> frozenset[str]:
        """Codons below a stated relative-frequency threshold.

        Low genomic frequency is not itself a measurement of elongation rate.
        """
        return frozenset(
            codon for codon, w in self.weights.items()
            if w < threshold and len(SYNONYMS[_CODON_TO_AA[codon]]) > 1
        )


_SNAPSHOT_PATH = Path(__file__).with_name("data") / "codon_usage_hive_2021.json"
_SNAPSHOT_BYTES = _SNAPSHOT_PATH.read_bytes()
CODON_DATA_SHA256 = hashlib.sha256(_SNAPSHOT_BYTES).hexdigest()
_SNAPSHOT = json.loads(_SNAPSHOT_BYTES)
CODON_DATA_VERSION = str(_SNAPSHOT["source"]["release"])
CODON_DATA_URL = str(_SNAPSHOT["source"]["url"])


def _weights_from_counts(name: str, counts: dict[str, int]) -> dict[str, float]:
    """Normalize raw codon counts within each synonymous family."""
    expected = set(_CODON_TO_AA)
    if set(counts) != expected:
        missing = ", ".join(sorted(expected - set(counts))) or "none"
        extra = ", ".join(sorted(set(counts) - expected)) or "none"
        raise ValueError(f"Invalid codon table for {name}: missing={missing}; extra={extra}")
    weights: dict[str, float] = {}
    for codons in SYNONYMS.values():
        maximum = max(counts[codon] for codon in codons)
        if maximum <= 0:
            raise ValueError(f"No observations for one synonymous family in {name}")
        for codon in codons:
            weights[codon] = counts[codon] / maximum
    return weights


def _load_bundled_tables() -> dict[str, CodonTable]:
    release = CODON_DATA_VERSION
    tables: dict[str, CodonTable] = {}
    for key, record in _SNAPSHOT["hosts"].items():
        counts = {codon: int(value) for codon, value in record["counts"].items()}
        dataset = str(record["dataset"])
        taxon_id = int(record["taxon_id"])
        coding_sequences = int(record["coding_sequences"])
        codon_count = int(record["codon_count"])
        source = (
            f"FDA HIVE-CUTs/CoCoPUTs {release} snapshot; {dataset} genomic "
            f"species aggregate; NCBI taxon {taxon_id}; "
            f"{coding_sequences:,} CDS; {codon_count:,} codons"
        )
        tables[key] = CodonTable(
            name=str(record["name"]),
            source=source,
            weights=_weights_from_counts(str(record["name"]), counts),
            category=str(record["category"]),
            taxon_id=taxon_id,
            dataset=dataset,
            dataset_release=release,
            data_scope=str(record["data_scope"]),
            coding_sequences=coding_sequences,
            codon_count=codon_count,
            gc_percent=float(record["gc_percent"]),
            source_url=CODON_DATA_URL,
        )
    return tables


DEFAULT_HOST = "ecoli"
TABLES = _load_bundled_tables()
ECOLI = TABLES[DEFAULT_HOST]


def build_table(sequences: list[str], *, name: str, source: str = "") -> CodonTable:
    """Derive a usage table from reference genes.

    A strict CAI is measured against genes chosen a priori, usually a set of
    highly expressed genes from the exact expression context. The caller must
    document that selection; G-Synth records its size and source.

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
    for codons in SYNONYMS.values():
        top = max((counts[c] for c in codons), default=0)
        for codon in codons:
            weights[codon] = (counts[codon] / top) if top else 0.0
    return CodonTable(
        name=name,
        source=source or "user reference set",
        weights=weights,
        category="Custom",
        coding_sequences=len(sequences),
        codon_count=sum(counts.values()),
    )


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
            str(ALL_ENZYMES[e]["recognition"])
            for e in self.avoid_enzymes
            if e in ALL_ENZYMES
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
    #: Breaches that cost profile fit rather than viability — a rare codon left
    #: in to satisfy the GC window, a repeat a supplier may charge more for.
    warnings: list[str] = field(default_factory=list)
    #: Protein exactly as supplied before any expression-start handling.
    input_protein: str | None = None
    #: Resolved biological role for a peptide input.
    protein_context: str | None = None
    initiator_methionine_added: bool = False
    recommended_design_is_coding: bool = False

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def is_clean(self) -> bool:
        """No *problems* — the gene can be built and cut as asked.

        Warnings are not consulted. A low-frequency codon left in to satisfy a
        GC window reduces profile fit and leaves this True; a
        restriction site that survived optimisation blocks the strategy and
        makes it False. Severity follows consequence.
        """
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
                found.append((i, 3, f"low-frequency codon {sequence[i:i + 3]}"))

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
    protein_context: str = "auto",
    keep_stop: bool = True,
    seed: int = 0,
    max_rounds: int = 40,
) -> OptimisationResult:
    """Rewrite a gene for the host, keeping the protein identical.

    Args:
        sequence: a coding sequence, or a protein when `is_protein` is set.
        protein_context: peptide-to-DNA start logic. ``auto`` treats an
            N-terminal methionine as a complete ORF and a peptide without one
            as a mature peptide. ``mature_peptide`` always preserves the
            supplied peptide exactly. ``complete_orf`` adds one initiator
            methionine only when it is absent.
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

    input_protein: str | None = None
    resolved_context: str | None = None
    initiator_methionine_added = False
    recommended_design_is_coding = False

    if is_protein:
        input_protein = "".join(sequence.split()).upper().rstrip("*")
        invalid = sorted(set(input_protein) - set(SYNONYMS))
        if invalid:
            raise SequenceError(
                "The protein contains characters that are not amino acids: "
                + ", ".join(invalid)
            )
        if not input_protein:
            raise SequenceError("The protein is empty.")

        allowed_contexts = {"auto", "mature_peptide", "complete_orf"}
        if protein_context not in allowed_contexts:
            raise SequenceError(
                "Protein context must be auto, mature_peptide or complete_orf."
            )
        resolved_context = protein_context
        if resolved_context == "auto":
            resolved_context = (
                "complete_orf"
                if input_protein.startswith("M")
                else "mature_peptide"
            )

        protein = input_protein
        if resolved_context == "complete_orf" and not protein.startswith("M"):
            protein = "M" + protein
            initiator_methionine_added = True
        recommended_design_is_coding = resolved_context == "complete_orf"
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
    # site left in cannot be cloned around; a low-frequency codon left in
    # reduces profile fit, and calling both "problems" would train the
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
        input_protein=input_protein,
        protein_context=resolved_context,
        initiator_methionine_added=initiator_methionine_added,
        recommended_design_is_coding=recommended_design_is_coding,
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
