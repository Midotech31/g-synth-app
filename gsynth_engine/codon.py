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

        return self.weights.get(codon.upper(), 0.0)

    def best(self, amino_acid: str) -> str:

        return max(SYNONYMS[amino_acid], key=self.weight)

    def ranked(self, amino_acid: str) -> list[str]:

        return sorted(SYNONYMS[amino_acid], key=self.weight, reverse=True)

    def rare(self, threshold: float = 0.1) -> frozenset[str]:

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


@dataclass
class Constraints:


    avoid_enzymes: tuple[str, ...] = ()

    avoid_motifs: tuple[str, ...] = ()
    max_homopolymer: int = 5
    gc_min: float = 30.0
    gc_max: float = 70.0

    gc_window: int = 50

    max_repeat: int = 15
    avoid_rare: bool = True
    rare_threshold: float = 0.1

    def motifs(self) -> tuple[str, ...]:

        sites = tuple(
            str(ALL_ENZYMES[e]["recognition"])
            for e in self.avoid_enzymes
            if e in ALL_ENZYMES
        )
        return sites + tuple(m.upper() for m in self.avoid_motifs if m)


@dataclass
class OptimisationResult:


    sequence: str
    protein: str
    table: str
    cai_before: float | None
    cai_after: float
    gc_before: float | None
    gc_after: float

    sites_removed: list[str] = field(default_factory=list)
    rare_codons_before: int = 0
    rare_codons_after: int = 0
    changed_codons: int = 0


    problems: list[str] = field(default_factory=list)


    warnings: list[str] = field(default_factory=list)

    input_protein: str | None = None

    protein_context: str | None = None
    initiator_methionine_added: bool = False
    recommended_design_is_coding: bool = False

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def is_clean(self) -> bool:

        return not self.problems


def _violations(
    sequence: str, constraints: Constraints, table: CodonTable,
) -> list[tuple[int, int, str]]:

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

    first = max(0, position // 3 - 2)
    last = min(len(codons) - 1, (position + width) // 3 + 2)


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


    codons = [table.best(aa) for aa in protein]


    problems: list[str] = []
    stalled = 0
    for _ in range(max_rounds):
        breaches = _violations("".join(codons), constraints, table)
        if not breaches:
            break
        position, width, reason = breaches[0]


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


    warnings: list[str] = []
    blocking = tuple(constraints.motifs())
    for position, _width, reason in _violations(optimised, constraints, table):
        message = f"{reason} at position {position + 1}."
        if any(motif in reason for motif in blocking):
            if message not in problems:
                problems.append(message)
        elif message not in warnings:
            warnings.append(message)


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


    if len(alternatives) > 1:
        codons[index] = rng.choice(alternatives)
        return True
    return False


def back_translate(protein: str, *, table: CodonTable = ECOLI) -> str:

    residues = "".join(protein.split()).upper()
    invalid = sorted(set(residues) - set(SYNONYMS))
    if invalid:
        raise SequenceError(
            "The protein contains characters that are not amino acids: "
            + ", ".join(invalid)
        )
    return "".join(table.best(aa) for aa in residues)
