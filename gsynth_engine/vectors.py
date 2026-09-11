from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

from gsynth_engine.constants import ALL_ENZYMES
from gsynth_engine.sequence import clean_dna

DATA = Path(__file__).parent / "vector_data"


T7_PROMOTER = "TAATACGACTCACTATAG"
LAC_OPERATOR = "GGAATTGTGAGCGGATAACAATT"
T7_TERMINATOR = "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG"


AMPR = "CGTTGTTGCCATTGCTGCAGGCATCGTGGTGTCACG"


@dataclass(frozen=True)
class Tag:


    name: str
    motif: str
    end: str
    note: str = ""


@dataclass(frozen=True)
class VectorSpec:


    key: str
    name: str
    length: int
    resistance: str
    promoter: str
    host: str = "E. coli"
    supplier: str = ""
    summary: str = ""

    unique_sites: tuple[str, ...] = ()

    recommended_pairs: tuple[str, ...] = ()
    tags: tuple[Tag, ...] = ()

    motifs: tuple[tuple[str, str], ...] = ()
    aliases: tuple[str, ...] = ()
    reference: str = ""
    notes: tuple[str, ...] = ()

    bundled: str = ""


    supplies_translation_start: bool = True


    expression_capable: bool = True

    @property
    def has_sequence(self) -> bool:

        return bool(self.bundled)

    @property
    def tag_summary(self) -> str:
        if not self.tags:
            return "no tags"
        return " · ".join(f"{tag.end}-terminal {tag.name}" for tag in self.tags)


HIS6 = "HHHHHH"
THROMBIN_PEPTIDE = "LVPRGS"
T7_TAG_PEPTIDE = "MASMTGGQQMG"


CATALOGUE: tuple[VectorSpec, ...] = (
    VectorSpec(
        key="pET-21a",
        name="pET-21a(+)",
        length=5443,
        resistance="Ampicillin",
        promoter="T7lac",
        supplier="Novagen / Merck",
        summary="T7 expression with a C-terminal His-tag. The default: it "
                "supplies the ribosome binding site and the ATG, inside an "
                "NdeI site — the pair the G-Synth cassette is built for.",
        unique_sites=("NdeI", "XhoI", "BamHI", "EcoRI", "HindIII", "NotI",
                      "SacI", "SalI", "XbaI", "BglII", "ApaI", "MluI", "EcoRV",
                      "PstI"),
        recommended_pairs=("NdeI / XhoI", "NdeI / EcoRI", "NdeI / HindIII",
                           "BamHI / XhoI", "NdeI / NotI"),
        tags=(
            Tag("T7·Tag", T7_TAG_PEPTIDE, "N",
                "Sits just downstream of the ATG, so cloning at NdeI replaces "
                "it with your own sequence."),
            Tag("His-tag", HIS6, "C",
                "Immediately after XhoI. It appears only if the insert reads "
                "in frame through XhoI without its own stop codon."),
        ),
        motifs=(("T7 promoter", T7_PROMOTER), ("lac operator", LAC_OPERATOR),
                ("AmpR", AMPR)),
        aliases=("pET21a", "pET21a(+)", "pET-21a(+)"),
        reference="https://www.snapgene.com/plasmids/pet_and_duet_vectors_(novagen)/pET-21a(%2B)",
        bundled="pET-21a.json",
        supplies_translation_start=True,
        notes=(
            (
                "The His-tag is C-terminal, so an insert carrying its own stop "
                "codon will not be tagged. Leave the stop off, or put the tag on "
                "the insert instead."
            ),
            "Cloning at NdeI uses the vector's ATG and removes the T7·Tag.",
        ),
    ),
    VectorSpec(
        key="pET-21",
        name="pET-21(+)",
        length=5369,
        resistance="Ampicillin",
        promoter="T7lac",
        supplier="Novagen / Merck",
        summary="T7 expression with a C-terminal His-tag. The original pET-21: "
                "its cloning region runs BamHI to XhoI and it carries no "
                "ribosome binding site or start codon of its own.",
        unique_sites=("BamHI", "EcoRI", "SacI", "SalI", "HindIII", "NotI", "XhoI",
                      "BglII", "ApaI", "MluI", "EcoRV", "PstI"),
        recommended_pairs=("BamHI / XhoI", "BamHI / NotI", "EcoRI / XhoI",
                           "BamHI / HindIII", "SacI / XhoI"),
        tags=(
            Tag("His-tag", HIS6, "C",
                "Sits immediately after XhoI, so it appears only if the "
                "insert reads in frame through XhoI without its own stop."),
        ),
        motifs=(("T7 promoter", T7_PROMOTER), ("lac operator", LAC_OPERATOR),
                ("AmpR", AMPR)),
        aliases=("pET21", "pET21(+)", "pET-21(+)"),
        reference="Bundled from a verified laboratory reference file.",
        bundled="pET-21.json",
        supplies_translation_start=False,
        notes=(
            (
                "There is no NdeI site: this is pET-21(+), not pET-21a(+). The "
                "a/b/c/d variants add 74 bp carrying NdeI, NheI and the T7·Tag."
            ),
            (
                "No ribosome binding site and no ATG between the lac operator and "
                "the cloning region. The insert has to bring its own, or nothing "
                "is translated."
            ),
        ),
    ),
    VectorSpec(
        key="pET-28a",
        name="pET-28a(+)",
        length=5369,
        resistance="Kanamycin",
        promoter="T7lac",
        supplier="Novagen / Merck",
        summary="T7 expression with an N-terminal His-tag and a thrombin site, "
                "plus an optional C-terminal His-tag.",
        unique_sites=("NdeI", "XhoI", "NcoI", "BamHI", "EcoRI", "HindIII", "NotI", "SacI", "SalI"),
        recommended_pairs=("NcoI / XhoI", "NdeI / XhoI", "BamHI / XhoI", "NdeI / HindIII"),
        tags=(
            Tag("His-tag", HIS6, "N", "Whether it survives depends on where you cut."),
            Tag("thrombin site", THROMBIN_PEPTIDE, "N",
                "Removes the N-terminal tag after purification."),
            Tag("His-tag", HIS6, "C", "Only without a stop codon in the insert."),
        ),
        motifs=(("T7 promoter", T7_PROMOTER), ("lac operator", LAC_OPERATOR)),
        aliases=("pET28a", "pET28a(+)", "pET-28a(+)"),
        reference="https://www.snapgene.com/plasmids/pet_and_duet_vectors_(novagen)/pET-28a(%2B)",
        notes=(
            (
                "NcoI / XhoI replaces the vector's N-terminal His-tag and thrombin "
                "cassette with the insert cassette. NdeI or BamHI retains the "
                "upstream vector leader; review the complete fusion and avoid "
                "duplicating tags."
            ),
        ),
    ),
    VectorSpec(
        key="pET-22b",
        name="pET-22b(+)",
        length=5493,
        resistance="Ampicillin",
        promoter="T7lac",
        supplier="Novagen / Merck",
        summary="T7 expression with an N-terminal pelB leader that exports the "
                "product to the periplasm, and a C-terminal His-tag.",
        unique_sites=("NdeI", "XhoI", "NcoI", "BamHI", "EcoRI", "HindIII", "NotI", "SacI", "SalI"),
        recommended_pairs=("NcoI / XhoI", "NdeI / XhoI", "BamHI / XhoI"),
        tags=(
            Tag("pelB leader", "MKYLLPTAAAGLLLLAAQPAMA", "N",
                "Cleaved by signal peptidase on export, so it is not on the "
                "mature protein."),
            Tag("His-tag", HIS6, "C", "Only without a stop codon in the insert."),
        ),
        motifs=(("T7 promoter", T7_PROMOTER), ("lac operator", LAC_OPERATOR)),
        aliases=("pET22b", "pET22b(+)", "pET-22b(+)"),
        reference="https://www.addgene.org/vector-database/2555/",
        notes=(
            (
                "Cloning at NdeI removes the pelB leader, and with it periplasmic "
                "export. Use NcoI to keep it."
            ),
        ),
    ),
    VectorSpec(
        key="pET-32a",
        name="pET-32a(+)",
        length=5900,
        resistance="Ampicillin",
        promoter="T7lac",
        supplier="Novagen / Merck",
        summary="Thioredoxin fusion, for proteins that are otherwise insoluble.",
        unique_sites=("NdeI", "XhoI", "NcoI", "BamHI", "EcoRI", "HindIII", "SacI", "SalI"),
        recommended_pairs=("NcoI / XhoI", "NdeI / XhoI", "BamHI / XhoI"),
        tags=(
            Tag("Trx·Tag", "SDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEY", "N",
                "Thioredoxin, ~12 kDa — improves solubility, and is large "
                "enough that most people cleave it off."),
            Tag("His-tag", HIS6, "N"),
            Tag("thrombin site", THROMBIN_PEPTIDE, "N"),
        ),
        motifs=(("T7 promoter", T7_PROMOTER), ("lac operator", LAC_OPERATOR)),
        aliases=("pET32a", "pET32a(+)"),
        reference="https://www.addgene.org/vector-database/2571/",
        notes=(
            (
                "Cloning at NdeI removes the thioredoxin fusion and the N-terminal "
                "tags with it."
            ),
        ),
    ),
    VectorSpec(
        key="pGEX-4T-1",
        name="pGEX-4T-1",
        length=4969,
        resistance="Ampicillin",
        promoter="tac",
        supplier="GE / Cytiva",
        summary="GST fusion with a thrombin site. IPTG-induced from tac, so it "
                "does not need a T7 expression strain.",
        unique_sites=("BamHI", "EcoRI", "SalI", "XhoI", "NotI"),
        recommended_pairs=("BamHI / EcoRI", "BamHI / XhoI", "EcoRI / SalI"),
        tags=(
            Tag("GST", "MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHNMLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALDVVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK",
                "N", "26 kDa — always cleaved off after purification."),
            Tag("thrombin site", THROMBIN_PEPTIDE, "N"),
        ),
        motifs=(),
        aliases=("pGEX4T1", "pGEX-4T1"),
        reference="https://www.addgene.org/vector-database/2610/",
        notes=(
            (
                "No NdeI site: the G-Synth default pair does not apply. Use "
                "BamHI / EcoRI, and keep the insert in frame with GST."
            ),
        ),
    ),
    VectorSpec(
        key="pUC19",
        name="pUC19",
        length=2686,
        resistance="Ampicillin",
        promoter="lac",
        supplier="NEB",
        summary="A cloning vector, not an expression one. For holding and "
                "sequencing a construct before moving it into an expression "
                "backbone.",
        unique_sites=("EcoRI", "SacI", "KpnI", "XbaI", "SalI", "PstI", "HindIII", "BamHI", "SmaI"),
        recommended_pairs=("EcoRI / HindIII", "BamHI / EcoRI", "EcoRI / SalI"),
        tags=(),
        motifs=(),
        aliases=("pUC-19",),
        reference="https://www.addgene.org/vector-database/2871/",
        expression_capable=False,
        notes=(
            (
                "No T7 promoter and no tags — a construct cloned here is stored, "
                "not expressed."
            ),
        ),
    ),
)

DEFAULT_VECTOR = CATALOGUE[0]


@cache
def sequence_of(key: str) -> dict | None:

    spec = get(key)
    if spec is None or not spec.bundled:
        return None

    path = DATA / spec.bundled
    if not path.exists():
        return None

    record = json.loads(path.read_text())
    if len(record["sequence"]) != spec.length:
        raise ValueError(
            f"{spec.name}: the bundled sequence is "
            f"{len(record['sequence'])} bp but the entry says {spec.length}."
        )
    return record

BY_KEY: dict[str, VectorSpec] = {spec.key: spec for spec in CATALOGUE}


def get(key: str) -> VectorSpec | None:

    wanted = _normalise(key)
    for spec in CATALOGUE:
        candidates = {spec.key, spec.name, *spec.aliases}
        if any(_normalise(c) == wanted for c in candidates):
            return spec
    return None


def _normalise(text: str) -> str:
    return text.lower().replace("-", "").replace("_", "").replace(" ", "")


@dataclass
class VectorCheck:


    spec: VectorSpec
    length: int
    problems: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    found_motifs: list[str] = field(default_factory=list)
    missing_motifs: list[str] = field(default_factory=list)

    @property
    def matches(self) -> bool:

        return not self.problems


def validate(sequence: str, spec: VectorSpec, *, length_tolerance: int = 0) -> VectorCheck:

    from gsynth_engine.cloning import find_sites

    seq = clean_dna(sequence)
    check = VectorCheck(spec=spec, length=len(seq))

    if not seq:
        check.problems.append("No sequence was supplied.")
        return check

    difference = len(seq) - spec.length
    if abs(difference) > length_tolerance:
        direction = "longer" if difference > 0 else "shorter"
        check.problems.append(
            f"This sequence is {len(seq):,} bp — {abs(difference):,} bp "
            f"{direction} than {spec.name} ({spec.length:,} bp). Check you "
            f"selected the right vector."
        )

    for name, motif in spec.motifs:
        if motif in seq or _reverse_complement(motif) in seq:
            check.found_motifs.append(name)
        else:
            check.missing_motifs.append(name)
            check.problems.append(
                f"The {name} was not found. {spec.name} carries one, so this "
                f"is either a different vector or a modified copy."
            )


    essential = {
        enzyme
        for pair in spec.recommended_pairs
        for enzyme in (part.strip() for part in pair.split("/"))
    }

    for enzyme in spec.unique_sites:
        if enzyme not in ALL_ENZYMES:
            continue
        count = len(find_sites(seq, enzyme, circular=True))
        if count == 0 and enzyme in essential:
            check.problems.append(
                f"{enzyme} does not cut this sequence, but {spec.name} has a "
                f"site for it and is cloned with it. This is a different vector."
            )
        elif count == 0:
            check.notes.append(f"{enzyme} does not cut this sequence.")
        elif count > 1:
            check.notes.append(
                f"{enzyme} cuts {count} times — it cannot be used for "
                f"directional cloning here."
            )

    return check


def identify(sequence: str) -> VectorSpec | None:

    seq = clean_dna(sequence)
    if not seq:
        return None

    for spec in CATALOGUE:
        record = sequence_of(spec.key)
        if record and record["sequence"] == seq:
            return spec

    candidates = [
        spec for spec in CATALOGUE
        if len(seq) == spec.length
        and spec.motifs
        and all(
            motif in seq or _reverse_complement(motif) in seq
            for _, motif in spec.motifs
        )
    ]
    return candidates[0] if len(candidates) == 1 else None


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
