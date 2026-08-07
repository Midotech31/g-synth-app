"""The vector catalogue — what you are cloning into.

A cloning strategy is only half a design. The other half is the backbone:
what promoter drives it, what tags it already carries, whether the enzymes
you picked cut it once, and what your protein will look like once the
vector's own contributions are translated too.

**Which entries carry a sequence.** A vector ships with one only when it came
from an authoritative file — a supplier's or the lab's own SnapGene/GenBank
export — and passes `validate` against its own entry. Everything else is
metadata, and the user supplies the sequence once from their own copy.

That restraint is deliberate. A transcription error in 5 000 bases is
invisible and would silently poison every design made against it. And
working labs run modified backbones: a pET-21 with a swapped promoter is
still the vector on their bench, and an "official" copy that quietly
disagreed with it would be worse than no copy at all.

What the catalogue always does is *check*. `validate` compares a supplied
sequence against the entry — length, the enzymes that must cut exactly once,
the promoter and operator that must be there — so pasting pET-28a while
pET-21 is selected is caught immediately. The two are both 5 369 bp, so
length alone would not have caught it.

**Tags are declared as peptides, not as DNA.** Whether a vector's C-terminal
His-tag ends up on your protein depends on the reading frame and on whether
your insert carries its own stop codon — questions answered by translating
the recombinant plasmid, not by looking for a codon string. Peptide motifs
also survive the codon differences between vector generations, and between a
supplier's vector and a lab's own version of it.

References for the metadata: Novagen/Merck pET System Manual (TB055),
supplier product pages, and the Addgene vector database.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

from gsynth_engine.constants import RESTRICTION_ENZYMES
from gsynth_engine.sequence import clean_dna

#: Where bundled sequences live. A vector earns a place here only when its
#: sequence came from an authoritative file — a supplier's or a lab's own
#: SnapGene/GenBank export — and passes `validate` against its own entry.
DATA = Path(__file__).parent / "vector_data"

#: Non-coding elements whose sequence is standard across every vector that
#: carries them, so they are safe to match at the DNA level.
T7_PROMOTER = "TAATACGACTCACTATAG"
LAC_OPERATOR = "GGAATTGTGAGCGGATAACAATT"
T7_TERMINATOR = "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG"

#: A stretch from the middle of bla, taken from the verified pET-21 file.
#: pET-21(+) and pET-28a(+) are both 5 369 bp; the resistance marker is what
#: tells them apart, so identification cannot rely on length.
AMPR = "CGTTGTTGCCATTGCTGCAGGCATCGTGGTGTCACG"


@dataclass(frozen=True)
class Tag:
    """Something the vector adds to your protein, as a peptide.

    `motif` is matched against the translated product, so it is independent
    of the vector's codon usage and of whatever the lab has changed. `end`
    says where the vector places it, which decides what has to be true for
    it to appear: an N-terminal tag needs the frame to start upstream of the
    insert, a C-terminal one needs the insert *not* to stop.
    """

    name: str
    motif: str
    end: str          #: "N" or "C"
    note: str = ""


@dataclass(frozen=True)
class VectorSpec:
    """One catalogue entry. Metadata and checks — never a sequence."""

    key: str
    name: str
    length: int
    resistance: str
    promoter: str
    host: str = "E. coli"
    supplier: str = ""
    summary: str = ""
    #: Enzymes that must cut exactly once. Anything else is not this vector.
    unique_sites: tuple[str, ...] = ()
    #: Enzyme pairs the vector is designed to be cloned with, best first.
    recommended_pairs: tuple[str, ...] = ()
    tags: tuple[Tag, ...] = ()
    #: DNA motifs that must be present, as (name, sequence).
    motifs: tuple[tuple[str, str], ...] = ()
    aliases: tuple[str, ...] = ()
    reference: str = ""
    notes: tuple[str, ...] = ()
    #: File under `vector_data/`, when the sequence ships with the engine.
    bundled: str = ""
    #: Whether the vector puts a ribosome binding site and a start codon
    #: upstream of the cloning region. False means the insert must.
    supplies_translation_start: bool = True

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


#: Ordered: the first entry is what the app offers by default.
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
        reference="Bundled from a SnapGene export.",
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
        reference="Bundled from a SnapGene export.",
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
        recommended_pairs=("NdeI / XhoI", "NcoI / XhoI", "BamHI / XhoI", "NdeI / HindIII"),
        tags=(
            Tag("His-tag", HIS6, "N", "Whether it survives depends on where you cut."),
            Tag("thrombin site", THROMBIN_PEPTIDE, "N",
                "Removes the N-terminal tag after purification."),
            Tag("His-tag", HIS6, "C", "Only without a stop codon in the insert."),
        ),
        motifs=(("T7 promoter", T7_PROMOTER), ("lac operator", LAC_OPERATOR)),
        aliases=("pET28a", "pET28a(+)", "pET-28a(+)"),
        reference="https://www.addgene.org/vector-database/2565/",
        notes=(
            (
                "This vector already supplies an N-terminal His-tag and a thrombin "
                "site. Adding the same cassette to the insert gives the protein two "
                "of each."
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
    """The bundled sequence and features for a vector, or None.

    Returns the same shape the file parser produces, so a bundled vector and
    an imported one are interchangeable everywhere downstream.
    """
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
    """Look a vector up by key, name or alias — case- and dash-insensitive."""
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
    """What a supplied sequence does and does not match."""

    spec: VectorSpec
    length: int
    problems: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    found_motifs: list[str] = field(default_factory=list)
    missing_motifs: list[str] = field(default_factory=list)

    @property
    def matches(self) -> bool:
        """No problems means this really is the vector that was selected."""
        return not self.problems


def validate(sequence: str, spec: VectorSpec, *, length_tolerance: int = 0) -> VectorCheck:
    """Check a supplied sequence against a catalogue entry.

    Length and the promoter/operator motifs are the load-bearing checks: a
    sequence of the right length carrying the right promoter is the vector it
    claims to be, and one that is 74 bases short is pET-28a pretending.

    Site counts are reported as notes rather than problems. A lab that has
    knocked out one site in its own working copy has not stopped using the
    vector, but it does need to know before it picks that enzyme.
    """
    from gsynth_engine.cloning import find_sites  # local: cloning imports us

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

    # The enzymes the vector is meant to be cloned with. One of those missing
    # entirely is not a modified copy — it is a different vector. This is what
    # separates pET-21(+) from pET-28a(+), which share a length and a promoter
    # and differ by having NdeI at all.
    essential = {
        enzyme
        for pair in spec.recommended_pairs
        for enzyme in (part.strip() for part in pair.split("/"))
    }

    for enzyme in spec.unique_sites:
        if enzyme not in RESTRICTION_ENZYMES:
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
    """Work out which catalogue vector a sequence is, or None.

    An exact match against a bundled sequence is certain, so it wins. Failing
    that, length and motifs together have to single one entry out: pET-21(+)
    and pET-28a(+) are both 5 369 bp with the same promoter, and returning
    the first of the two would be a guess dressed up as an answer. When more
    than one entry fits, this returns None and lets the user choose.
    """
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
