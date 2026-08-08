"""In-silico cloning — what plasmid you actually end up with.

Everything upstream of this module tells you what to order. This one tells
you what you get: cut the vector with the two enzymes, drop the insert into
the gap, and hand back the recombinant plasmid — sequence, junctions, and
the vector's own annotations moved to their new coordinates.

It is also where the checks live that a design cannot make on its own,
because they depend on the vector:

* **Site count.** The strategy needs exactly one site per enzyme in the
  vector. A second NdeI site somewhere in the backbone means the digest
  produces three fragments and the ligation is a lottery. This is the most
  common way a perfectly good design fails at the bench.
* **Compatible ends.** The insert's sticky ends must match the backbone's,
  in the right orientation. With two different enzymes the orientation is
  forced, which is the whole reason for using a pair.
* **Reading frame.** If the vector supplies a promoter and the insert an
  ORF, the junction has to keep the frame.

**Coordinates.** 0-based, half-open, top-strand, matching the parser and the
rest of the engine. A circular vector is indexed from its own position 0;
sites that straddle that point are found by searching the sequence doubled,
which is what circular means in practice.

**Ends.** Described the same way as everywhere else in the engine: an end is
a pair (sequence in top-strand sense, which strand carries it). Two ends
ligate when their top-sense sequences are equal and both are the same
polarity — the top strand of one anneals to the bottom strand of the other.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from gsynth_engine.constants import RESTRICTION_ENZYMES, STOP_CODONS
from gsynth_engine.constants import overhang as enzyme_overhang
from gsynth_engine.sequence import (
    SequenceError,
    clean_dna,
    gc_content,
    reverse_complement,
    validate_dna,
)

if TYPE_CHECKING:                       # vectors imports us, so keep it lazy
    from gsynth_engine.vectors import VectorSpec

#: Genetic code, for the frame check. Only what this module needs.
_CODONS: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C",
    "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S",
    "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G", "TAA": "*", "TAG": "*", "TGA": "*",
}


def translate(sequence: str) -> str:
    """Translate a coding sequence from base 0. Partial final codon ignored."""
    seq = clean_dna(sequence)
    return "".join(
        _CODONS.get(seq[i : i + 3], "X") for i in range(0, len(seq) - 2, 3)
    )


# ── Ends ────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class End:
    """One end of a double-stranded fragment.

    `sequence` is the single-stranded overhang read in top-strand sense, so
    two ends anneal when their sequences are equal — regardless of which
    strand physically carries each one. `strand` says which does, because
    that is what decides whether the join is possible at all: a 5' overhang
    cannot ligate to a 3' overhang of the same sequence.
    """

    sequence: str
    strand: str        #: "top", "bottom" or "blunt"
    side: str = "left"  #: which end of the fragment this is

    @property
    def kind(self) -> str:
        """The polarity a catalogue would quote: 5', 3' or blunt.

        Which strand carries the overhang does not settle this on its own,
        because the two strands run in opposite directions. A protruding top
        strand is a 5' overhang at the fragment's left end and a 3' overhang
        at its right end; for the bottom strand it is the other way round.
        Deriving polarity from the strand alone reports NdeI — a textbook 5'
        cutter — as leaving a 3' overhang on the vector side of the junction.
        """
        if self.strand == "blunt":
            return "blunt"
        protrudes_at_its_own_five_prime = (
            (self.side == "left" and self.strand == "top")
            or (self.side == "right" and self.strand == "bottom")
        )
        return "5'" if protrudes_at_its_own_five_prime else "3'"

    def anneals_to(self, other: End) -> bool:
        """True when these two ends can be ligated together.

        One end's overhang is on the strand that runs into the join; the
        partner's is on the other strand. Two 5' overhangs meet as top-then-
        bottom, two 3' overhangs as bottom-then-top; a 5' facing a 3' leaves
        both strands unpaired at the seam and cannot close.
        """
        if self.strand == "blunt" or other.strand == "blunt":
            return self.strand == other.strand
        return self.sequence == other.sequence and self.strand != other.strand


# ── Digestion ───────────────────────────────────────────────────────────────


def find_sites(sequence: str, enzyme: str, *, circular: bool = True) -> list[int]:
    """Start positions of every recognition site, on either strand.

    A palindromic site is found once; a non-palindromic one is reported at
    the position where its top-strand match begins, which is where the cut
    offsets are measured from. On a circular sequence, sites spanning
    position 0 are found too — the doubled search is what "circular" means
    once you have to actually look for something.
    """
    if enzyme not in RESTRICTION_ENZYMES:
        raise SequenceError(f"Unknown enzyme: {enzyme}.")

    seq = clean_dna(sequence)
    site: str = RESTRICTION_ENZYMES[enzyme]["recognition"]  # type: ignore[index]
    patterns = {site, reverse_complement(site)}

    haystack = seq + seq[: len(site) - 1] if circular and len(seq) >= len(site) else seq
    positions: set[int] = set()
    for pattern in patterns:
        start = haystack.find(pattern)
        while start != -1:
            if start < len(seq):
                positions.add(start)
            start = haystack.find(pattern, start + 1)
    return sorted(positions)


@dataclass(frozen=True)
class Backbone:
    """The vector fragment that keeps the origin and the marker.

    `top` is the fragment's top strand, 5'→3', starting at the cut that will
    receive the insert's left end.
    """

    top: str
    left_end: End
    right_end: End
    #: Where `top` starts in the original vector, so annotations can follow.
    vector_start: int
    #: The stretch of vector removed by the double digest.
    removed_start: int
    removed_end: int
    removed_length: int
    circular_source: bool
    #: True when the vector had to be flipped for the insert's left enzyme to
    #: precede its right one — the cassette reads on the minus strand of the
    #: numbering the user supplied.
    reversed_insert: bool = False

    @property
    def length(self) -> int:
        return len(self.top)


def _cut_positions(enzyme: str, site_start: int) -> tuple[int, int]:
    """Top and bottom cut positions, in top-strand coordinates."""
    info = RESTRICTION_ENZYMES[enzyme]
    return (
        site_start + int(info["cut_top"]),      # type: ignore[arg-type]
        site_start + int(info["cut_bottom"]),   # type: ignore[arg-type]
    )


def _end_at(sequence: str, top_cut: int, bottom_cut: int, *, side: str) -> End:
    """The end a cut leaves on one side of itself.

    `side` is "downstream" for the fragment that begins at the cut and
    "upstream" for the one that ends there.
    """
    # The downstream fragment begins at the cut, so this is its left end;
    # the upstream one ends there, so it is its right end.
    fragment_side = "left" if side == "downstream" else "right"

    if top_cut == bottom_cut:
        return End("", "blunt", fragment_side)

    lo, hi = min(top_cut, bottom_cut), max(top_cut, bottom_cut)
    # Read base by base: on a circular vector the overhang can straddle
    # position 0, where a plain slice would silently return nothing.
    overhang = "".join(sequence[i % len(sequence)] for i in range(lo, hi))

    # A 5' overhang (top cut before bottom cut) sits on the top strand of the
    # downstream fragment and on the bottom strand of the upstream one.
    if top_cut < bottom_cut:
        return End(overhang, "top" if side == "downstream" else "bottom", fragment_side)
    return End(overhang, "bottom" if side == "downstream" else "top", fragment_side)


def linearise(
    vector: str,
    *,
    left_enzyme: str,
    right_enzyme: str,
    circular: bool = True,
) -> Backbone:
    """Cut a vector with two enzymes and return the backbone.

    The insert goes in where the removed stretch was, so the backbone's left
    end is the one the right-hand enzyme leaves and vice versa — the piece
    between the two sites is what comes out.

    **Orientation.** The two enzymes are named as they sit on the *insert*,
    which says nothing about how they sit in the vector's own numbering. In
    pET-21a(+) the expression cassette reads on the minus strand: NdeI is at
    position 236 and XhoI at 157, so the left-hand enzyme cuts *after* the
    right-hand one. Assuming otherwise picks the wrong arc of the circle and
    keeps the 80 bp cloning stuffer while discarding the origin, the marker
    and everything else — a plasmid that is arithmetically consistent and
    biologically nonsense. When that happens the vector is flipped, and the
    result records that the insert lands on the minus strand of the numbering
    the user supplied.

    Raises:
        SequenceError: when either enzyme does not cut exactly once, with a
            message naming the count. Nothing downstream can rescue a vector
            that cuts twice, so it fails here rather than producing a plan
            that will not work at the bench.
    """
    seq = validate_dna(vector, field="vector")
    if left_enzyme == right_enzyme:
        raise SequenceError(
            "The two enzymes must differ — with one enzyme the insert could "
            "go in either orientation."
        )

    if not circular:
        raise SequenceError(
            "The vector must be circular. Cutting a linear vector twice leaves "
            "the backbone in two separate pieces, which cannot receive an "
            "insert as one molecule."
        )

    for enzyme in (left_enzyme, right_enzyme):
        sites = find_sites(seq, enzyme, circular=circular)
        if len(sites) != 1:
            where = "does not cut this vector" if not sites else (
                f"cuts it {len(sites)} times (positions "
                + ", ".join(str(p + 1) for p in sites) + ")"
            )
            raise SequenceError(
                f"{enzyme} {where}. Cloning with this pair needs exactly one "
                f"site for each enzyme; otherwise the digest produces extra "
                f"fragments and the ligation cannot be directed."
            )

    length = len(seq)

    def arc(sequence: str) -> int:
        """How much of the circle the insert would replace, in this sense."""
        left = _cut_positions(left_enzyme, find_sites(sequence, left_enzyme)[0])[0]
        right = _cut_positions(right_enzyme, find_sites(sequence, right_enzyme)[0])[0]
        return (right - left) % length

    flipped = reverse_complement(seq)
    # Both arcs are geometrically valid cuts; only one leaves a backbone that
    # still carries the origin and the marker. That is the larger one.
    reversed_insert = arc(flipped) < arc(seq)
    working = flipped if reversed_insert else seq

    left_site = find_sites(working, left_enzyme, circular=True)[0]
    right_site = find_sites(working, right_enzyme, circular=True)[0]
    left_top, left_bottom = _cut_positions(left_enzyme, left_site)
    right_top, right_bottom = _cut_positions(right_enzyme, right_site)

    # The backbone runs from the right-hand enzyme's cut, around the origin,
    # back to the left-hand enzyme's cut.
    start = right_top % length
    stop = left_top % length
    doubled = working + working
    span = (stop - start) % length or length
    top = doubled[start : start + span]

    return Backbone(
        top=top,
        left_end=_end_at(working, right_top, right_bottom, side="downstream"),
        right_end=_end_at(working, left_top, left_bottom, side="upstream"),
        vector_start=start,
        removed_start=stop,
        removed_end=start,
        removed_length=(start - stop) % length,
        reversed_insert=reversed_insert,
        circular_source=circular,
    )


# ── Ligation ────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Junction:
    """One of the two seams between vector and insert."""

    name: str
    enzyme: str
    overhang: str
    kind: str                 #: 5', 3' or blunt
    position: int             #: in the recombinant plasmid, 0-based
    context: str              #: sequence either side, for eyeballing
    site_regenerated: bool    #: whether the enzyme can cut here again


@dataclass(frozen=True)
class TagOutcome:
    """Whether one of the vector's tags ends up on the protein.

    Answered by looking at the translated product rather than by a rule
    about the enzyme, because the answer depends on the reading frame and on
    whether the insert stops — and because a lab's own copy of a vector may
    not place its tags where the catalogue says.
    """

    name: str
    end: str            #: "N" or "C", where the vector places it
    present: bool
    position: int | None = None   #: residue index in the protein, 1-based
    note: str = ""


@dataclass
class CloningResult:
    """The plasmid you get, and why you should believe it."""

    plasmid: str                       #: recombinant top strand, circular
    name: str
    insert_start: int                  #: in the plasmid, 0-based
    insert_end: int
    backbone_length: int
    removed_length: int
    left_enzyme: str
    right_enzyme: str
    junctions: list[Junction] = field(default_factory=list)
    annotations: list[dict] = field(default_factory=list)
    protein: str = ""
    #: True when the insert reads on the minus strand of the vector's own
    #: numbering, which is how every pET expression cassette is arranged.
    reversed_insert: bool = False
    tags: list[TagOutcome] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.plasmid)

    @property
    def gc(self) -> float:
        return round(gc_content(self.plasmid), 1)

    @property
    def insert_length(self) -> int:
        return self.insert_end - self.insert_start

    @property
    def is_clonable(self) -> bool:
        """Empty problems means these two molecules really do join."""
        return not self.problems


def _flip_annotations(annotations: list[dict], length: int) -> list[dict]:
    """Mirror features onto the reverse complement of a circular sequence."""
    flipped: list[dict] = []
    for feature in annotations:
        entry = dict(feature)
        start, end = int(feature.get("start", 0)), int(feature.get("end", 0))
        entry["start"] = length - end
        entry["end"] = length - start
        entry["direction"] = -int(feature.get("direction", 1) or 1)
        flipped.append(entry)
    return flipped


def _remap_annotations(
    annotations: list[dict], backbone: Backbone, plasmid_length: int
) -> list[dict]:
    """Move the vector's features onto the recombinant plasmid.

    Features inside the removed stretch are dropped — they are not in the
    molecule any more, and drawing them would be a lie. Features that merely
    straddle a junction are kept and flagged, because a truncated promoter is
    something the user needs to see rather than have silently deleted.

    When the vector was flipped so the insert reads left to right, the
    features are mirrored first — otherwise every one of them lands on the
    wrong side of the plasmid, which looks like a plausible map and is not.
    """
    if not backbone.circular_source:
        return []

    if backbone.reversed_insert:
        annotations = _flip_annotations(
            annotations, backbone.length + backbone.removed_length
        )

    moved: list[dict] = []
    span = backbone.length
    for feature in annotations:
        start = int(feature.get("start", 0))
        end = int(feature.get("end", 0))
        # Position relative to the backbone's own start.
        new_start = (start - backbone.vector_start) % (span + backbone.removed_length)
        new_end = new_start + (end - start)

        if new_start >= span:
            continue                       # entirely in the removed stretch
        entry = dict(feature)
        entry["start"] = new_start
        if new_end > span:
            entry["end"] = span
            entry["truncated"] = True
        else:
            entry["end"] = new_end
        entry["end"] = min(entry["end"], plasmid_length)
        moved.append(entry)
    return moved


def clone(
    vector: str,
    insert: str,
    *,
    left_enzyme: str,
    right_enzyme: str,
    circular: bool = True,
    name: str = "recombinant",
    vector_annotations: list[dict] | None = None,
    vector_spec: VectorSpec | None = None,
    insert_reverse: str | None = None,
    insert_left_end: End | None = None,
    insert_right_end: End | None = None,
    orf_start: int | None = None,
) -> CloningResult:
    """Ligate an insert into a digested vector.

    Args:
        insert: the insert's top strand, including the bases the enzymes
            leave behind. For a G-Synth design this is `SSDResult.forward`
            or `AssemblyPlan.construct_forward` — they are built to fit.
        insert_reverse: the insert's reverse strand, 5'→3'. Supply it
            whenever you have it: both sticky ends can then be read off the
            actual duplex rather than assumed from the enzyme, which is the
            only way the compatibility check catches anything.
        insert_left_end / insert_right_end: the ends, when they are already
            known. Overrides both of the above.
        vector_spec: the catalogue entry for this vector, if it is a known
            one. Its declared tags are then checked against the protein that
            comes out, which is the only way to answer whether a C-terminal
            His-tag actually appears on this particular construct.
        orf_start: where the reading frame starts inside the insert, so the
            protein can be translated and the junction frame checked.

    Returns:
        A :class:`CloningResult`. Check `is_clonable` before believing the
        plasmid: incompatible ends produce problems rather than an exception,
        so the caller can show the user what does not fit.
    """
    insert_top = validate_dna(insert, field="insert")
    backbone = linearise(
        vector, left_enzyme=left_enzyme, right_enzyme=right_enzyme, circular=circular
    )

    problems: list[str] = []
    warnings: list[str] = []

    if insert_left_end and insert_right_end:
        left, right = insert_left_end, insert_right_end
    elif insert_reverse:
        left, right = _observed_insert_ends(insert_top, insert_reverse, left_enzyme)
    else:
        # Only the forward strand was supplied, so half the geometry is
        # unobservable: at a 5' cut the left overhang is on the top strand and
        # the right one is not. Check what can be checked and say so.
        left = _expected_insert_end(left_enzyme, side="left")
        right = _expected_insert_end(right_enzyme, side="right")
        for end, enzyme, side in ((left, left_enzyme, "left"), (right, right_enzyme, "right")):
            if end.strand == "blunt":
                continue
            visible = insert_top[: len(end.sequence)] if side == "left" \
                else insert_top[-len(end.sequence):]
            if end.strand == "top" and visible != end.sequence:
                problems.append(
                    f"The insert's {side} end reads {visible} where "
                    f"{enzyme} leaves {end.sequence}. This insert was not cut "
                    f"with {enzyme}."
                )
        warnings.append(
            "Only the forward strand was supplied, so one of the two ends was "
            "assumed from the enzyme rather than read from the molecule. "
            "Pass the reverse strand to check both."
        )

    if not backbone.right_end.anneals_to(left):
        problems.append(
            f"The insert's left end ({left.kind} {left.sequence or 'blunt'}) does "
            f"not match the vector's {left_enzyme} end "
            f"({backbone.right_end.kind} {backbone.right_end.sequence or 'blunt'})."
        )
    if not backbone.left_end.anneals_to(right):
        problems.append(
            f"The insert's right end ({right.kind} {right.sequence or 'blunt'}) "
            f"does not match the vector's {right_enzyme} end "
            f"({backbone.left_end.kind} {backbone.left_end.sequence or 'blunt'})."
        )

    # The plasmid, read from the backbone's start so the vector keeps its
    # familiar layout: backbone, then insert, then round to the beginning.
    plasmid = backbone.top + insert_top
    insert_start = backbone.length
    insert_end = insert_start + len(insert_top)

    # The backbone runs from the right-hand enzyme's cut round to the
    # left-hand one's, so the seam where the insert begins is the *left*
    # enzyme's, and the seam where it ends is the right enzyme's.
    junctions = [
        _junction(
            "vector → insert", left_enzyme, backbone.right_end,
            plasmid, insert_start,
        ),
        _junction(
            "insert → vector", right_enzyme, backbone.left_end,
            plasmid, insert_end % len(plasmid),
        ),
    ]

    for enzyme in (left_enzyme, right_enzyme):
        internal = [
            p for p in find_sites(insert_top, enzyme, circular=False)
        ]
        if internal:
            warnings.append(
                f"The insert contains {len(internal)} internal {enzyme} site"
                f"{'s' if len(internal) > 1 else ''}. Merzoug assembly never "
                f"digests the insert, so the build is unaffected — but a "
                f"diagnostic digest with {enzyme} will cut inside the gene."
            )

    protein = ""
    if orf_start is not None and 0 <= orf_start < len(insert_top):
        # Read the frame in the *plasmid*, not in the insert. A C-terminal
        # fusion runs off the end of the insert and into the vector, so
        # stopping at the insert's last base would report a protein nobody
        # will ever purify.
        protein, stop_at = _translate_in_plasmid(plasmid, insert_start + orf_start)

        if protein and protein[0] != "M":
            warnings.append(
                "The reading frame does not start with a methionine — check "
                "the ATG and the left-hand enzyme."
            )

        # Distances are measured along the frame, not in plasmid coordinates:
        # the frame can run off the end of the sequence and round to the
        # start, where a raw subtraction produces a negative residue number.
        frame_start = insert_start + orf_start
        to_insert_end = insert_end - frame_start
        last_two_codons = max(0, to_insert_end - 6)

        if stop_at is None:
            warnings.append(
                "No stop codon was found in frame anywhere in the plasmid. "
                "The construct would read around the whole molecule."
            )
        else:
            reached = (stop_at - frame_start) % len(plasmid)
            if reached < last_two_codons:
                problems.append(
                    f"A stop codon appears at residue {reached // 3 + 1}, "
                    f"{to_insert_end - reached} nt before the end of the "
                    f"insert. The protein would be truncated there."
                )
            elif reached < to_insert_end:
                # The gene's own terminus. Worth saying out loud: it means a
                # C-terminal tag offered by the vector is never translated.
                warnings.append(
                    "The insert supplies its own stop codon, so nothing "
                    "downstream in the vector is translated — a C-terminal "
                    "tag on the vector would not appear on the protein."
                )
            elif reached > to_insert_end:
                extra = reached - to_insert_end
                residues = extra // 3
                warnings.append(
                    f"The reading frame runs {extra} nt past the insert into "
                    f"the vector before stopping, adding "
                    f"{residues} vector-encoded residue"
                    f"{'' if residues == 1 else 's'} to your protein."
                )

    # Residues the insert encodes, so a vector tag is only credited when it
    # sits in the stretch the vector actually contributed.
    frame_start = insert_start + (orf_start or 0)
    insert_residues = max(0, (insert_end - frame_start)) // 3
    upstream_residues = max(0, (frame_start - insert_start)) // 3 if orf_start else 0

    tags = (
        _tag_outcomes(
            protein, vector_spec,
            insert_residues=insert_residues,
            upstream_residues=upstream_residues,
        )
        if vector_spec else []
    )
    warnings.extend(_tag_warnings(tags, vector_spec, protein))

    return CloningResult(
        plasmid=plasmid,
        name=name,
        insert_start=insert_start,
        insert_end=insert_end,
        backbone_length=backbone.length,
        removed_length=backbone.removed_length,
        left_enzyme=left_enzyme,
        right_enzyme=right_enzyme,
        junctions=junctions,
        annotations=_remap_annotations(
            vector_annotations or [], backbone, len(plasmid)
        ),
        protein=protein,
        reversed_insert=backbone.reversed_insert,
        tags=tags,
        warnings=warnings,
        problems=problems,
    )


def _tag_outcomes(
    protein: str, spec: VectorSpec, *, insert_residues: int, upstream_residues: int,
) -> list[TagOutcome]:
    """Which of the vector's tags actually made it onto the protein.

    Read from the translated product, so the answer accounts for the reading
    frame, for an insert that stops early, and for a lab copy whose tags are
    not quite where the catalogue says.

    **Where the match has to be** is the load-bearing part. Searching the
    whole protein finds the *insert's* own 6×His and reports it as the
    vector's C-terminal tag — a construct that will not bind the column,
    described as one that will. A vector tag counts only if it lies in the
    stretch the vector actually contributed: after the insert for a
    C-terminal tag, before it for an N-terminal one.
    """
    outcomes: list[TagOutcome] = []
    for tag in spec.tags:
        region = (
            protein[insert_residues:] if tag.end == "C" else protein[:upstream_residues]
        )
        offset = insert_residues if tag.end == "C" else 0
        at = region.find(tag.motif) if region else -1
        outcomes.append(
            TagOutcome(
                name=tag.name,
                end=tag.end,
                present=at >= 0,
                position=offset + at + 1 if at >= 0 else None,
                note=tag.note,
            )
        )
    return outcomes


def _tag_warnings(
    outcomes: list[TagOutcome], spec: VectorSpec | None, protein: str,
) -> list[str]:
    """Say what the vector did and did not contribute, and where it collides.

    Two things routinely surprise people. A C-terminal tag that is silently
    absent because the insert carries its own stop codon or leaves the frame
    — the construct looks right and does not bind the column. And a vector
    tag that duplicates one the G-Synth cassette already added, giving a
    protein with two His-tags and two protease sites.
    """
    if spec is None:
        return []

    messages: list[str] = []
    for outcome in outcomes:
        end = "C-terminal" if outcome.end == "C" else "N-terminal"
        motif = _motif_of(spec, outcome.name)
        if not outcome.present:
            messages.append(
                f"{spec.name}'s {end} {outcome.name} is not on this protein."
                + (f" {outcome.note}" if outcome.note else "")
            )
        elif motif and protein.count(motif) > 1:
            messages.append(
                f"{outcome.name} appears more than once: the insert already "
                f"carries one and {spec.name} adds its own. Turn the cassette "
                f"option off, or use a vector that does not supply it."
            )
    return messages


def _motif_of(spec: VectorSpec, name: str) -> str:
    for tag in spec.tags:
        if tag.name == name:
            return tag.motif
    return ""


def _translate_in_plasmid(plasmid: str, start: int) -> tuple[str, int | None]:
    """Translate from `start`, around the circle, to the first in-frame stop.

    Returns the protein (stop excluded) and the plasmid position where the
    stop codon begins, or None when the frame reads the whole way round
    without one.
    """
    length = len(plasmid)
    residues: list[str] = []
    for step in range(length // 3):
        at = (start + 3 * step) % length
        codon = "".join(plasmid[(at + i) % length] for i in range(3))
        if codon in STOP_CODONS:
            return "".join(residues), at
        residues.append(_CODONS.get(codon, "X"))
    return "".join(residues), None


def _expected_insert_end(enzyme: str, *, side: str) -> End:
    """The sticky end an insert cut with this enzyme *should* present.

    The sequence comes from the enzyme rather than from the insert's top
    strand, because only one of the two ends has its overhang on that strand.
    At the right-hand end of a 5'-overhang cut the single-stranded bases live
    on the reverse oligo; reading them off the forward strand returns four
    unrelated bases that look like a plausible overhang and are not one.

    This is what the insert is *expected* to have. `_observed_insert_ends`
    reads what it actually has.
    """
    sequence, kind = enzyme_overhang(enzyme)
    if kind == "blunt":
        return End("", "blunt", side)

    # The insert is the downstream piece of the left-hand cut and the
    # upstream piece of the right-hand one.
    if side == "left":
        return End(sequence, "top" if kind == "5'" else "bottom", "left")
    return End(sequence, "bottom" if kind == "5'" else "top", "right")


def _observed_insert_ends(
    top: str, reverse: str, left_enzyme: str,
) -> tuple[End, End]:
    """Read both ends off the actual duplex the two oligos form.

    This is the only way the end check means anything: derived ends always
    agree with the enzyme they were derived from, so comparing them to the
    backbone proves nothing. Reading the molecule catches an insert that was
    built for a different pair, or truncated, or pasted in by hand.

    The stagger between the strands comes from the left-hand enzyme's own cut
    geometry, and everything else follows from the two lengths.
    """
    info = RESTRICTION_ENZYMES[left_enzyme]
    offset = int(info["cut_bottom"]) - int(info["cut_top"])  # type: ignore[arg-type]
    bottom = reverse_complement(clean_dna(reverse))

    if offset > 0:
        left = End(top[:offset], "top", "left")
    elif offset < 0:
        left = End(bottom[:-offset], "bottom", "left")
    else:
        left = End("", "blunt", "left")

    tail = (offset + len(bottom)) - len(top)
    if tail > 0:
        right = End(bottom[-tail:], "bottom", "right")
    elif tail < 0:
        right = End(top[tail:], "top", "right")
    else:
        right = End("", "blunt", "right")

    return left, right


def _junction(
    name: str, enzyme: str, end: End, plasmid: str, position: int,
) -> Junction:
    """One seam, with enough sequence either side to read it."""
    window = 12
    length = len(plasmid)
    context = "".join(
        plasmid[(position + offset) % length]
        for offset in range(-window, window)
    )
    site: str = RESTRICTION_ENZYMES[enzyme]["recognition"]  # type: ignore[index]
    return Junction(
        name=name,
        enzyme=enzyme,
        overhang=end.sequence,
        kind=end.kind,
        position=position,
        context=context,
        # A regenerated site means the construct can be cut back out, which
        # is how most people verify a clone.
        site_regenerated=site in context or reverse_complement(site) in context,
    )


def open_reading_frames(
    sequence: str, *, minimum_codons: int = 30, circular: bool = True,
) -> list[dict]:
    """Every ORF on the top strand, longest first.

    Used to sanity-check a recombinant plasmid: the insert's ORF should be
    there, at the length the design predicted. On a circular molecule the
    scan wraps — a gene cloned near the end of the sequence still has its
    stop codon, it is simply on the other side of position 0, and a linear
    scan would report the construct as having no ORF at all.
    """
    seq = clean_dna(sequence)
    if len(seq) < 6:
        return []

    # Scanning a doubled sequence lets a frame run past the end and round to
    # the start; anything beginning past the original length is a repeat.
    scan = seq + seq if circular else seq
    found: list[dict] = []

    for frame in range(3):
        start = None
        for i in range(frame, len(scan) - 2, 3):
            codon = scan[i : i + 3]
            if start is None and codon == "ATG":
                if i >= len(seq):
                    break              # past the join: the rest is a repeat
                start = i
            elif start is not None and codon in STOP_CODONS:
                codons = (i + 3 - start) // 3
                if codons >= minimum_codons:
                    found.append({
                        "start": start,
                        "end": (i + 3) % len(seq) if circular else i + 3,
                        "frame": frame,
                        "codons": codons,
                        "wraps": circular and i + 3 > len(seq),
                        "protein": translate(scan[start : i + 3]),
                    })
                start = None
    return sorted(found, key=lambda orf: orf["codons"], reverse=True)
