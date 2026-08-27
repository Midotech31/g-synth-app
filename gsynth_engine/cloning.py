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

from gsynth_engine.constants import ALL_ENZYMES, STOP_CODONS
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
    if enzyme not in ALL_ENZYMES:
        raise SequenceError(f"Unknown enzyme: {enzyme}.")

    seq = clean_dna(sequence)
    site: str = ALL_ENZYMES[enzyme]["recognition"]  # type: ignore[index]
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
    info = ALL_ENZYMES[enzyme]
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


@dataclass(frozen=True)
class Digest:
    """A linear fragment cut at both ends, and what was trimmed away.

    This is what a PCR product becomes on the bench once the tails carrying
    the restriction sites have done their job: a shorter duplex whose ends
    are sticky. `top` and `bottom` are both 5'→3' in their own directions, so
    they are what a supplier would print, not one string and its reverse.
    """

    top: str
    bottom: str
    left_end: End
    right_end: End
    #: Bases lost from each end — the clamp plus the part of the site that
    #: stays behind on the discarded stub.
    trimmed_left: int
    trimmed_right: int

    @property
    def length(self) -> int:
        """Length of the top strand, overhangs included."""
        return len(self.top)


def digest_linear(
    fragment: str,
    *,
    left_enzyme: str,
    right_enzyme: str,
) -> Digest:
    """Cut a linear fragment with one enzyme near each end.

    The middle piece is kept: the two short stubs outside the cuts are what
    fall away, which on a cloning PCR product is the clamp bases plus part of
    each recognition site. What survives is the insert, carrying the sticky
    ends the enzymes leave.

    The ends are read off the cut molecule rather than looked up from the
    enzyme table, for the same reason `terminal_ends` is measured elsewhere:
    a value copied from the table agrees with the table whatever the bases
    actually spell.

    Raises:
        SequenceError: when either enzyme fails to cut, or when the left
            enzyme's site does not lie to the left of the right enzyme's —
            which means the primers or the enzymes have been swapped, and
            the "insert" that would come back is the piece meant to be
            thrown away.
    """
    working = clean_dna(fragment)

    left_sites = find_sites(working, left_enzyme, circular=False)
    right_sites = find_sites(working, right_enzyme, circular=False)

    for enzyme, sites in ((left_enzyme, left_sites), (right_enzyme, right_sites)):
        if not sites:
            raise SequenceError(
                f"{enzyme} does not cut this fragment. Check that the primer "
                f"tail carries the {enzyme} site and that the sequence is the "
                f"PCR product rather than the template."
            )

    # With the same enzyme at both ends there is one site per end; take the
    # outermost of each so the whole insert is kept.
    left_start = min(left_sites)
    right_start = max(right_sites)
    if left_enzyme == right_enzyme and len(left_sites) < 2:
        raise SequenceError(
            f"{left_enzyme} cuts this fragment only once, so it cannot open "
            f"both ends. Use a different enzyme at one end."
        )

    left_top, left_bottom = _cut_positions(left_enzyme, left_start)
    right_top, right_bottom = _cut_positions(right_enzyme, right_start)

    if min(left_top, left_bottom) >= min(right_top, right_bottom):
        raise SequenceError(
            f"The {l