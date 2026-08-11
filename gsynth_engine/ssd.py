"""Small Sequence Design — the heart of G-Synth.

Turns a gene or peptide-coding sequence into the forward and reverse oligos
that, once annealed, form a duplex whose sticky ends match a chosen
restriction pair in the target vector.

The logic below is a faithful extraction of the validated G-Synth
implementation. It reproduces the worked examples in the specification base
for base (see tests/test_ssd_golden.py) — treat those tests as the
definition of correct, not this code.

Two paths:

**Coding** — the input already carries its own ATG. For NdeI the leading ATG
is dropped, because NdeI's site CATATG supplies it: the forward oligo starts
"TATG", the T completing CATATG on ligation. The stop codon can be removed
so the gene can be read through into a C-terminal tag.

**Non-coding** — the input has no initiator. G-Synth prepends the standard
cassette: ATG (unless the enzyme supplies one) → linker → 6×His → linker →
protease site → insert.

In both cases the reverse oligo is built so that after annealing the duplex
presents exactly the two enzyme overhangs — e.g. 5'-TA for NdeI on the left
and 5'-TCGA for XhoI on the right — and nothing else.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from gsynth_engine.constants import (
    ALL_ENZYMES,
    CLEAVAGE_SITES,
    HIS_TAG,
    LEFT_LINKER,
    RESTRICTION_ENZYMES,
    RIGHT_LINKER,
    STOP_CODONS,
    left_remainders,
    overhang,
    right_remainders,
)
from gsynth_engine.sequence import (
    SequenceError,
    gc_content,
    reverse_complement,
    validate_dna,
)
from gsynth_engine.thermo import ANNEALING, melting_temperature


@dataclass(frozen=True)
class Segment:
    """One labelled stretch of the forward oligo, for display and reports."""
    name: str
    start: int
    end: int
    sequence: str


@dataclass
class SSDResult:
    """One small single-stranded DNA design: the two oligos and what they mean.

    `forward` and `reverse` are what you order — each written 5'→3' in its
    own direction, so `reverse` is *not* the reverse complement of `forward`
    read backwards; annealed, they form the cassette with the sticky ends the
    chosen enzyme pair leaves at either end.

    `segments` records which part of the design each stretch of bases came
    from — site, start codon, tag, linker, protease site, insert — so the
    interface can colour the construct without re-deriving the layout.
    """

    forward: str
    reverse: str
    is_coding: bool
    left_enzyme: str
    right_enzyme: str
    cleavage_site: str | None
    segments: list[Segment] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    # ── derived properties, computed on demand ──────────────────────────────
    @property
    def forward_length(self) -> int:
        return len(self.forward)

    @property
    def reverse_length(self) -> int:
        return len(self.reverse)

    @property
    def forward_gc(self) -> float:
        return round(gc_content(self.forward), 1)

    @property
    def reverse_gc(self) -> float:
        return round(gc_content(self.reverse), 1)

    @property
    def forward_tm(self) -> float:
        """Tm under the annealing reaction, not a generic primer dilution."""
        return round(melting_temperature(self.forward, conditions=ANNEALING), 1)

    @property
    def reverse_tm(self) -> float:
        return round(melting_temperature(self.reverse, conditions=ANNEALING), 1)

    @property
    def left_overhang(self) -> str:
        """The single-stranded end the duplex presents on the left."""
        return overhang(self.left_enzyme)[0]

    @property
    def right_overhang(self) -> str:
        """The 5' single-stranded end the duplex presents on the right.

        Given in bottom-strand 5'→3' sense, which is how it anneals to the
        cut vector.
        """
        return overhang(self.right_enzyme)[0]

    @property
    def orf_start(self) -> int:
        """Index in `forward` where the reading frame begins.

        Worth stating explicitly, because NdeI is counter-intuitive: its
        forward remainder is TATG, and the ATG occupies indices 1–3. The A of
        the 5'-TA overhang is *also* the A of the start codon — the overhang
        and the initiator overlap by one base. Every other enzyme gets an
        explicit ATG appended after its remainder, so the frame starts right
        after it.
        """
        cut, _ = left_remainders(self.left_enzyme)
        if cut.endswith("ATG"):
            return len(cut) - 3
        return len(cut)

    @property
    def coding_region(self) -> str:
        """`forward` from the start codon onwards."""
        return self.forward[self.orf_start :]


def _split_pair(enzyme_pair: str) -> tuple[str, str]:
    parts = [p.strip() for p in enzyme_pair.replace(" ", "").split("/")]
    if len(parts) != 2 or not all(parts):
        raise SequenceError(
            f"Enzyme pair must look like 'NdeI / XhoI' — got {enzyme_pair!r}."
        )
    for name in parts:
        if name not in ALL_ENZYMES:
            # 109 names is not a useful error message; list the ones this
            # lab actually keeps, and say how many others are known.
            known = (", ".join(sorted(RESTRICTION_ENZYMES))
                     + f", and {len(ALL_ENZYMES) - len(RESTRICTION_ENZYMES)} more")
            raise SequenceError(f"Unknown restriction enzyme {name!r}. Known: {known}")
    return parts[0], parts[1]


def _strip_stop_codon(sequence: str) -> tuple[str, bool]:
    """Truncate at the first in-frame stop codon. Returns (sequence, found)."""
    for i in range(0, len(sequence) - 2, 3):
        if sequence[i : i + 3] in STOP_CODONS:
            return sequence[:i], True
    return sequence, False


def design_small_sequence(
    sequence: str,
    *,
    enzyme_pair: str = "NdeI / XhoI",
    is_coding: bool = False,
    remove_stop: bool = False,
    cleavage_site: str | None = "Thrombin",
    include_his_tag: bool = True,
    include_linkers: bool = True,
) -> SSDResult:
    """Design the forward and reverse oligos for one insert.

    Args:
        sequence: the insert, A/C/G/T only.
        enzyme_pair: e.g. "NdeI / XhoI" — left enzyme first.
        is_coding: the insert already starts with its own ATG.
        remove_stop: truncate at the first in-frame stop (coding only), so a
            C-terminal tag stays in frame.
        cleavage_site: protease site inserted before the insert
            (non-coding only). None or "" for no site.
        include_his_tag / include_linkers: non-coding only.

    Raises:
        SequenceError: with a message meant for the person at the bench.
    """
    seq = validate_dna(sequence, field="insert sequence")
    left, right = _split_pair(enzyme_pair)
    warnings: list[str] = []

    if cleavage_site and cleavage_site not in CLEAVAGE_SITES:
        raise SequenceError(
            f"Unknown cleavage site {cleavage_site!r}. Known: "
            + ", ".join(sorted(CLEAVAGE_SITES))
        )

    fwd_cut_l, rev_cut_l = left_remainders(left)
    fwd_cut_r, rev_cut_r = right_remainders(right)

    segments: list[Segment] = []

    def add(name: str, part: str, cursor: int) -> int:
        if part:
            segments.append(Segment(name, cursor, cursor + len(part), part))
        return cursor + len(part)

    if is_coding:
        # NdeI supplies the ATG through its own recognition site; keeping the
        # insert's ATG as well would duplicate the start codon.
        if left == "NdeI" and seq.startswith("ATG"):
            seq = seq[3:]
            warnings.append(
                "The insert's ATG was removed: NdeI's site CATATG already "
                "provides the start codon."
            )
        if remove_stop:
            seq, found = _strip_stop_codon(seq)
            if not found:
                warnings.append("No in-frame stop codon was found to remove.")
            elif not seq:
                raise SequenceError(
                    "Removing the stop codon left nothing — check the reading frame."
                )

        forward = fwd_cut_l + seq + fwd_cut_r
        reverse = rev_cut_r + reverse_complement(seq) + rev_cut_l

        cursor = add(f"{left} overhang", fwd_cut_l, 0)
        cursor = add("insert", seq, cursor)
        add(f"{right} overhang", fwd_cut_r, cursor)

    else:
        # Non-coding: build the standard expression cassette in front.
        atg = "" if left == "NdeI" else "ATG"
        if left == "NdeI":
            warnings.append(
                "No separate ATG was added: NdeI's site CATATG provides the "
                "start codon."
            )

        cassette = ""
        if include_linkers and include_his_tag:
            cassette = LEFT_LINKER + HIS_TAG + RIGHT_LINKER
        elif include_his_tag:
            cassette = HIS_TAG
        elif include_linkers:
            cassette = LEFT_LINKER + RIGHT_LINKER

        cleavage_seq = CLEAVAGE_SITES[cleavage_site] if cleavage_site else ""

        forward = fwd_cut_l + atg + cassette + cleavage_seq + seq + fwd_cut_r

        # The reverse oligo mirrors the forward, element by element, in
        # reverse order — so that the annealed duplex is flush everywhere
        # except the two enzyme overhangs.
        reverse = rev_cut_r + reverse_complement(seq)
        if cleavage_seq:
            reverse += reverse_complement(cleavage_seq)
        if include_linkers and include_his_tag:
            reverse += reverse_complement(RIGHT_LINKER)
            reverse += reverse_complement(HIS_TAG)
            reverse += reverse_complement(LEFT_LINKER)
        elif include_his_tag:
            reverse += reverse_complement(HIS_TAG)
        elif include_linkers:
            reverse += reverse_complement(RIGHT_LINKER)
            reverse += reverse_complement(LEFT_LINKER)
        reverse += reverse_complement(atg) + rev_cut_l

        cursor = add(f"{left} overhang", fwd_cut_l, 0)
        cursor = add("start codon", atg, cursor)
        if include_linkers and include_his_tag:
            cursor = add("linker", LEFT_LINKER, cursor)
            cursor = add("6×His tag", HIS_TAG, cursor)
            cursor = add("linker", RIGHT_LINKER, cursor)
        elif include_his_tag:
            cursor = add("6×His tag", HIS_TAG, cursor)
        elif include_linkers:
            cursor = add("linker", LEFT_LINKER, cursor)
            cursor = add("linker", RIGHT_LINKER, cursor)
        if cleavage_seq:
            cursor = add(f"{cleavage_site} site", cleavage_seq, cursor)
        cursor = add("insert", seq, cursor)
        add(f"{right} overhang", fwd_cut_r, cursor)

    # Warn about anything that would bite at the bench.
    for name, enzyme in ((left, left), (right, right)):
        site = str(ALL_ENZYMES[enzyme]["recognition"])
        # The site is expected once at each end after ligation; an extra copy
        # inside the insert would be cut during cloning.
        if site and site in seq:
            warnings.append(
                f"The insert contains an internal {name} site ({site}) — it "
                f"will be cut during cloning. Choose another enzyme or "
                f"silently mutate the site."
            )

    return SSDResult(
        forward=forward,
        reverse=reverse,
        is_coding=is_coding,
        left_enzyme=left,
        right_enzyme=right,
        cleavage_site=cleavage_site or None,
        segments=segments,
        warnings=warnings,
    )
