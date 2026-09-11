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
    supplies_start_codon,
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

    name: str
    start: int
    end: int
    sequence: str


@dataclass
class SSDResult:


    forward: str
    reverse: str
    is_coding: bool
    left_enzyme: str
    right_enzyme: str
    cleavage_site: str | None
    segments: list[Segment] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


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

        return round(melting_temperature(self.forward, conditions=ANNEALING), 1)

    @property
    def reverse_tm(self) -> float:
        return round(melting_temperature(self.reverse, conditions=ANNEALING), 1)

    @property
    def left_overhang(self) -> str:

        return overhang(self.left_enzyme)[0]

    @property
    def right_overhang(self) -> str:

        return overhang(self.right_enzyme)[0]

    @property
    def orf_start(self) -> int:

        cut, _ = left_remainders(self.left_enzyme)
        if supplies_start_codon(self.left_enzyme):
            return len(cut) - 3
        return len(cut)

    @property
    def coding_region(self) -> str:

        return self.forward[self.orf_start :]


def _split_pair(enzyme_pair: str) -> tuple[str, str]:
    parts = [p.strip() for p in enzyme_pair.replace(" ", "").split("/")]
    if len(parts) != 2 or not all(parts):
        raise SequenceError(
            f"Enzyme pair must look like 'NdeI / XhoI' — got {enzyme_pair!r}."
        )
    for name in parts:
        if name not in ALL_ENZYMES:


            known = (", ".join(sorted(RESTRICTION_ENZYMES))
                     + f", and {len(ALL_ENZYMES) - len(RESTRICTION_ENZYMES)} more")
            raise SequenceError(f"Unknown restriction enzyme {name!r}. Known: {known}")
    return parts[0], parts[1]


def _strip_stop_codon(sequence: str) -> tuple[str, bool]:

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
        if not seq.startswith("ATG"):
            raise SequenceError(
                "A coding insert must begin with ATG. Turn off 'already has "
                "its own ATG' if this is a non-coding sequence, or add the "
                "correct start codon."
            )


        if supplies_start_codon(left):
            seq = seq[3:]
            warnings.append(
                f"The insert's ATG was removed: {left}'s retained site "
                "already provides the start codon."
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
        atg = "" if supplies_start_codon(left) else "ATG"
        if supplies_start_codon(left):
            warnings.append(
                f"No separate ATG was added: {left}'s retained site provides "
                "the start codon."
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


    for name, enzyme in ((left, left), (right, right)):
        site = str(ALL_ENZYMES[enzyme]["recognition"])


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
