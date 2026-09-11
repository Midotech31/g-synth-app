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

if TYPE_CHECKING:
    from gsynth_engine.vectors import VectorSpec


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

    seq = clean_dna(sequence)
    return "".join(
        _CODONS.get(seq[i : i + 3], "X") for i in range(0, len(seq) - 2, 3)
    )


@dataclass(frozen=True)
class End:


    sequence: str
    strand: str
    side: str = "left"

    @property
    def kind(self) -> str:

        if self.strand == "blunt":
            return "blunt"
        protrudes_at_its_own_five_prime = (
            (self.side == "left" and self.strand == "top")
            or (self.side == "right" and self.strand == "bottom")
        )
        return "5'" if protrudes_at_its_own_five_prime else "3'"

    def anneals_to(self, other: End) -> bool:

        if self.strand == "blunt" or other.strand == "blunt":
            return self.strand == other.strand
        return self.sequence == other.sequence and self.strand != other.strand


def find_sites(sequence: str, enzyme: str, *, circular: bool = True) -> list[int]:

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


    top: str
    left_end: End
    right_end: End

    vector_start: int

    removed_start: int
    removed_end: int
    removed_length: int
    circular_source: bool


    reversed_insert: bool = False

    @property
    def length(self) -> int:
        return len(self.top)


def _cut_positions(enzyme: str, site_start: int) -> tuple[int, int]:

    info = ALL_ENZYMES[enzyme]
    return (
        site_start + int(info["cut_top"]),      # type: ignore[arg-type]
        site_start + int(info["cut_bottom"]),   # type: ignore[arg-type]
    )


def _end_at(sequence: str, top_cut: int, bottom_cut: int, *, side: str) -> End:


    fragment_side = "left" if side == "downstream" else "right"

    if top_cut == bottom_cut:
        return End("", "blunt", fragment_side)

    lo, hi = min(top_cut, bottom_cut), max(top_cut, bottom_cut)


    overhang = "".join(sequence[i % len(sequence)] for i in range(lo, hi))


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

        left = _cut_positions(left_enzyme, find_sites(sequence, left_enzyme)[0])[0]
        right = _cut_positions(right_enzyme, find_sites(sequence, right_enzyme)[0])[0]
        return (right - left) % length

    flipped = reverse_complement(seq)


    reversed_insert = arc(flipped) < arc(seq)
    working = flipped if reversed_insert else seq

    left_site = find_sites(working, left_enzyme, circular=True)[0]
    right_site = find_sites(working, right_enzyme, circular=True)[0]
    left_top, left_bottom = _cut_positions(left_enzyme, left_site)
    right_top, right_bottom = _cut_positions(right_enzyme, right_site)


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


    top: str
    bottom: str
    left_end: End
    right_end: End


    trimmed_left: int
    trimmed_right: int

    @property
    def length(self) -> int:

        return len(self.top)


def digest_linear(
    fragment: str,
    *,
    left_enzyme: str,
    right_enzyme: str,
) -> Digest:

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
            f"The {left_enzyme} site is not to the left of the {right_enzyme} "
            f"site in this product. The enzymes are the wrong way round for "
            f"these primers — swap them, or swap the tails."
        )


    top = working[left_top:right_top]
    bottom = reverse_complement(working[left_bottom:right_bottom])

    return Digest(
        top=top,
        bottom=bottom,


        left_end=_end_at(working, left_top, left_bottom, side="downstream"),
        right_end=_end_at(working, right_top, right_bottom, side="upstream"),


        trimmed_left=left_top,
        trimmed_right=len(working) - right_top,
    )


@dataclass(frozen=True)
class Junction:


    name: str
    enzyme: str
    overhang: str
    kind: str
    position: int
    context: str
    site_regenerated: bool


@dataclass(frozen=True)
class TagOutcome:


    name: str
    end: str
    present: bool
    position: int | None = None
    note: str = ""


@dataclass(frozen=True)
class ReadingFrameCheck:


    code: str
    label: str
    status: str
    detail: str


@dataclass(frozen=True)
class ReadingFrameAssessment:


    status: str
    summary: str
    translation_start: int | None = None
    start_codon: str | None = None
    start_source: str | None = None
    rbs_name: str | None = None
    rbs_start: int | None = None
    rbs_end: int | None = None
    rbs_spacing_nt: int | None = None
    rbs_source: str | None = None
    promoter_name: str | None = None
    promoter_start: int | None = None
    promoter_source: str | None = None
    stop_codon: str | None = None
    stop_position: int | None = None
    stop_context: str | None = None
    left_junction_offset: int | None = None
    right_junction_phase: int | None = None
    protein_length: int = 0
    checks: tuple[ReadingFrameCheck, ...] = ()

    @property
    def confirmed(self) -> bool:
        return self.status == "pass"


@dataclass
class CloningResult:


    plasmid: str
    name: str
    insert_start: int
    insert_end: int
    backbone_length: int
    removed_length: int
    left_enzyme: str
    right_enzyme: str
    junctions: list[Junction] = field(default_factory=list)
    annotations: list[dict] = field(default_factory=list)
    protein: str = ""


    reversed_insert: bool = False
    tags: list[TagOutcome] = field(default_factory=list)
    reading_frame: ReadingFrameAssessment = field(
        default_factory=lambda: ReadingFrameAssessment(
            status="review",
            summary="The reading frame has not been assessed.",
        )
    )
    warnings: list[str] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)


    translation_start: int | None = None

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

        return not self.problems


def _flip_annotations(annotations: list[dict], length: int) -> list[dict]:

    flipped: list[dict] = []
    for feature in annotations:
        entry = dict(feature)
        start, end = int(feature.get("start", 0)), int(feature.get("end", 0))
        mirrored = (length - end) % length
        shift = mirrored - (length - end)
        entry["start"] = mirrored
        entry["end"] = length - start + shift
        entry["direction"] = -int(feature.get("direction", 1))
        if feature.get('translation_start') is not None and feature.get('translation_end') is not None:
            entry['translation_start'] = length - int(feature['translation_end']) + shift
            entry['translation_end'] = length - int(feature['translation_start']) + shift
        flipped.append(entry)
    return flipped


def _remap_annotations(
    annotations: list[dict], backbone: Backbone, plasmid_length: int
) -> list[dict]:

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

        new_start = (start - backbone.vector_start) % (span + backbone.removed_length)
        new_end = new_start + (end - start)

        if new_start >= span:
            continue
        entry = dict(feature)
        entry["start"] = new_start
        if new_end > span:
            entry["end"] = span
            entry["truncated"] = True
        else:
            entry["end"] = new_end
        entry["end"] = min(entry["end"], plasmid_length)
        if feature.get('translation_start') is not None and feature.get('translation_end') is not None:
            translated_start = new_start + int(feature['translation_start']) - start
            translated_end = new_start + int(feature['translation_end']) - start
            clipped_end = min(translated_end, entry['end'])
            if entry.get('direction') == -1:

                clipped_end -= (clipped_end - translated_end) % 3
            if translated_start < clipped_end:
                entry['translation_start'], entry['translation_end'] = translated_start, clipped_end
            else:
                entry.pop('translation_start', None)
                entry.pop('translation_end', None)
                if entry.get('type') == 'CDS':
                    entry['type'] = 'misc_feature'
                    entry['basis'] = 'The coding region was removed; only a flanking fragment remains.'
        moved.append(entry)
    return moved


def _infer_insert_orf_start(insert: str, *, search_limit: int = 30) -> int | None:

    candidates: list[int] = []
    limit = min(max(0, search_limit), max(0, len(insert) - 2))
    for position in range(limit):
        if insert[position : position + 3] != "ATG":
            continue
        translated = translate(insert[position:])
        first_stop = translated.find("*")
        if first_stop == -1 or first_stop >= max(0, len(translated) - 2):
            candidates.append(position)
    return candidates[0] if len(candidates) == 1 else None


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
    auto_detect_frame: bool = False,
) -> CloningResult:

    insert_top = validate_dna(insert, field="insert")
    start_source = "declared" if orf_start is not None else None
    if orf_start is None and auto_detect_frame:
        orf_start = _infer_insert_orf_start(insert_top)
        if orf_start is not None:
            start_source = "sequence_candidate"
    backbone = linearise(
        vector, left_enzyme=left_enzyme, right_enzyme=right_enzyme, circular=circular
    )

    problems: list[str] = []
    warnings: list[str] = []

    if insert_left_end and insert_right_end:
        left, right = insert_left_end, insert_right_end
    elif insert_reverse:
        left, right = _observed_insert_ends(insert_top, insert_reverse, left_enzyme)
        mismatches = insert_duplex_mismatches(insert_top, insert_reverse, left_enzyme)
        if mismatches:
            preview = ", ".join(str(position + 1) for position in mismatches[:5])
            suffix = "…" if len(mismatches) > 5 else ""
            problems.append(
                f"The two supplied insert strands do not pair at {len(mismatches)} "
                f"position(s): {preview}{suffix}."
            )
    else:


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


    plasmid = backbone.top + insert_top
    insert_start = backbone.length
    insert_end = insert_start + len(insert_top)
    remapped_annotations = _remap_annotations(
        vector_annotations or [], backbone, len(plasmid)
    )


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
                f"{'s' if len(internal) > 1 else ''}. Extended Sequence Design never "
                f"digests the insert, so the build is unaffected — but a "
                f"diagnostic digest with {enzyme} will cut inside the gene."
            )

    protein = ""
    stop_at: int | None = None
    translation_start: int | None = None
    if orf_start is not None and 0 <= orf_start < len(insert_top):

        translation_start = insert_start + orf_start
        protein, stop_at = _translate_in_plasmid(plasmid, translation_start)

        if protein and protein[0] != "M":
            warnings.append(
                "The reading frame does not start with a methionine — check "
                "the ATG and the left-hand enzyme."
            )


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


    if translation_start is not None:
        from gsynth_engine.annotations import detect_common_features

        context = [*remapped_annotations, {
            "name": name, "type": "CDS", "start": translation_start,
            "end": insert_end, "direction": 1,
        }]
        for match in detect_common_features(plasmid, circular=True, existing=context):
            annotation = match["annotation"]
            if (annotation["type"] == "RBS" and annotation["direction"] == 1
                    and 4 <= (translation_start - annotation["end"]) % len(plasmid) <= 14):
                remapped_annotations.append(annotation)

    reading_frame = _assess_reading_frame(
        plasmid,
        insert_start=insert_start,
        insert_end=insert_end,
        orf_start=orf_start,
        protein=protein,
        stop_at=stop_at,
        annotations=remapped_annotations,
        vector_spec=vector_spec,
        tags=tags,
        existing_problems=problems,
        start_source=start_source,
    )

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
        annotations=remapped_annotations,
        protein=protein,
        reversed_insert=backbone.reversed_insert,
        translation_start=translation_start,
        tags=tags,
        reading_frame=reading_frame,
        warnings=warnings,
        problems=problems,
    )


def _nearest_upstream_feature(
    annotations: list[dict],
    position: int,
    plasmid_length: int,
    *,
    feature_types: set[str],
    name_terms: tuple[str, ...],
    max_distance: int,
) -> tuple[dict | None, int | None]:

    candidates: list[tuple[int, dict]] = []
    for feature in annotations:
        feature_type = str(feature.get("type", "")).lower().replace("-", "_")
        if feature_type == "regulatory":
            feature_type = str(feature.get("regulatory_class", "")).lower()
            if feature_type == "ribosome_binding_site":
                feature_type = "rbs"
        name = str(feature.get("name", "")).lower()
        if feature_type not in feature_types and not any(term in name for term in name_terms):
            continue
        if int(feature.get("direction", 0) or 0) == -1:
            continue
        end = int(feature.get("end", 0)) % plasmid_length
        distance = (position - end) % plasmid_length
        if distance <= max_distance:
            candidates.append((distance, feature))
    if not candidates:
        return None, None
    distance, feature = min(candidates, key=lambda candidate: candidate[0])
    return feature, distance


def _assess_reading_frame(
    plasmid: str,
    *,
    insert_start: int,
    insert_end: int,
    orf_start: int | None,
    protein: str,
    stop_at: int | None,
    annotations: list[dict],
    vector_spec: VectorSpec | None,
    tags: list[TagOutcome],
    existing_problems: list[str],
    start_source: str | None,
) -> ReadingFrameAssessment:

    if vector_spec is not None and not vector_spec.expression_capable:
        check = ReadingFrameCheck(
            "FRAME_NOT_APPLICABLE",
            "Expression context",
            "pass",
            f"{vector_spec.name} is a cloning backbone; insert expression is not claimed.",
        )
        return ReadingFrameAssessment(
            status="not_applicable",
            summary="Reading-frame validation is not applicable to this cloning backbone.",
            checks=(check,),
        )

    if orf_start is None or not 0 <= orf_start < insert_end - insert_start:
        check = ReadingFrameCheck(
            "FRAME_START_UNDECLARED",
            "Translation start",
            "review",
            "No coding start was declared, so G-Synth cannot assign a reading frame.",
        )
        return ReadingFrameAssessment(
            status="review",
            summary="Reading frame not confirmed: a coding start must be supplied.",
            checks=(check,),
        )

    length = len(plasmid)
    start = insert_start + orf_start
    start_codon = "".join(plasmid[(start + offset) % length] for offset in range(3))
    checks: list[ReadingFrameCheck] = []

    start_ok = start_codon == "ATG"
    start_status = "pass" if start_ok and start_source == "declared" else (
        "review" if start_ok else "block"
    )
    checks.append(ReadingFrameCheck(
        "FRAME_START_CODON",
        "Initiating codon",
        start_status,
        (
            f"ATG at recombinant base {start + 1} fixes the downstream triplet phase."
            if start_ok and start_source == "declared" else
            f"A single ATG candidate was detected at recombinant base {start + 1}; confirm it as the intended start."
            if start_ok else
            f"{start_codon or 'No complete codon'} occurs at the declared start; ATG is required by this design."
        ),
    ))

    rbs, rbs_spacing = _nearest_upstream_feature(
        annotations,
        start,
        length,
        feature_types={"rbs", "ribosome_binding_site"},
        name_terms=("rbs", "shine-dalgarno", "shine dalgarno"),
        max_distance=40,
    )
    rbs_status = "review"
    if rbs is not None and rbs_spacing is not None:
        if 4 <= rbs_spacing <= 14:
            rbs_status = (
                "review"
                if rbs.get("inferred") or int(rbs.get("direction", 0) or 0) == 0
                else "pass"
            )
            rbs_detail = (
                f"{rbs.get('name', 'RBS')} ends {rbs_spacing} nt before the initiating ATG."
                + (" It was detected by exact sequence matching and requires confirmation." if rbs.get("inferred") else "")
            )
        elif 1 <= rbs_spacing <= 25:
            rbs_status = "review"
            rbs_detail = (
                f"The nearest annotated RBS ends {rbs_spacing} nt before the ATG; "
                "confirm this host- and RBS-specific spacing experimentally."
            )
        else:
            rbs_status = "block"
            rbs_detail = (
                f"The nearest annotated RBS is {rbs_spacing} nt from the initiating ATG; "
                "it is outside the accepted initiation-context window."
            )
    else:
        if vector_spec is not None and not vector_spec.supplies_translation_start:
            rbs_status = "block"
            rbs_detail = (
                f"{vector_spec.name} does not supply an RBS/start module, and the generated insert "
                "contains an ATG but no annotated RBS. Add a validated translation-initiation module."
            )
        else:
            rbs_detail = (
                "No RBS is annotated in the 40 nt upstream of the initiating ATG; "
                "translation initiation cannot be confirmed from this vector record."
            )
    checks.append(ReadingFrameCheck(
        "FRAME_RBS_CONTEXT",
        "Ribosome-binding site",
        rbs_status,
        rbs_detail,
    ))

    promoter, promoter_distance = _nearest_upstream_feature(
        annotations,
        start,
        length,
        feature_types={"promoter"},
        name_terms=("promoter",),
        max_distance=500,
    )
    checks.append(ReadingFrameCheck(
        "FRAME_PROMOTER_CONTEXT",
        "Upstream promoter",
        (
            "review"
            if promoter is not None and (
                promoter.get("inferred") or int(promoter.get("direction", 0) or 0) == 0
            )
            else "pass" if promoter is not None else "review"
        ),
        (
            f"{promoter.get('name', 'Promoter')} is annotated {promoter_distance} nt upstream of the start."
            + (" The feature was detected by exact sequence matching and requires confirmation." if promoter.get("inferred") else "")
            if promoter is not None else
            "No promoter is annotated within 500 nt upstream; transcriptional context is unconfirmed."
        ),
    ))

    left_offset = orf_start
    right_phase = (insert_end - start) % 3
    checks.append(ReadingFrameCheck(
        "FRAME_JUNCTION_PHASE",
        "Cloning-junction phase",
        start_status,
        (
            f"The ATG begins {left_offset} nt after the left junction; the right junction "
            f"occurs {right_phase} nt into its codon relative to that ATG."
        ),
    ))

    stop_context: str | None = None
    stop_codon: str | None = None
    terminus_status = "pass"
    if stop_at is None:
        terminus_status = "block"
        terminus_detail = "No in-frame stop codon occurs before translation would traverse the circular plasmid."
    else:
        stop_codon = "".join(plasmid[(stop_at + offset) % length] for offset in range(3))
        reached = (stop_at - start) % length
        to_insert_end = insert_end - start
        if reached < max(0, to_insert_end - 6):
            stop_context = "premature_insert"
            terminus_status = "block"
            terminus_detail = f"{stop_codon} terminates translation prematurely at residue {reached // 3 + 1}."
        elif reached < to_insert_end:
            stop_context = "insert"
            terminus_detail = f"The insert terminates with {stop_codon} after {len(protein)} residues."
        else:
            stop_context = "vector"
            c_terminal_tags = [tag for tag in tags if tag.end == "C"]
            missing_expected_tag = c_terminal_tags and not any(tag.present for tag in c_terminal_tags)
            if missing_expected_tag:
                terminus_status = "block"
                terminus_detail = (
                    "Translation reaches vector sequence, but no declared C-terminal vector tag is in frame."
                )
            else:
                terminus_detail = (
                    f"Translation crosses the right junction and terminates at vector base {stop_at + 1}"
                    f" after {len(protein)} residues."
                )
    if any("truncated" in problem.lower() for problem in existing_problems):
        terminus_status = "block"
    checks.append(ReadingFrameCheck(
        "FRAME_TRANSLATED_TERMINUS",
        "Translated terminus",
        terminus_status,
        terminus_detail,
    ))

    statuses = {check.status for check in checks}
    status = "block" if "block" in statuses else "review" if "review" in statuses else "pass"
    summary = {
        "pass": "Expression reading frame confirmed from the annotated RBS and ATG through the translated terminus.",
        "review": "Reading frame is plausible but not confirmed because expression-context evidence is incomplete.",
        "block": "Expression reading frame is invalid; correct the flagged initiation or translation evidence before use.",
    }[status]
    return ReadingFrameAssessment(
        status=status,
        summary=summary,
        translation_start=start,
        start_codon=start_codon,
        start_source=start_source,
        rbs_name=str(rbs.get("name")) if rbs is not None else None,
        rbs_start=int(rbs.get("start", 0)) if rbs is not None else None,
        rbs_end=int(rbs.get("end", 0)) if rbs is not None else None,
        rbs_spacing_nt=rbs_spacing,
        rbs_source="sequence_motif" if rbs is not None and rbs.get("inferred") else "annotation" if rbs is not None else None,
        promoter_name=str(promoter.get("name")) if promoter is not None else None,
        promoter_start=int(promoter.get("start", 0)) if promoter is not None else None,
        promoter_source="sequence_motif" if promoter is not None and promoter.get("inferred") else "annotation" if promoter is not None else None,
        stop_codon=stop_codon,
        stop_position=stop_at,
        stop_context=stop_context,
        left_junction_offset=left_offset,
        right_junction_phase=right_phase,
        protein_length=len(protein),
        checks=tuple(checks),
    )


def _tag_outcomes(
    protein: str, spec: VectorSpec, *, insert_residues: int, upstream_residues: int,
) -> list[TagOutcome]:

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

    sequence, kind = enzyme_overhang(enzyme)
    if kind == "blunt":
        return End("", "blunt", side)


    if side == "left":
        return End(sequence, "top" if kind == "5'" else "bottom", "left")
    return End(sequence, "bottom" if kind == "5'" else "top", "right")


def _observed_insert_ends(
    top: str, reverse: str, left_enzyme: str,
) -> tuple[End, End]:

    info = ALL_ENZYMES[left_enzyme]
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


def insert_duplex_mismatches(top: str, reverse: str, left_enzyme: str) -> list[int]:

    top_clean = clean_dna(top)
    bottom_sense = reverse_complement(clean_dna(reverse))
    info = ALL_ENZYMES[left_enzyme]
    offset = int(info["cut_bottom"]) - int(info["cut_top"])  # type: ignore[arg-type]
    overlap_start = max(0, offset)
    overlap_end = min(len(top_clean), offset + len(bottom_sense))
    return [
        position
        for position in range(overlap_start, overlap_end)
        if top_clean[position] != bottom_sense[position - offset]
    ]


def _junction(
    name: str, enzyme: str, end: End, plasmid: str, position: int,
) -> Junction:

    window = 12
    length = len(plasmid)
    context = "".join(
        plasmid[(position + offset) % length]
        for offset in range(-window, window)
    )
    site: str = ALL_ENZYMES[enzyme]["recognition"]  # type: ignore[index]
    return Junction(
        name=name,
        enzyme=enzyme,
        overhang=end.sequence,
        kind=end.kind,
        position=position,
        context=context,


        site_regenerated=site in context or reverse_complement(site) in context,
    )


def open_reading_frames(
    sequence: str, *, minimum_codons: int = 30, circular: bool = True,
) -> list[dict]:

    seq = clean_dna(sequence)
    if len(seq) < 6:
        return []


    scan = seq + seq if circular else seq
    found: list[dict] = []

    for frame in range(3):
        start = None
        for i in range(frame, len(scan) - 2, 3):
            codon = scan[i : i + 3]
            if start is None and codon == "ATG":
                if i >= len(seq):
                    break
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
