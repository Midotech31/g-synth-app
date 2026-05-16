"""
SnapGene `.dna` file reader.

Uses the third-party `snapgene-reader` package when available. The format
is proprietary and only partially documented; we extract:
- the bare sequence
- the topology (circular / linear)
- features (label, location, type)

For writing SnapGene files, we use GenBank as an intermediate (SnapGene
itself opens GenBank natively).
"""

from __future__ import annotations

from pathlib import Path

from gsynth_core.io.records import Feature, FeatureLocation, SeqRecord


def _try_snapgene_reader():
    try:
        import snapgene_reader  # type: ignore
        return snapgene_reader
    except ImportError:
        return None


def parse_snapgene(path: str | Path) -> SeqRecord:
    """
    Parse a SnapGene `.dna` file into a :class:`SeqRecord`.

    Requires the optional `snapgene-reader` package (`pip install
    snapgene-reader`). If unavailable, raises ImportError with a helpful
    message.
    """
    reader = _try_snapgene_reader()
    if reader is None:
        raise ImportError(
            "Reading SnapGene `.dna` files requires the optional `snapgene-reader` package. "
            "Install it with: pip install snapgene-reader"
        )

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    raw = reader.snapgene_file_to_dict(str(p))
    seq = raw.get("seq", "").upper()
    dna_info = raw.get("dna")
    topology = "linear"
    if isinstance(dna_info, dict):
        if dna_info.get("topology") == "circular":
            topology = "circular"
    elif isinstance(dna_info, str) and dna_info.strip().lower() == "circular":
        topology = "circular"

    features: list[Feature] = []
    for f in raw.get("features", []):
        try:
            # snapgene-reader returns features with 'start', 'end' (1-based, inclusive),
            # 'name', 'type', and 'strand' fields.
            start = int(f.get("start", 0)) - 1
            end = int(f.get("end", 0))
            strand = 1 if str(f.get("strand", "1")) in ("1", "+") else -1
            ftype = f.get("type") or "misc_feature"
            features.append(Feature(
                type=ftype,
                location=FeatureLocation(start=max(0, start), end=end, strand=strand),
                qualifiers={
                    "label": f.get("name") or "",
                    "color": f.get("color", ""),
                    "note": f.get("note", ""),
                },
            ))
        except Exception:
            continue

    return SeqRecord(
        sequence=seq,
        name=p.stem,
        description=raw.get("notes", {}).get("Description", "") if isinstance(raw.get("notes"), dict) else "",
        features=features,
        topology=topology,
        molecule_type="DNA",
    )
