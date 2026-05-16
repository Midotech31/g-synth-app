"""FASTA reader and writer."""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, TextIO

from gsynth_core.io.records import SeqRecord
from gsynth_core.sequence.ops import clean_dna


def _parse_stream(stream: Iterable[str]) -> list[SeqRecord]:
    records: list[SeqRecord] = []
    name, desc, buf = "", "", []

    def flush() -> None:
        if name or buf:
            records.append(SeqRecord(
                sequence=clean_dna("".join(buf), allow_ambiguous=True),
                name=name, description=desc,
            ))

    for raw in stream:
        line = raw.rstrip("\r\n")
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            parts = header.split(None, 1)
            name = parts[0] if parts else ""
            desc = parts[1] if len(parts) > 1 else ""
            buf = []
        else:
            buf.append(line.strip())
    flush()
    return records


def parse_fasta(source: str | Path | TextIO) -> list[SeqRecord]:
    """Parse a FASTA file or string."""
    if isinstance(source, (str, Path)):
        p = Path(source)
        if p.exists() and p.is_file():
            with p.open("r", encoding="utf-8") as f:
                return _parse_stream(f)
        if isinstance(source, str) and ">" in source:
            return _parse_stream(source.splitlines(keepends=True))
        raise FileNotFoundError(f"Not a valid FASTA file or string: {source!r}")
    return _parse_stream(source)


def write_fasta(records: SeqRecord | list[SeqRecord],
                destination: str | Path | TextIO | None = None,
                *, line_width: int = 60) -> str:
    """Serialize records as FASTA. Returns the string (and writes if destination given)."""
    if isinstance(records, SeqRecord):
        records = [records]
    chunks: list[str] = []
    for rec in records:
        header = f">{rec.name}"
        if rec.description:
            header += f" {rec.description}"
        chunks.append(header)
        for i in range(0, len(rec.sequence), line_width):
            chunks.append(rec.sequence[i:i+line_width])
    out = "\n".join(chunks) + "\n"
    if destination is not None:
        if isinstance(destination, (str, Path)):
            Path(destination).write_text(out, encoding="utf-8")
        else:
            destination.write(out)
    return out
