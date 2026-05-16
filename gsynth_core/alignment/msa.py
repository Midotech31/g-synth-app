"""Multiple Sequence Alignment.

Tries MAFFT → MUSCLE → Clustal Omega on PATH; falls back to a pure-Python
star alignment for ≤50 sequences.
"""
from __future__ import annotations
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from gsynth_core.alignment.pairwise import align_global


@dataclass(frozen=True, slots=True)
class MSAResult:
    names: list[str]
    aligned: list[str]
    method: str
    consensus: str = ""

    @property
    def length(self) -> int:
        return len(self.aligned[0]) if self.aligned else 0

    def pretty(self, width: int = 60) -> str:
        if not self.aligned:
            return "(empty MSA)"
        L = len(self.aligned[0])
        name_w = max(12, max(len(n) for n in self.names))
        lines = [f"MSA ({self.method}, {len(self.aligned)} sequences, {L} columns)"]
        for i in range(0, L, width):
            for name, seq in zip(self.names, self.aligned):
                lines.append(f"  {name:{name_w}s}  {seq[i:i+width]}")
            if self.consensus:
                lines.append(f"  {'consensus':{name_w}s}  {self.consensus[i:i+width]}")
            lines.append("")
        return "\n".join(lines)


def _parse_fasta(text: str) -> tuple[list[str], list[str]]:
    names, seqs, cur = [], [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
            names.append(line[1:].split()[0] if line[1:].strip() else f"s{len(names)+1}")
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return names, seqs


def _compute_consensus(aligned: list[str]) -> str:
    if not aligned:
        return ""
    L = len(aligned[0]); cons: list[str] = []
    for col in range(L):
        cb = [s[col] for s in aligned if s[col] != "-"]
        if not cb:
            cons.append("-"); continue
        counts: dict[str, int] = {}
        for b in cb:
            counts[b] = counts.get(b, 0) + 1
        top, n = max(counts.items(), key=lambda x: x[1])
        cons.append(top if n > len(cb) / 2 else "N")
    return "".join(cons)


def _run_external(tool: str, sequences: list[tuple[str, str]]) -> list[str] | None:
    if not shutil.which(tool):
        return None
    fasta_in = "\n".join(f">{n}\n{s}" for n, s in sequences)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        f.write(fasta_in)
        in_path = Path(f.name)
    try:
        if tool == "mafft":
            cmd = ["mafft", "--auto", "--quiet", str(in_path)]
        elif tool == "muscle":
            cmd = ["muscle", "-in", str(in_path)]
        elif tool == "clustalo":
            cmd = ["clustalo", "-i", str(in_path), "--outfmt=fa", "--force"]
        else:
            return None
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            return None
        _, out = _parse_fasta(r.stdout)
        return out
    except Exception:
        return None
    finally:
        in_path.unlink(missing_ok=True)


def _star_alignment(sequences: list[tuple[str, str]]) -> list[str]:
    """Star MSA: align all to the center sequence then propagate gaps."""
    if len(sequences) < 2:
        return [s for _, s in sequences]
    n = len(sequences); seqs = [s for _, s in sequences]
    best_center = 0; best_id = -1.0
    for i in range(n):
        total = 0.0
        for j in range(n):
            if i == j: continue
            total += align_global(seqs[i], seqs[j]).identity_percent
        m = total / (n - 1)
        if m > best_id:
            best_id = m; best_center = i
    center = seqs[best_center]
    pa: list[tuple[str, str]] = []
    for i in range(n):
        if i == best_center:
            pa.append((center, center))
        else:
            aln = align_global(center, seqs[i])
            pa.append((aln.aligned_b, aln.aligned_a))
    cents = [p[1] for p in pa]; oths = [p[0] for p in pa]
    out: list[list[str]] = [[] for _ in range(n)]
    pos = [0] * n; center_pos = 0
    while center_pos <= len(center):
        any_gap = any(pos[k] < len(cents[k]) and cents[k][pos[k]] == "-" for k in range(n))
        if any_gap:
            col = []
            for k in range(n):
                if pos[k] < len(cents[k]) and cents[k][pos[k]] == "-":
                    col.append(oths[k][pos[k]]); pos[k] += 1
                else:
                    col.append("-")
            for k in range(n):
                out[k].append(col[k])
            continue
        if center_pos == len(center): break
        col = []
        for k in range(n):
            if pos[k] < len(oths[k]):
                col.append(oths[k][pos[k]]); pos[k] += 1
            else:
                col.append("-")
        for k in range(n):
            out[k].append(col[k])
        center_pos += 1
    return ["".join(o) for o in out]


def multiple_alignment(sequences: list[tuple[str, str]], *,
                       method: Literal["auto", "mafft", "muscle", "clustalo", "star"] = "auto"
                       ) -> MSAResult:
    """Align multiple sequences. method='auto' picks the best available tool."""
    if not sequences:
        return MSAResult([], [], "empty", "")
    if len(sequences) == 1:
        n, s = sequences[0]
        return MSAResult([n], [s], "identity", s)
    if method == "auto":
        for tool in ("mafft", "muscle", "clustalo"):
            aligned = _run_external(tool, sequences)
            if aligned is not None:
                names = [n for n, _ in sequences]
                return MSAResult(names, aligned, tool, _compute_consensus(aligned))
        aligned = _star_alignment(sequences)
        names = [n for n, _ in sequences]
        return MSAResult(names, aligned, "star", _compute_consensus(aligned))
    if method in ("mafft", "muscle", "clustalo"):
        aligned = _run_external(method, sequences)
        if aligned is None:
            raise RuntimeError(f"{method} is not available on PATH")
        names = [n for n, _ in sequences]
        return MSAResult(names, aligned, method, _compute_consensus(aligned))
    if method == "star":
        aligned = _star_alignment(sequences)
        names = [n for n, _ in sequences]
        return MSAResult(names, aligned, "star", _compute_consensus(aligned))
    raise ValueError(f"Unknown method: {method}")
