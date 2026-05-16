"""Pairwise alignment: Needleman-Wunsch global + Smith-Waterman local (Gotoh affine)."""
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Alignment:
    aligned_a: str
    aligned_b: str
    score: float
    start_a: int
    end_a: int
    start_b: int
    end_b: int
    mode: str
    identity_percent: float

    def pretty(self, width: int = 60) -> str:
        a, b = self.aligned_a, self.aligned_b
        match = "".join("|" if x == y and x != "-" else (" " if x == "-" or y == "-" else ".")
                        for x, y in zip(a, b))
        lines: list[str] = []
        for i in range(0, len(a), width):
            lines.append(f"Query  {i+1:5d} {a[i:i+width]}")
            lines.append(f"              {match[i:i+width]}")
            lines.append(f"Subjct {i+1:5d} {b[i:i+width]}\n")
        return "\n".join(lines)


def _gotoh(a: str, b: str, match: float, mismatch: float,
           gap_open: float, gap_extend: float, mode: str):
    n, m = len(a), len(b)
    NEG = -1e12
    M = [[0.0] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    trace = [[-1] * (m + 1) for _ in range(n + 1)]
    if mode == "global":
        for i in range(1, n + 1):
            X[i][0] = gap_open + (i - 1) * gap_extend
            M[i][0] = X[i][0]; trace[i][0] = 1
        for j in range(1, m + 1):
            Y[0][j] = gap_open + (j - 1) * gap_extend
            M[0][j] = Y[0][j]; trace[0][j] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if a[i-1] == b[j-1] else mismatch
            diag = M[i-1][j-1] + s
            X[i][j] = max(M[i-1][j] + gap_open, X[i-1][j] + gap_extend)
            Y[i][j] = max(M[i][j-1] + gap_open, Y[i][j-1] + gap_extend)
            if mode == "local":
                best = max(0.0, diag, X[i][j], Y[i][j])
                if best == 0.0:
                    M[i][j] = 0.0; trace[i][j] = -1
                elif best == diag:
                    M[i][j] = diag; trace[i][j] = 0
                elif best == X[i][j]:
                    M[i][j] = X[i][j]; trace[i][j] = 1
                else:
                    M[i][j] = Y[i][j]; trace[i][j] = 2
            else:
                best = max(diag, X[i][j], Y[i][j])
                if best == diag:
                    M[i][j] = diag; trace[i][j] = 0
                elif best == X[i][j]:
                    M[i][j] = X[i][j]; trace[i][j] = 1
                else:
                    M[i][j] = Y[i][j]; trace[i][j] = 2
    return M, trace


def _traceback(a, b, trace, mode, start):
    aa: list[str] = []; bb: list[str] = []
    i, j = start
    end_i = i
    while i > 0 and j > 0:
        t = trace[i][j]
        if t == -1: break
        if t == 0:
            aa.append(a[i-1]); bb.append(b[j-1]); i -= 1; j -= 1
        elif t == 1:
            aa.append(a[i-1]); bb.append("-"); i -= 1
        elif t == 2:
            aa.append("-"); bb.append(b[j-1]); j -= 1
        else:
            break
    if mode == "global":
        while i > 0:
            aa.append(a[i-1]); bb.append("-"); i -= 1
        while j > 0:
            aa.append("-"); bb.append(b[j-1]); j -= 1
    return "".join(reversed(aa)), "".join(reversed(bb)), i, j


def align_global(a: str, b: str, *, match: float = 2.0, mismatch: float = -1.0,
                 gap_open: float = -10.0, gap_extend: float = -0.5) -> Alignment:
    """Needleman-Wunsch global alignment."""
    a, b = a.upper(), b.upper()
    if not a or not b:
        raise ValueError("Cannot align empty sequences")
    M, trace = _gotoh(a, b, match, mismatch, gap_open, gap_extend, "global")
    aa, bb, si, sj = _traceback(a, b, trace, "global", (len(a), len(b)))
    matches = sum(1 for x, y in zip(aa, bb) if x == y and x != "-")
    identity = matches / max(1, len(aa)) * 100.0
    return Alignment(aa, bb, M[len(a)][len(b)], 0, len(a), 0, len(b), "global", identity)


def align_local(a: str, b: str, *, match: float = 2.0, mismatch: float = -1.0,
                gap_open: float = -10.0, gap_extend: float = -0.5) -> Alignment:
    """Smith-Waterman local alignment."""
    a, b = a.upper(), b.upper()
    if not a or not b:
        raise ValueError("Cannot align empty sequences")
    M, trace = _gotoh(a, b, match, mismatch, gap_open, gap_extend, "local")
    bi = bj = 0; bs = 0.0
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            if M[i][j] > bs:
                bs = M[i][j]; bi, bj = i, j
    if bs == 0.0:
        return Alignment("", "", 0.0, 0, 0, 0, 0, "local", 0.0)
    aa, bb, si, sj = _traceback(a, b, trace, "local", (bi, bj))
    matches = sum(1 for x, y in zip(aa, bb) if x == y and x != "-")
    identity = matches / max(1, len(aa)) * 100.0
    return Alignment(aa, bb, bs, si, bi, sj, bj, "local", identity)
