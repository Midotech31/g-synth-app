"""
Remote sequence retrieval.

NCBI Entrez:
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
UniProt:
    https://rest.uniprot.org/uniprotkb/search?query=...&format=fasta

Fixes the bug in G-Synth 2.x where the UniProt URL was malformed
(missing query parameter encoding, used the obsolete `sequence:` filter).
"""

from __future__ import annotations

import urllib.parse
import urllib.request
from typing import Literal

from gsynth_core.io.fasta import parse_fasta
from gsynth_core.io.genbank import parse_genbank
from gsynth_core.io.records import SeqRecord


def _http_get(url: str, *, timeout: int = 30) -> str:
    """Fetch a URL with a polite User-Agent. Raises on HTTP errors."""
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "G-Synth/3.0 (https://github.com/gsynth-team/gsynth; merzoug.mohamed@essbo.dz)"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def fetch_ncbi(
    accession: str,
    *,
    db: Literal["nuccore", "nucleotide", "protein", "gene"] = "nuccore",
    rettype: Literal["fasta", "gb", "genbank"] = "gb",
    email: str | None = None,
    api_key: str | None = None,
    timeout: int = 30,
) -> SeqRecord:
    """
    Fetch a single record from NCBI Entrez.

    Args:
        accession: NCBI accession (e.g. "NM_000546" for TP53 mRNA, "U13369"
            for E. coli rrnB, "P04637" → use fetch_uniprot for proteins).
        db: NCBI database name. `nuccore` covers most genomic / mRNA records.
        rettype: `gb` returns GenBank (with features); `fasta` returns raw seq.
        email: optional contact email — strongly recommended by NCBI ToS.
        api_key: optional NCBI API key for higher rate limits.

    Returns:
        :class:`SeqRecord` with features populated when GenBank format is used.

    Raises:
        ValueError: if the accession is empty or NCBI returns no record.
        urllib.error.URLError: for network failures.
    """
    if not accession or not accession.strip():
        raise ValueError("Accession cannot be empty")
    if rettype == "genbank":
        rettype = "gb"

    params = {
        "db": db,
        "id": accession.strip(),
        "rettype": rettype,
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
    text = _http_get(url, timeout=timeout)
    if not text.strip():
        raise ValueError(f"NCBI returned empty record for {accession!r}")

    if rettype == "gb":
        records = parse_genbank(text)
    else:
        records = parse_fasta(text)

    if not records:
        raise ValueError(f"Failed to parse NCBI response for {accession!r}")
    return records[0]


def fetch_uniprot(
    query: str,
    *,
    timeout: int = 30,
    format: Literal["fasta", "txt"] = "fasta",
) -> SeqRecord:
    """
    Fetch the top-matching UniProt entry.

    Args:
        query: an accession (e.g. "P04637") OR a free-text query
            (e.g. "p53 human"). Accessions of length 6 or 10 are looked up
            via the entry endpoint directly; otherwise we use the search API.
        format: "fasta" (default) or "txt" (UniProt flat-file).

    Returns:
        :class:`SeqRecord` for the top hit.
    """
    q = query.strip()
    if not q:
        raise ValueError("Empty UniProt query")

    looks_like_accession = (
        len(q) in (6, 10) and q[0].isalpha() and q[1:].replace("-", "").isalnum() and " " not in q
    )

    if looks_like_accession:
        # Direct entry endpoint — no malformed query strings, no `sequence:` filter
        url = f"https://rest.uniprot.org/uniprotkb/{urllib.parse.quote(q)}.{format}"
    else:
        params = {
            "query": q,
            "format": format,
            "size": "1",
        }
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)

    text = _http_get(url, timeout=timeout)
    if not text.strip():
        raise ValueError(f"UniProt returned no record for {q!r}")

    if format == "fasta":
        records = parse_fasta(text)
        if not records:
            raise ValueError(f"Could not parse UniProt response for {q!r}")
        rec = records[0]
        # Mark as protein
        rec.molecule_type = "protein"
        return rec

    # Flat-file: trivial parse for now (just take the sequence and name)
    name = q
    seq_lines: list[str] = []
    in_sq = False
    for line in text.splitlines():
        if line.startswith("ID"):
            parts = line.split()
            if len(parts) > 1:
                name = parts[1]
        elif line.startswith("SQ"):
            in_sq = True
            continue
        elif line.startswith("//"):
            in_sq = False
        elif in_sq:
            seq_lines.append("".join(c for c in line if c.isalpha()))
    return SeqRecord(sequence="".join(seq_lines), name=name, molecule_type="protein")
