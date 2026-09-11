#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import io
import json
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

SERVICE_ID = 537
SERVICE_RELEASE = "September 2021"
SERVICE_MODIFIED = "2021-10-01"
BASE_URL = "https://dnahive.fda.gov/dna.cgi"

HOSTS = (
    ("ecoli", "Escherichia coli", "Bacteria", 562),
    ("b_subtilis", "Bacillus subtilis", "Bacteria", 1423),
    ("p_putida", "Pseudomonas putida", "Bacteria", 303),
    ("l_lactis", "Lactococcus lactis", "Bacteria", 1358),
    ("c_glutamicum", "Corynebacterium glutamicum", "Bacteria", 1718),
    ("s_coelicolor", "Streptomyces coelicolor", "Bacteria", 100226),
    ("s_cerevisiae", "Saccharomyces cerevisiae", "Yeasts", 4932),
    ("k_phaffii", "Komagataella phaffii (Pichia pastoris)", "Yeasts", 4922),
    ("k_lactis", "Kluyveromyces lactis", "Yeasts", 28985),
    ("y_lipolytica", "Yarrowia lipolytica", "Yeasts", 4952),
    ("h_sapiens", "Homo sapiens (species-level HEK proxy)", "Mammalian", 9606),
    ("c_griseus", "Cricetulus griseus (species-level CHO proxy)", "Mammalian", 10029),
    ("s_frugiperda", "Spodoptera frugiperda (Sf9/Sf21 proxy)", "Insect", 7108),
    ("d_melanogaster", "Drosophila melanogaster (S2 proxy)", "Insect", 7227),
    ("n_benthamiana", "Nicotiana benthamiana", "Plant", 4100),
)


def _fetch(taxon_id: int, dataset: str) -> tuple[dict[str, str], dict[str, int]]:
    parameters = {
        "cmd": "ionTaxidCollapseExt",
        "svcType": "svc-codon-usage",
        "objId": str(SERVICE_ID),
        "fileSource": f"{dataset}_species.tsv",
        "plen": "3",
        "taxid": str(taxon_id),
        "filterInColName": '["Organelle"]',
        "filterIn": '["genomic"]',
        "searchDeep": "true",
        "raw": "1",
    }
    request = Request(
        f"{BASE_URL}?{urlencode(parameters)}",
        headers={"User-Agent": "G-Synth codon-table updater/1.0"},
    )
    with urlopen(request, timeout=60) as response:
        text = response.read().decode("utf-8-sig")

    metadata: dict[str, str] = {}
    counts: dict[str, int] = {}
    for key, value, *_ in csv.reader(io.StringIO(text)):
        if key == "id":
            continue
        if len(key) == 3 and set(key) <= set("ACGT"):
            counts[key] = int(value)
        else:
            metadata[key] = value
    return metadata, counts


def _host_record(name: str, category: str, taxon_id: int) -> dict[str, object]:
    selected: tuple[str, dict[str, str], dict[str, int]] | None = None
    for dataset in ("Refseq", "genbank"):
        metadata, counts = _fetch(taxon_id, dataset)
        if int(metadata.get("#codon", "0")) > 0 and len(counts) == 64:
            returned_taxon = int(metadata.get("taxid", "-1"))
            if returned_taxon != taxon_id:
                raise RuntimeError(
                    f"Taxon mismatch: requested {taxon_id}, received {returned_taxon}"
                )
            if any(count < 0 for count in counts.values()):
                raise RuntimeError(f"Negative codon count for NCBI taxon {taxon_id}")
            selected = dataset, metadata, counts
            break
    if selected is None:
        raise RuntimeError(f"No complete HIVE-CUTs table for NCBI taxon {taxon_id}")

    dataset, metadata, counts = selected
    codon_count = int(metadata["#codon"])
    if sum(counts.values()) != codon_count:
        raise RuntimeError(
            f"Codon-count mismatch for taxon {taxon_id}: "
            f"metadata={codon_count}, rows={sum(counts.values())}"
        )
    return {
        "name": name,
        "category": category,
        "taxon_id": taxon_id,
        "dataset": "RefSeq" if dataset == "Refseq" else "GenBank",
        "data_scope": "genomic species aggregate including descendant taxa",
        "coding_sequences": int(metadata["#CDS"]),
        "codon_count": codon_count,
        "gc_percent": float(metadata["GC%"]),
        "counts": dict(sorted(counts.items())),
    }


def build_snapshot(retrieved: str) -> dict[str, object]:
    return {
        "schema_version": 1,
        "source": {
            "name": "FDA HIVE-CUTs / CoCoPUTs",
            "service_id": SERVICE_ID,
            "release": SERVICE_RELEASE,
            "modified": SERVICE_MODIFIED,
            "retrieved": retrieved,
            "url": "https://dnahive.fda.gov/dna.cgi?cmd=cuts_main",
            "method": (
                "Genomic species aggregate with descendant taxa; RefSeq preferred, "
                "GenBank fallback only when RefSeq was unavailable"
            ),
        },
        "hosts": {
            key: _host_record(name, category, taxon_id)
            for key, name, category, taxon_id in HOSTS
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("gsynth_engine/data/codon_usage_hive_2021.json"),
    )
    parser.add_argument("--retrieved", default=date.today().isoformat())
    arguments = parser.parse_args()
    snapshot = build_snapshot(arguments.retrieved)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(snapshot, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
