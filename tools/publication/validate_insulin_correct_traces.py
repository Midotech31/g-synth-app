#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import Bio
from Bio.Align import PairwiseAligner

from gsynth_engine.chromatogram import Chromatogram, read_trace, summarise
from gsynth_engine.sequence import reverse_complement
from gsynth_engine.verify import ConsensusReport, VerificationReport, assemble_consensus, verify

TRACE_NAMES = {
    "A": {"forward": "A Forward Seq.ab1", "reverse": "A Reverse Seq.ab1"},
    "B": {"forward": "B Forward Seq.ab1", "reverse": "B Reverse Seq.ab1"},
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text().splitlines()
        if line and not line.startswith(">")
    ).upper()


def report_dict(report: VerificationReport) -> dict:
    return {
        "design_length": report.design_length,
        "coverage_percent": report.coverage,
        "fully_covered": report.fully_covered,
        "is_verified": report.is_verified,
        "gaps_zero_based_half_open": report.gaps,
        "differences": [
            {
                "kind": difference.kind,
                "position_zero_based": difference.position,
                "expected": difference.expected,
                "found": difference.found,
                "quality": difference.quality,
                "confident": difference.confident,
            }
            for difference in report.differences
        ],
        "reads": [
            {
                "name": read.name,
                "reference_start_zero_based": read.start,
                "reference_end_zero_based_half_open": read.end,
                "reverse_complemented": read.reverse_complemented,
                "identity_percent": read.identity,
                "matched_bases": read.matched,
                "trimmed_from_read_start": read.trimmed_start,
                "trimmed_from_read_end": read.trimmed_end,
                "difference_count": len(read.differences),
                "mean_quality": read.mean_quality,
            }
            for read in report.reads
        ],
        "warnings": report.warnings,
    }


def independent_biopython(reference: str, trace: Chromatogram) -> dict:

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    candidates = []
    for reverse, sequence in (
        (False, trace.sequence),
        (True, reverse_complement(trace.sequence)),
    ):
        alignment = aligner.align(reference, sequence)[0]
        counts = alignment.counts()
        aligned_bases = counts.identities + counts.mismatches
        candidates.append((alignment.score, reverse, alignment, counts, aligned_bases))
    _score, reverse, alignment, counts, aligned_bases = max(candidates, key=lambda item: item[0])
    target_blocks = alignment.aligned[0]
    start = int(target_blocks[0][0])
    stop = int(target_blocks[-1][1])
    identity = 100.0 * counts.identities / max(1, aligned_bases + counts.gaps)
    return {
        "reverse_complemented": reverse,
        "reference_start_zero_based": start,
        "reference_end_zero_based_half_open": stop,
        "identities": int(counts.identities),
        "mismatches": int(counts.mismatches),
        "gaps": int(counts.gaps),
        "identity_percent": round(identity, 2),
    }


def union_coverage(length: int, reads: list[dict]) -> float:
    covered: set[int] = set()
    for read in reads:
        covered.update(
            range(
                read["reference_start_zero_based"],
                read["reference_end_zero_based_half_open"],
            )
        )
    return round(100.0 * len(covered) / length, 1)


def consensus_dict(report: ConsensusReport) -> dict:
    return {
        "reference_length": report.reference_length,
        "sequence": report.sequence,
        "coverage_percent": report.coverage,
        "identity_percent": report.identity,
        "fully_covered": report.fully_covered,
        "bidirectional_overlap_percent": report.bidirectional_overlap,
        "bidirectional_overlap_agreement_percent": report.bidirectional_agreement,
        "gaps_zero_based_half_open": report.gaps,
        "differences": [
            {
                "position_zero_based": position.position,
                "reference": position.reference,
                "consensus_call": position.call,
                "supporting_reads": position.supporting_reads,
                "combined_quality": position.combined_quality,
            }
            for position in report.differences
        ],
        "warnings": report.warnings,
    }


def plot_consensus_evidence(
    output: Path,
    chains: dict[str, dict],
    raw_reports: dict,
    consensus_reports: dict[str, ConsensusReport],
) -> None:

    import matplotlib.pyplot as plt

    colours = {"A": "#15803d", "C": "#1d4ed8", "G": "#111827", "T": "#b91c1c"}
    panels = []
    for chain in ("A", "B"):
        length = len(chains[chain]["reference"])
        split = (length + 1) // 2
        panels.extend(((chain, 1, split), (chain, split + 1, length)))

    figure, axes = plt.subplots(4, 1, figsize=(12.0, 9.2))
    for axis, (chain, segment_start, segment_stop) in zip(axes, panels, strict=True):
        reference = chains[chain]["reference"]
        consensus = consensus_reports[chain]
        report = raw_reports[chain]
        axis.set_xlim(segment_start - 0.5, segment_stop + 0.5)
        axis.set_ylim(0.0, 3.0)
        axis.set_yticks([0.65, 1.25, 2.15, 2.65], ["Reverse", "Forward", "Consensus", "Reference"])
        axis.set_xlabel("Reference position (bp)")
        axis.grid(axis="x", color="#e5e7eb", linewidth=0.55)
        for spine in ("top", "right", "left"):
            axis.spines[spine].set_visible(False)

        suffix = ""
        if segment_start == 1:
            suffix = (
                f" · consensus {consensus.coverage:.1f}% covered, "
                f"{consensus.identity:.1f}% identical · F/R overlap "
                f"{consensus.bidirectional_overlap:.1f}% · overlap agreement "
                f"{consensus.bidirectional_agreement:.1f}%"
            )
        axis.set_title(
            f"{chain}-chain · bp {segment_start}–{segment_stop}{suffix}",
            loc="left", fontsize=9.2, fontweight="bold",
        )


        forward = [read for read in report.reads if not read.reverse_complemented]
        reverse = [read for read in report.reads if read.reverse_complemented]
        for fwd in forward:
            for rev in reverse:
                low, high = max(fwd.start, rev.start), min(fwd.end, rev.end)
                if high > low:
                    axis.axvspan(
                        max(segment_start - 0.5, low + 0.5),
                        min(segment_stop + 0.5, high + 0.5),
                        color="#d8f0e5", alpha=0.9, zorder=0,
                    )

        for read in report.reads:
            y = 0.65 if read.reverse_complemented else 1.25
            visible_start = max(segment_start, read.start + 1)
            visible_stop = min(segment_stop, read.end)
            if visible_start <= visible_stop:
                axis.plot(
                    [visible_start, visible_stop], [y, y],
                    color="#0f766e" if not read.reverse_complemented else "#72518d",
                    linewidth=6, solid_capstyle="butt",
                )
                axis.text(
                    (visible_start + visible_stop) / 2, y - 0.17,
                    read.name, ha="center", va="top", fontsize=6.8,
                    color="#334155",
                )

        for position in range(segment_start, segment_stop + 1):
            ref_base = reference[position - 1]
            call = consensus.sequence[position - 1]
            axis.text(
                position, 2.65, ref_base, ha="center", va="center",
                fontsize=7.8, family="DejaVu Sans Mono", fontweight="bold",
                color=colours.get(ref_base, "#64748b"),
            )
            axis.text(
                position, 2.15, call, ha="center", va="center",
                fontsize=7.8, family="DejaVu Sans Mono", fontweight="bold",
                color=colours.get(call, "#64748b"),
            )

    figure.suptitle(
        "G-Synth forward/reverse consensus validation",
        fontsize=14, fontweight="bold", color="#0f172a",
    )
    figure.text(
        0.5, 0.012,
        "Green shading marks positions supported by both read orientations; every overlapping call agreed.\n"
        "Unshaded terminal positions are supported by one oriented read and remain part of the complete assembled consensus.",
        ha="center", va="bottom", fontsize=8.0, color="#475569",
    )
    figure.subplots_adjust(left=0.10, right=0.99, top=0.91, bottom=0.105, hspace=0.72)
    figure.savefig(output, dpi=300, facecolor="white")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True,
                        help="Directory containing the two reference FASTA files")
    parser.add_argument("--trace-dir", type=Path, required=True,
                        help="Directory containing exactly the four approved trace files")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    trace_paths = {
        chain: {
            direction: args.trace_dir / filename
            for direction, filename in directions.items()
        }
        for chain, directions in TRACE_NAMES.items()
    }
    missing = [
        str(path)
        for directions in trace_paths.values()
        for path in directions.values()
        if not path.is_file()
    ]
    if missing:
        raise FileNotFoundError("Missing approved trace file(s): " + ", ".join(missing))

    references = {
        chain: args.root / f"Synthetic Insulin Glargine {chain}-Chain.fasta"
        for chain in ("A", "B")
    }
    evidence = {
        "analysis": {
            "software": "G-Synth",
            "version": "1.0.0",
            "biopython_version": Bio.__version__,
            "primary_consensus_call_cutoff": 0,
            "trace_set_policy": (
                "Only the four author-designated files named A Forward Seq.ab1, "
                "A Reverse Seq.ab1, B Forward Seq.ab1 and B Reverse Seq.ab1 "
                "were admitted to this validation."
            ),
            "interpretation": (
                "Forward and reverse reads are oriented and assembled before "
                "consensus coverage is calculated. Consensus coverage, identity, "
                "bidirectional overlap and overlap agreement are the primary "
                "retrospective endpoints. No post hoc threshold-specific coverage "
                "is reported as an article result."
            ),
        },
        "chains": {},
    }
    figure_chains = {}
    raw_reports = {}
    consensus_reports = {}

    for chain in ("A", "B"):
        reference = read_fasta(references[chain])
        traces = {
            direction: read_trace(path.read_bytes(), name=direction)
            for direction, path in trace_paths[chain].items()
        }
        reads = {name: trace.sequence for name, trace in traces.items()}
        raw = verify(
            reference,
            reads,
            circular=False,
            traces=traces,
            region=(0, len(reference)),
            trim_quality=0,
        )
        raw_consensus = assemble_consensus(
            reference,
            traces,
            circular=False,
            region=(0, len(reference)),
            trim_quality=0,
        )
        independent = {
            name: independent_biopython(reference, trace)
            for name, trace in traces.items()
        }
        evidence["chains"][chain] = {
            "reference": {
                "path": str(references[chain]),
                "sha256": digest(references[chain]),
                "length_bp": len(reference),
                "sequence": reference,
            },
            "trace_files": {
                direction: {
                    "approved_filename": path.name,
                    "path": str(path),
                    "sha256": digest(path),
                    "magic": path.read_bytes()[:4].decode("ascii"),
                    "filename_extension": path.suffix,
                    "summary": summarise(traces[direction]).__dict__,
                }
                for direction, path in trace_paths[chain].items()
            },
            "gsynth_raw_base_calls_q0": report_dict(raw),
            "gsynth_forward_reverse_consensus_raw_q0": consensus_dict(raw_consensus),
            f"biopython_{Bio.__version__.replace('.', '_')}_local_alignment": {
                "reads": independent,
                "combined_reference_coverage_percent": union_coverage(
                    len(reference), list(independent.values())
                ),
            },
        }
        figure_chains[chain] = {"reference": reference, "traces": traces}
        raw_reports[chain] = raw
        consensus_reports[chain] = raw_consensus

    plot_consensus_evidence(
        args.output_dir / "Figure_GSynth_Sanger_Approved_Traces.png",
        figure_chains,
        raw_reports,
        consensus_reports,
    )
    output = args.output_dir / "glargine_approved_trace_validation.json"
    output.write_text(json.dumps(evidence, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
