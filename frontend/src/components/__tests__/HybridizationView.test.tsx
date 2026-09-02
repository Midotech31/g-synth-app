import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import type { HybridizationResult } from "../../api/client";
import HybridizationView from "../HybridizationView";

const RESULT: HybridizationResult = {
  first: "AATTATGCGT",
  second: "ACGCAT",
  top: "AATTATGCGT",
  marks: "    ||||||",
  bottom: "    TACGCA",
  rows: [{
    start: 0,
    stop: 10,
    top: "AATTATGCGT",
    marks: "    ||||||",
    bottom: "    TACGCA",
    top_start: 1,
    top_end: 10,
    bottom_start: 6,
    bottom_end: 1,
  }],
  width: 10,
  offset: 4,
  overlap_start: 4,
  overlap_end: 10,
  overlap_length: 6,
  paired_bases: 6,
  paired_percent: 100,
  mismatches: 0,
  longest_perfect_run: 6,
  complementarity: "exact",
  predicted_state: "favourable_at_temperature",
  overhangs: [{
    end: "left",
    strand: "first",
    polarity: "5′",
    sequence: "AATT",
    length: 4,
    start: 0,
    end_position: 4,
  }],
  left_end: {
    end: "left",
    strand: "first",
    polarity: "5′",
    sequence: "AATT",
    length: 4,
    start: 0,
    end_position: 4,
  },
  right_end: { end: "right", kind: "blunt", length: 0 },
  alternative_placements: 0,
  tm_c: 32.4,
  tm_margin_c: 7.4,
  delta_h_kcal_mol: -40.1,
  delta_s_cal_mol_k: -110.2,
  analysis_temperature_c: 25,
  conditions: {
    name: "annealing",
    oligo_nM: 50_000,
    na_mM: 50,
    mg_mM: 0,
    dntp_mM: 0,
    summary: "50 µM total strand, 50 mM Na⁺",
  },
  warnings: [],
};

describe("HybridizationView", () => {
  it("shows a cohesive end separately from the paired duplex", () => {
    render(<HybridizationView result={RESULT} detail="simple" />);

    expect(screen.getByText("5′-AATT-3′")).toBeInTheDocument();
    expect(screen.getByText("6 complementary base pairs")).toBeInTheDocument();
    expect(screen.getByText("Blunt")).toBeInTheDocument();
  });

  it("draws the lower strand antiparallel and leaves the overhang unpaired", () => {
    const { container } = render(<HybridizationView result={RESULT} detail="detailed" />);

    const duplex = screen.getByLabelText("Nucleotide-level antiparallel double-strand visualization");
    expect(duplex).toHaveTextContent("5′AATTATGCGT3′");
    expect(duplex).toHaveTextContent(/3′\s*TACGCA5′/);
    expect(container.querySelectorAll(".hybrid-overhang-base")).toHaveLength(4);
    expect(container.querySelectorAll(".hybrid-mismatch-base")).toHaveLength(0);
  });

  it("marks an internal mismatch instead of hiding it as an overhang", () => {
    const mismatched: HybridizationResult = {
      ...RESULT,
      marks: "    ||×|||",
      mismatches: 1,
      complementarity: "partial",
      rows: [{ ...RESULT.rows[0], marks: "    ||×|||" }],
    };

    const { container } = render(<HybridizationView result={mismatched} detail="detailed" />);

    expect(container.querySelectorAll(".hybrid-mismatch-base")).toHaveLength(2);
    expect(screen.getByText(/× marks a mismatch/i)).toBeInTheDocument();
  });
});
