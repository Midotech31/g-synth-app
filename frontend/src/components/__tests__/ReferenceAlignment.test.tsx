import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import type { TraceTrack, VerifyReport } from "../../api/client";
import ReferenceAlignment from "../ReferenceAlignment";

const reference = "AACCGGTTAACCGGTT";

function track(overrides: Partial<TraceTrack> = {}): TraceTrack {
  return {
    read: "forward.ab1",
    reference_start: 2,
    reference_end: 10,
    reverse_complemented: false,
    sequence: "CCGGT TAA".replace(" ", ""),
    qualities: [35, 34, 33, 32, 18, 31, 30, 29],
    peaks: [2, 5, 8, 11, 14, 17, 20, 23],
    sample_count: 26,
    traces: {
      A: [0, 1, 3, 8, 3, 1, 0],
      C: [0, 4, 9, 4, 1, 0, 0],
      G: [0, 1, 2, 7, 2, 1, 0],
      T: [0, 0, 1, 5, 10, 4, 0],
    },
    ...overrides,
  };
}

function report(tracks: TraceTrack[]): VerifyReport {
  return {
    design_length: reference.length,
    coverage: 75,
    gaps: [[0, 2]],
    fully_covered: false,
    is_verified: false,
    region_start: 0,
    region_end: reference.length,
    differences: [],
    reads: [],
    warnings: [],
    trace_tracks: tracks,
  };
}

describe("ReferenceAlignment", () => {
  it("renders reference, consensus, strand labels, and chromatogram evidence", () => {
    const reverse = track({
      read: "reverse.ab1",
      reference_start: 6,
      reference_end: 14,
      reverse_complemented: true,
      sequence: "TTAACCGG",
    });

    const { container } = render(
      <ReferenceAlignment reference={reference} report={report([track(), reverse])} />,
    );

    expect(screen.getByRole("heading", { name: "Reference-aligned chromatograms" })).toBeInTheDocument();
    expect(screen.getByText("Consensus")).toBeInTheDocument();
    expect(screen.getByText("Reference")).toBeInTheDocument();
    expect(screen.getByText(/FWD/)).toBeInTheDocument();
    expect(screen.getByText(/REV/)).toBeInTheDocument();
    expect(screen.getByRole("img", { name: /forward\.ab1, forward chromatogram/ })).toBeInTheDocument();
    expect(screen.getByRole("img", { name: /reverse\.ab1, reverse chromatogram/ })).toBeInTheDocument();
    expect(container.querySelectorAll("polyline")).toHaveLength(8);
    expect(container.querySelector(".trace-quality-low")).toBeInTheDocument();
  });

  it("exposes keyboard-scrollable evidence and deterministic zoom controls", () => {
    render(<ReferenceAlignment reference={reference} report={report([track()])} />);

    const scroller = screen.getByLabelText("Scrollable sequencing alignment");
    expect(scroller).toHaveAttribute("tabindex", "0");

    const standard = screen.getByRole("button", { name: "Standard" });
    const compact = screen.getByRole("button", { name: "Compact" });
    expect(standard).toHaveAttribute("aria-pressed", "true");
    fireEvent.click(compact);
    expect(compact).toHaveAttribute("aria-pressed", "true");
    expect(standard).toHaveAttribute("aria-pressed", "false");
  });

  it("renders nothing when no admitted trace evidence is available", () => {
    const { container } = render(
      <ReferenceAlignment reference={reference} report={report([])} />,
    );
    expect(container).toBeEmptyDOMElement();
  });
});
