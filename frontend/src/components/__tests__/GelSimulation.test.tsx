import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import type { GelSimulation as GelData } from "../../api/client";
import GelSimulation, { gelBandPosition } from "../GelSimulation";

const simulation: GelData = {
  title: "Predicted PCR product",
  prediction_only: true,
  notice: "In-silico size prediction; confirm experimentally.",
  recommended_ladder: "100-bp",
  ladders: [
    { key: "100-bp", name: "100 bp reference ladder", bands: [100, 500, 1000, 1500], range: "100–1,500 bp" },
    { key: "1-kb", name: "1 kb reference ladder", bands: [500, 1000, 5000, 10000], range: "500–10,000 bp" },
  ],
  lanes: [
    { name: "PCR", description: "Expected amplicon", bands: [{ size_bp: 740, label: "740 bp amplicon" }] },
    { name: "NTC", description: "No-template control", bands: [] },
  ],
};

describe("predicted gel visualization", () => {
  it("places larger DNA above smaller DNA on a logarithmic gel", () => {
    expect(gelBandPosition(5000, 100, 10000)).toBeLessThan(gelBandPosition(500, 100, 10000));
  });

  it("labels the plot as in silico and lists exact marker sizes", () => {
    render(<GelSimulation simulation={simulation} />);
    expect(screen.getByText("In silico prediction")).toBeInTheDocument();
    expect(screen.getByText(/100, 500, 1,000, 1,500 bp/i)).toBeInTheDocument();
    expect(screen.getByText("No specific band expected")).toBeInTheDocument();
  });

  it("lets the user choose a different size marker", () => {
    render(<GelSimulation simulation={simulation} />);
    fireEvent.change(screen.getByLabelText("DNA size marker"), { target: { value: "1-kb" } });
    expect(screen.getByText(/500, 1,000, 5,000, 10,000 bp/i)).toBeInTheDocument();
  });
});
