import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import type { Annotation } from "../../api/client";
import AnnotatedSequenceView, {
  chooseAnnotationWindow,
  readableTextColour,
  visibleFeaturesForRow,
} from "../AnnotatedSequenceView";

const annotations: Annotation[] = [
  {
    name: "T7 promoter", type: "promoter", start: 4, end: 18,
    direction: 1, color: "#c97634",
  },
  {
    name: "Insulin cassette", type: "CDS", start: 21, end: 42,
    direction: 1, color: "#0e6e77", translation_start: 21, translation_end: 42,
  },
  {
    name: "6×His tag", type: "misc_feature", start: 24, end: 36,
    direction: 1, color: "#6a4c93",
  },
];

describe("coordinate-level annotation view", () => {
  it("selects an exact range across sequence rows and passes it to the editor", () => {
    const onAnnotateRange = vi.fn();
    render(<AnnotatedSequenceView sequence={"ACGT".repeat(40)} annotations={[]} selected={null}
      onSelect={vi.fn()} onAnnotateRange={onAnnotateRange} />);
    fireEvent.keyDown(screen.getByRole("button", { name: "Base 59: G" }), { key: "Enter" });
    fireEvent.keyDown(screen.getByRole("button", { name: "Base 63: G" }), { key: "Enter", shiftKey: true });
    expect(screen.getByText("59–63 · 5 nt selected")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "Annotate selection" }));
    expect(onAnnotateRange).toHaveBeenCalledWith({ start: 58, end: 63 });
    fireEvent.click(screen.getByRole("button", { name: "Clear range" }));
    expect(screen.getByRole("button", { name: "Annotate selection" })).toBeDisabled();
  });

  it("normalizes a range selected across the circular origin", () => {
    const onAnnotateRange = vi.fn();
    render(<AnnotatedSequenceView sequence={"A".repeat(1000)} circular
      annotations={[{ name: "join", type: "misc_feature", start: 950, end: 1010, direction: 1, color: "#3F7A52" }]}
      selected={null} preferredName="join" onSelect={vi.fn()} onAnnotateRange={onAnnotateRange} />);
    fireEvent.keyDown(screen.getByRole("button", { name: "Base 999: A" }), { key: "Enter" });
    fireEvent.keyDown(screen.getByRole("button", { name: "Base 3: A" }), { key: "Enter", shiftKey: true });
    fireEvent.click(screen.getByRole("button", { name: "Annotate selection" }));
    expect(onAnnotateRange).toHaveBeenCalledWith({ start: 998, end: 1003 });
  });

  it("keeps the preferred construct and nearby regulatory features in view", () => {
    expect(chooseAnnotationWindow(annotations, 120, null, "Insulin cassette"))
      .toEqual({ start: 0, end: 120 });
  });

  it("keeps downstream features visible when a circular cassette crosses the origin", () => {
    const circularFeatures: Annotation[] = [
      { name: "cassette", type: "CDS", start: 930, end: 1000, direction: 1, color: "#0e6e77" },
      { name: "terminator", type: "terminator", start: 35, end: 75, direction: 1, color: "#9e3d3d" },
    ];
    const window = chooseAnnotationWindow(circularFeatures, 1000, null, "cassette", true);
    expect(window.start).toBeLessThan(930);
    expect(window.end).toBeGreaterThan(1075);
  });

  it("assigns overlapping features to separate lanes", () => {
    const visible = visibleFeaturesForRow(annotations, 0, 60);
    const cassette = visible.find((feature) => feature.annotation.name === "Insulin cassette");
    const tag = visible.find((feature) => feature.annotation.name === "6×His tag");
    expect(cassette?.lane).not.toBe(tag?.lane);
  });

  it("chooses readable text for both light and dark feature colours", () => {
    expect(readableTextColour("#f7e8bc")).toBe("#0b1f3b");
    expect(readableTextColour("#0e6e77")).toBe("#ffffff");
  });

  it("renders coordinates, feature tracks and an in-frame translation", () => {
    const onSelect = vi.fn();
    const sequence = `${"A".repeat(21)}ATGGGTTCCTAA${"C".repeat(87)}`;
    render(
      <AnnotatedSequenceView
        sequence={sequence}
        annotations={annotations}
        selected={null}
        preferredName="Insulin cassette"
        onSelect={onSelect}
      />,
    );

    expect(screen.getByRole("heading", { name: "Annotated sequence" })).toBeInTheDocument();
    expect(screen.getByText("Met")).toBeInTheDocument();
    expect(screen.getByText("Gly")).toBeInTheDocument();
    expect(screen.getByText("*")).toBeInTheDocument();

    fireEvent.click(screen.getByRole("button", { name: /6×His tag, misc_feature/ }));
    expect(onSelect).toHaveBeenCalledWith(annotations[2]);
  });
});


describe("whole-record navigation", () => {
  it("reaches the final bases of a 200 kb record without rendering every base", () => {
    const { container } = render(<AnnotatedSequenceView sequence={"A".repeat(200000)} annotations={[]} selected={null} circular onSelect={vi.fn()} />);
    fireEvent.click(screen.getByRole("button", { name: "Whole plasmid" }));
    expect(screen.getByText(/Showing bases 1–200,000/)).toBeInTheDocument();
    expect(container.querySelectorAll(".annotation-base").length).toBeLessThan(1500);
    fireEvent.change(screen.getByLabelText("Go to base"), { target: { value: "200000" } });
    fireEvent.click(screen.getByRole("button", { name: "Go" }));
    expect(screen.getByTitle("200,000: A")).toBeInTheDocument();
    expect(screen.queryByTitle("200,001: A")).not.toBeInTheDocument();
  });

  it("rejects invalid coordinates and keeps sequence navigation usable", () => {
    render(<AnnotatedSequenceView sequence={"A".repeat(1500)} annotations={[]} selected={null} onSelect={vi.fn()} />);
    fireEvent.change(screen.getByLabelText("Go to base"), { target: { value: "1600" } });
    fireEvent.click(screen.getByRole("button", { name: "Go" }));
    expect(screen.getByRole("alert")).toHaveTextContent("1 to 1,500");
    fireEvent.click(screen.getByRole("button", { name: "Next region" }));
    expect(screen.getByText(/Showing bases 301–900/)).toBeInTheDocument();
    fireEvent.change(screen.getByLabelText("Bases per row"), { target: { value: "30" } });
    expect(screen.getByRole("combobox")).toHaveValue("30");
  });

  it("does not repeat a small circular molecule or invent strandedness", () => {
    const feature = { name: "unstranded", type: "misc_feature", start: 0, end: 6, direction: 0, color: "#ffffff" };
    render(<AnnotatedSequenceView sequence="AAACCCGGG" annotations={[feature]} selected={null} circular onSelect={vi.fn()} />);
    expect(screen.getByRole("button", { name: /unstranded strand/ })).toBeInTheDocument();
    expect(screen.getByText(/Showing bases 1–9/)).toBeInTheDocument();
  });

  it("keeps reverse-strand codons aligned when a row splits a codon", () => {
    const feature = { name: "reverse CDS", type: "CDS", start: 53, end: 62, direction: -1, color: "#0e6e77" };
    render(<AnnotatedSequenceView sequence={"A".repeat(53) + "TTAGGGCAT"} annotations={[feature]} selected={null} onSelect={vi.fn()} />);
    expect(screen.getAllByTitle("Residue 1: Met (ATG)")).toHaveLength(2);
    expect(screen.getByTitle("Residue 2: Pro (CCC)")).toBeInTheDocument();
  });
});
