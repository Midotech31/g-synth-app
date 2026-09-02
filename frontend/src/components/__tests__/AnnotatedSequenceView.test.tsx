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
