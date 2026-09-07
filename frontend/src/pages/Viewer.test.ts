import { describe, expect, it } from "vitest";

import type { Annotation } from "../api/client";
import { annotationAt, sequenceForAnnotation } from "./Viewer";

const wrapped: Annotation = {
  name: "origin feature",
  type: "CDS",
  start: 8,
  end: 14,
  direction: 1,
  color: "#000000",
};

describe("wrapped circular annotations", () => {
  it("extracts the bases on both sides of the origin", () => {
    expect(sequenceForAnnotation("AACCGGTTAA", wrapped)).toBe("AAAACC");
  });

  it("extracts a wrapped reverse-strand feature in biological orientation", () => {
    expect(sequenceForAnnotation("AACCGGTTAA", { ...wrapped, direction: -1 }))
      .toBe("GGTTTT");
  });

  it("finds the feature when SeqViz reports a click after the origin", () => {
    expect(annotationAt([wrapped], 1, 2, 10)).toBe(wrapped);
  });
});


it("does not select a wrapped feature for a span extending past its end", () => {
  expect(annotationAt([wrapped], 1, 7, 10)).toBeNull();
});
it("complements ambiguous IUPAC bases on a reverse-strand feature", () => {
  expect(sequenceForAnnotation("ARYKBDHV", { ...wrapped, start: 0, end: 8, direction: -1 })).toBe("BDHVMRYT");
});
