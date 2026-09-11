import { describe, expect, it } from "vitest";

import { SEGMENT_COLOURS, fragmentColour, segmentColour } from "../segmentColour";

describe("segmentColour", () => {
  it("colours a segment named in context, not only the bare key", () => {
    expect(segmentColour("NdeI overhang")).toBe(SEGMENT_COLOURS.overhang);
    expect(segmentColour("XhoI overhang")).toBe(SEGMENT_COLOURS.overhang);
    expect(segmentColour("left linker")).toBe(SEGMENT_COLOURS.linker);
  });

  it("ignores case, because the engine capitalises enzyme and site names", () => {
    expect(segmentColour("Thrombin Site")).toBe(SEGMENT_COLOURS.site);
    expect(segmentColour("START CODON")).toBe(SEGMENT_COLOURS["start codon"]);
  });

  it("gives the two ends of the cassette the same colour at both ends", () => {


    expect(segmentColour("BamHI overhang")).toBe(segmentColour("EcoRI overhang"));
  });

  it("falls back to a real colour for a name it does not know", () => {

    expect(segmentColour("unnamed spacer")).toMatch(/^#[0-9a-f]{6}$/);
    expect(segmentColour("")).toMatch(/^#[0-9a-f]{6}$/);
  });
});

describe("fragmentColour", () => {
  it("gives the first four fragments four different colours", () => {


    const shades = [0, 1, 2, 3].map(fragmentColour);
    expect(new Set(shades).size).toBe(4);
  });

  it("repeats only after four, so neighbours never share a colour", () => {
    expect(fragmentColour(4)).toBe(fragmentColour(0));
    expect(fragmentColour(5)).not.toBe(fragmentColour(4));
  });

  it("stays inside the palette for an index below zero", () => {

    expect(fragmentColour(-1)).toBe(fragmentColour(0));
  });
});
