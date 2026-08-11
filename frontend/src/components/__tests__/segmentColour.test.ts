/**
 * The construct map and the hybridisation view colour the same cassette from
 * this one table. A part that is teal in one drawing and grey in the other
 * reads as two different parts, and the reader trusts the drawing over the
 * sequence. These checks pin the matching rule the engine's segment names
 * depend on: the engine names segments in context ("NdeI overhang", "left
 * linker"), never by the bare key.
 */
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
    // The overhangs are the pair a reader checks against the enzyme table;
    // they have to be findable at a glance on either side.
    expect(segmentColour("BamHI overhang")).toBe(segmentColour("EcoRI overhang"));
  });

  it("falls back to a real colour for a name it does not know", () => {
    // An undefined fill paints black in SVG, which looks like a part of its
    // own rather than like one nobody has classified.
    expect(segmentColour("unnamed spacer")).toMatch(/^#[0-9a-f]{6}$/);
    expect(segmentColour("")).toMatch(/^#[0-9a-f]{6}$/);
  });
});

describe("fragmentColour", () => {
  it("gives the first four fragments four different colours", () => {
    // With two, F1 and F3 come out the same and a reader counting fragments
    // off the drawing has to count rather than see.
    const shades = [0, 1, 2, 3].map(fragmentColour);
    expect(new Set(shades).size).toBe(4);
  });

  it("repeats only after four, so neighbours never share a colour", () => {
    expect(fragmentColour(4)).toBe(fragmentColour(0));
    expect(fragmentColour(5)).not.toBe(fragmentColour(4));
  });

  it("stays inside the palette for an index below zero", () => {
    // A negative modulo would index off the end and return undefined.
    expect(fragmentColour(-1)).toBe(fragmentColour(0));
  });
});
