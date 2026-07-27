/**
 * One colour per part of the cassette, shared by every view that draws it.
 *
 * The construct map and the hybridisation view have to agree: a linker that
 * is grey in one and teal in the other is worse than no colour at all.
 * Matching is by substring because the engine names segments in context —
 * "NdeI overhang", "left linker", "Thrombin site".
 */
export const SEGMENT_COLOURS: Record<string, string> = {
  overhang: "#c97634",
  "start codon": "#9e3d3d",
  linker: "#78889b",
  "6×His tag": "#0e6e77",
  site: "#6a4c93",
  insert: "#3f7a52",
};

const FALLBACK = "#78889b";

export function segmentColour(name: string): string {
  const key = Object.keys(SEGMENT_COLOURS).find((candidate) =>
    name.toLowerCase().includes(candidate.toLowerCase()),
  );
  return key ? SEGMENT_COLOURS[key] : FALLBACK;
}

/** Rotating colours for fragment-by-fragment shading.
 *
 * Four rather than two: with two, F1 and F3 come out identical, and a reader
 * counting fragments off the drawing has to count rather than see.
 */
const FRAGMENT_COLOURS = ["#0e6e77", "#c97634", "#6a4c93", "#3f7a52"];

export function fragmentColour(index: number): string {
  return FRAGMENT_COLOURS[Math.max(0, index) % FRAGMENT_COLOURS.length];
}
