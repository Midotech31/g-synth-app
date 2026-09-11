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


const FRAGMENT_COLOURS = ["#0e6e77", "#c97634", "#6a4c93", "#3f7a52"];

export function fragmentColour(index: number): string {
  return FRAGMENT_COLOURS[Math.max(0, index) % FRAGMENT_COLOURS.length];
}
