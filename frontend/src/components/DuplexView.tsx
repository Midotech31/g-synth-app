import { useMemo, useState } from "react";

import type { Duplex, DuplexSpan } from "../api/client";
import { fragmentColour, segmentColour } from "./segmentColour";


function colourIndex(width: number, segments: DuplexSpan[]): (string | null)[] {
  const colours: (string | null)[] = new Array(width).fill(null);
  for (const span of segments) {
    const colour = segmentColour(span.name);
    for (let i = span.start; i < span.end && i < width; i += 1) colours[i] = colour;
  }
  return colours;
}


function fragmentIndex(width: number, spans: DuplexSpan[]): number[] {
  const owner: number[] = new Array(width).fill(-1);
  spans.forEach((span, index) => {
    for (let i = span.start; i < span.end && i < width; i += 1) owner[i] = index;
  });
  return owner;
}

type Props = {
  duplex: Duplex;

  width?: number;
};

export default function DuplexView({ duplex, width = 60 }: Props) {
  const [colourBy, setColourBy] = useState<"cassette" | "fragment">("cassette");

  const colours = useMemo(
    () => colourIndex(duplex.width, duplex.segments),
    [duplex.width, duplex.segments],
  );
  const topOwner = useMemo(
    () => fragmentIndex(duplex.width, duplex.top_fragments),
    [duplex.width, duplex.top_fragments],
  );
  const bottomOwner = useMemo(
    () => fragmentIndex(duplex.width, duplex.bottom_fragments),
    [duplex.width, duplex.bottom_fragments],
  );


  const topCuts = useMemo(() => new Set(duplex.junctions), [duplex.junctions]);
  const bottomCuts = useMemo(
    () => new Set(duplex.bottom_fragments.slice(0, -1).map((span) => span.end)),
    [duplex.bottom_fragments],
  );


  const overhangColumns = useMemo(() => {
    const columns = new Set<number>();
    const bottom = [...bottomCuts].sort((a, b) => a - b);
    for (const cut of duplex.junctions) {
      const partner = bottom.find((position) => position > cut);
      if (partner === undefined) continue;
      for (let i = cut; i < partner; i += 1) columns.add(i);
    }
    return columns;
  }, [duplex.junctions, bottomCuts]);

  const rows = useMemo(() => {
    const out: { start: number; topFrom: number | null; bottomFrom: number | null }[] = [];
    let topSeen = 0;
    let bottomSeen = 0;

    for (let start = 0; start < duplex.width; start += width) {
      const stop = Math.min(start + width, duplex.width);
      const topBases = countBases(duplex.top, start, stop);
      const bottomBases = countBases(duplex.bottom, start, stop);

      out.push({
        start,
        topFrom: topBases > 0 ? topSeen + 1 : null,
        bottomFrom: bottomBases > 0 ? bottomSeen + 1 : null,
      });
      topSeen += topBases;
      bottomSeen += bottomBases;
    }
    return out;
  }, [duplex, width]);

  const numberWidth = String(duplex.width).length;

  function base(strand: string, owner: number[], column: number, kind: "top" | "bottom") {
    const character = strand[column];
    if (character === " ") return <span key={column} className="dx-gap"> </span>;

    const paired = duplex.pairs[column] === "|";
    const mismatch = duplex.pairs[column] === "x";
    const which = owner[column];

    const style =
      colourBy === "cassette"
        ? { color: colours[column] ?? "var(--ink-soft)" }
        : { color: fragmentColour(which) };

    const classes = ["dx-base"];
    if (!paired) classes.push("dx-single");
    if (mismatch) classes.push("dx-mismatch");
    if ((kind === "top" ? topCuts : bottomCuts).has(column)) classes.push("dx-cut");

    return (
      <span key={column} className={classes.join(" ")} style={style}>
        {character}
      </span>
    );
  }

  return (
    <div className="duplex">
      <div className="duplex-legend">
        <span className="label">Colour by</span>
        <div className="seg-toggle">
          <button
            type="button"
            className={colourBy === "cassette" ? "on" : ""}
            onClick={() => setColourBy("cassette")}
          >
            Cassette
          </button>
          <button
            type="button"
            className={colourBy === "fragment" ? "on" : ""}
            onClick={() => setColourBy("fragment")}
          >
            Fragment
          </button>
        </div>
        <span className="grow" />
        {colourBy === "cassette" ? (
          <span className="keys">
            {dedupe(duplex.segments.map((s) => s.name)).map((name) => (
              <span key={name} className="key">
                <i style={{ background: segmentColour(name) }} />
                {name}
              </span>
            ))}
          </span>
        ) : (
          <span className="keys">
            {duplex.top_fragments.map((span, index) => (
              <span key={span.name} className="key">
                <i style={{ background: fragmentColour(index) }} />
                {span.name}
              </span>
            ))}
          </span>
        )}
      </div>

      <div className="duplex-scroll">
        {rows.map((row) => {
          const stop = Math.min(row.start + width, duplex.width);
          const columns = range(row.start, stop);
          return (
            <div className="duplex-row" key={row.start}>
              <span className="dx-num" style={{ minWidth: `${numberWidth}ch` }}>
                {row.topFrom ?? ""}
              </span>
              <span className="dx-end">5&apos;</span>
              <span className="dx-seq">
                {columns.map((c) => base(duplex.top, topOwner, c, "top"))}
              </span>

              <span className="dx-num" style={{ minWidth: `${numberWidth}ch` }} />
              <span className="dx-end" />
              <span className="dx-seq dx-ticks">
                {columns.map((c) => (
                  <span key={c} className={overhangColumns.has(c) ? "dx-stagger" : ""}>
                    {duplex.pairs[c] === "|" ? "|" : " "}
                  </span>
                ))}
              </span>

              <span className="dx-num" style={{ minWidth: `${numberWidth}ch` }}>
                {row.bottomFrom ?? ""}
              </span>
              <span className="dx-end">3&apos;</span>
              <span className="dx-seq">
                {columns.map((c) => base(duplex.bottom, bottomOwner, c, "bottom"))}
              </span>
            </div>
          );
        })}
      </div>

      <p className="note duplex-foot">
        Top strand 5&apos;→3&apos;, bottom strand 3&apos;→5&apos;. Faded bases are
        single-stranded: the sticky ends the two cloning enzymes leave —{" "}
        <b>{duplex.left_overhang || "blunt"}</b> on the left,{" "}
        <b>{duplex.right_overhang || "blunt"}</b> on the right.
        {duplex.junctions.length > 0 && (
          <>
            {" "}
            The {duplex.junctions.length} pair
            {duplex.junctions.length === 1 ? "" : "s"} of vertical rules mark where
            each strand is cut; the offset between them is the overhang that holds
            the next fragment in place.
          </>
        )}
      </p>
    </div>
  );
}

function countBases(strand: string, start: number, stop: number): number {
  let count = 0;
  for (let i = start; i < stop; i += 1) if (strand[i] !== " ") count += 1;
  return count;
}

function range(start: number, stop: number): number[] {
  return Array.from({ length: stop - start }, (_, i) => start + i);
}

function dedupe(names: string[]): string[] {
  return [...new Set(names)];
}
