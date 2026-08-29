import { useEffect, useMemo, useRef, useState } from "react";

import type { TraceTrack, VerifyReport } from "../api/client";

const TRACE_COLOUR: Record<string, string> = {
  A: "#2f8a52",
  C: "#2869b0",
  G: "#20262e",
  T: "#bd3f3f",
};

const BASE_CLASS: Record<string, string> = {
  A: "base-a", C: "base-c", G: "base-g", T: "base-t", N: "base-n",
};

type Props = {
  reference: string;
  report: VerifyReport;
};

function consensus(reference: string, tracks: TraceTrack[], start: number, end: number) {
  const result: { base: string; covered: boolean; agrees: boolean }[] = [];
  for (let position = start; position < end; position += 1) {
    const votes = new Map<string, number>();
    for (const track of tracks) {
      const index = position - track.reference_start;
      if (index < 0 || index >= track.sequence.length) continue;
      const base = track.sequence[index] ?? "N";
      votes.set(base, (votes.get(base) ?? 0) + Math.max(1, track.qualities[index] ?? 1));
    }
    const called = [...votes].sort((a, b) => b[1] - a[1])[0]?.[0] ?? "·";
    result.push({
      base: called,
      covered: votes.size > 0,
      agrees: !votes.size || called === reference[position],
    });
  }
  return result;
}

function BaseStrip({
  sequence, start, cell, qualities, reference,
}: {
  sequence: string;
  start: number;
  cell: number;
  qualities?: number[];
  reference?: string;
}) {
  return (
      <div className="alignment-sequence" style={{ width: sequence.length * cell }}>
        {[...sequence].map((base, index) => {
          const quality = qualities?.[index];
          const mismatch = reference !== undefined && base !== "·" && base !== reference[index];
          return (
            <span
              key={`${start + index}-${base}`}
              className={`alignment-base ${BASE_CLASS[base] ?? "base-n"}${mismatch ? " is-mismatch" : ""}${base === "·" ? " is-gap" : ""}`}
              style={{ width: cell, opacity: quality === undefined ? 1 : Math.max(0.38, Math.min(1, quality / 30)) }}
              title={quality === undefined
                ? `Position ${start + index + 1}: ${base === "·" ? "not covered" : base}`
                : `Position ${start + index + 1}: ${base}, Q${quality}`}
            >
              {base}
            </span>
          );
        })}
      </div>
  );
}

function SequenceRow({ label, ...strip }: Parameters<typeof BaseStrip>[0] & { label: string }) {
  return (
    <div className="alignment-row">
      <div className="alignment-label">{label}</div>
      <BaseStrip {...strip} />
    </div>
  );
}

function TracePlot({ track, cell }: { track: TraceTrack; cell: number }) {
  const width = Math.max(cell, track.sequence.length * cell);
  const height = 72;
  const peak = Math.max(1, ...Object.values(track.traces).flatMap((values) => values));
  const points = (values: number[]) => values.map((value, index) => {
    const x = track.sample_count > 1 ? (index / (track.sample_count - 1)) * width : 0;
    const y = height - 4 - (value / peak) * (height - 10);
    return `${x.toFixed(2)},${y.toFixed(2)}`;
  }).join(" ");

  return (
    <svg
      className="alignment-trace"
      width={width}
      height={height}
      viewBox={`0 0 ${width} ${height}`}
      role="img"
      aria-label={`${track.read}, ${track.reverse_complemented ? "reverse" : "forward"} chromatogram, ${track.sequence.length} quality-trimmed bases`}
    >
      {track.qualities.map((quality, index) => (
        <rect
          key={index}
          x={index * cell}
          y={0}
          width={cell}
          height={height}
          className={quality >= 20 ? "trace-quality-high" : "trace-quality-low"}
        />
      ))}
      {Object.entries(track.traces).map(([base, values]) => values.length ? (
        <polyline
          key={base}
          points={points(values)}
          fill="none"
          stroke={TRACE_COLOUR[base] ?? "#607185"}
          strokeWidth="1.25"
          vectorEffect="non-scaling-stroke"
        />
      ) : null)}
    </svg>
  );
}

function clipTrack(track: TraceTrack, firstBase: number, baseCount: number): TraceTrack {
  if (firstBase === 0 && baseCount === track.sequence.length) return track;

  const afterLastBase = firstBase + baseCount;
  const firstSample = firstBase === 0
    ? 0
    : Math.floor(((track.peaks[firstBase - 1] ?? 0) + (track.peaks[firstBase] ?? 0)) / 2);
  const lastSample = afterLastBase >= track.sequence.length
    ? track.sample_count
    : Math.ceil(((track.peaks[afterLastBase - 1] ?? track.sample_count) + (track.peaks[afterLastBase] ?? track.sample_count)) / 2);

  return {
    ...track,
    sequence: track.sequence.slice(firstBase, afterLastBase),
    qualities: track.qualities.slice(firstBase, afterLastBase),
    peaks: track.peaks.slice(firstBase, afterLastBase).map((peak) => peak - firstSample),
    sample_count: Math.max(0, lastSample - firstSample),
    traces: Object.fromEntries(
      Object.entries(track.traces).map(([base, values]) => [base, values.slice(firstSample, lastSample)]),
    ),
  };
}

/** Reference, consensus, oriented base calls and real Sanger peaks in one coordinate system. */
export default function ReferenceAlignment({ reference, report }: Props) {
  const tracks = report.trace_tracks ?? [];
  const [cell, setCell] = useState(18);
  const scrollRef = useRef<HTMLDivElement>(null);
  const start = report.region_start ?? 0;
  const end = report.region_end ?? reference.length;
  const span = Math.max(0, end - start);
  const referenceSlice = reference.slice(start, end);
  const calls = useMemo(
    () => consensus(reference, tracks, start, end),
    [reference, tracks, start, end],
  );
  const consensusSlice = calls.map((call) => call.base).join("");

  useEffect(() => {
    const firstEvidence = Math.min(...tracks.map((track) => track.reference_start));
    if (scrollRef.current && Number.isFinite(firstEvidence)) {
      scrollRef.current.scrollLeft = Math.max(0, (firstEvidence - start) * cell - cell * 3);
    }
  }, [cell, start, tracks]);

  if (!tracks.length || span === 0) return null;

  return (
    <section className="reference-alignment" aria-labelledby="reference-alignment-title">
      <div className="alignment-toolbar">
        <div>
          <h3 id="reference-alignment-title">Reference-aligned chromatograms</h3>
          <p>Only the quality-trimmed bases used by verification are mapped. Peaks are shown in reference orientation.</p>
        </div>
        <div className="alignment-zoom" role="group" aria-label="Alignment zoom">
          <button type="button" onClick={() => setCell(12)} aria-pressed={cell === 12}>Compact</button>
          <button type="button" onClick={() => setCell(18)} aria-pressed={cell === 18}>Standard</button>
          <button type="button" onClick={() => setCell(24)} aria-pressed={cell === 24}>Large</button>
        </div>
      </div>

      <div ref={scrollRef} className="alignment-scroll" tabIndex={0} aria-label="Scrollable sequencing alignment">
        <div className="alignment-canvas" style={{ width: 150 + span * cell }}>
          <div className="alignment-row alignment-ruler-row" aria-hidden="true">
            <div className="alignment-label">Position</div>
            <div className="alignment-ruler" style={{ width: span * cell }}>
              {Array.from({ length: span }, (_, index) => {
                const position = start + index + 1;
                return (position === start + 1 || position % 10 === 0) ? (
                  <span key={position} style={{ left: index * cell }}>{position}</span>
                ) : null;
              })}
            </div>
          </div>
          <SequenceRow label="Consensus" sequence={consensusSlice} start={start} cell={cell} reference={referenceSlice} />
          <SequenceRow label="Reference" sequence={referenceSlice} start={start} cell={cell} />

          {tracks.map((track) => {
            const offset = track.reference_start - start;
            const visibleStart = Math.max(0, offset);
            const clip = Math.max(0, -offset);
            const available = Math.max(0, Math.min(track.sequence.length - clip, span - visibleStart));
            if (!available) return null;
            const sequence = track.sequence.slice(clip, clip + available);
            const qualities = track.qualities.slice(clip, clip + available);
            const plotTrack = clipTrack(track, clip, available);
            return (
              <div className="alignment-read" key={track.read}>
                <div className="alignment-row alignment-read-head">
                  <div className="alignment-label" title={track.read}>
                    <strong>{track.read}</strong>
                    <span>{track.reverse_complemented ? "REV" : "FWD"} · Q{Math.round(track.qualities.reduce((sum, q) => sum + q, 0) / Math.max(1, track.qualities.length))}</span>
                  </div>
                  <div className="alignment-track" style={{ width: span * cell }}>
                    <div className="alignment-read-arrow" style={{ left: visibleStart * cell, width: available * cell }}>
                      <span>{track.reverse_complemented ? "reverse read" : "forward read"}</span>
                    </div>
                  </div>
                </div>
                <div className="alignment-row">
                  <div className="alignment-label">Base calls</div>
                  <div className="alignment-track" style={{ width: span * cell }}>
                    <div style={{ position: "absolute", left: visibleStart * cell }}>
                      <BaseStrip sequence={sequence} qualities={qualities} start={start + visibleStart} cell={cell} reference={reference.slice(start + visibleStart, start + visibleStart + available)} />
                    </div>
                  </div>
                </div>
                <div className="alignment-row alignment-plot-row">
                  <div className="alignment-label">Chromatogram</div>
                  <div className="alignment-track" style={{ width: span * cell }}>
                    <div style={{ position: "absolute", left: visibleStart * cell }}>
                      <TracePlot track={plotTrack} cell={cell} />
                    </div>
                  </div>
                </div>
              </div>
            );
          })}
        </div>
      </div>

      <div className="alignment-legend" aria-label="Alignment legend">
        <span><i className="legend-covered" />Q20 or higher</span>
        <span><i className="legend-low" />below Q20</span>
        <span><i className="legend-mismatch" />difference from reference</span>
        <span>A green · C blue · G black · T red</span>
      </div>
    </section>
  );
}
