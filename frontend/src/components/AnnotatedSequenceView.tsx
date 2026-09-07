import { useEffect, useMemo, useRef, useState, type CSSProperties } from "react";

import type { Annotation } from "../api/client";
import { featureLabel } from "./featureLabel";

const ROW_BASES = 60;
const WINDOW_PADDING = 160;
const MAX_WINDOW = 720;

const COMPLEMENT: Record<string, string> = {
  A: "T", T: "A", G: "C", C: "G", N: "N",
};

const CODON_TO_AA: Record<string, string> = {
  TTT: "Phe", TTC: "Phe", TTA: "Leu", TTG: "Leu",
  TCT: "Ser", TCC: "Ser", TCA: "Ser", TCG: "Ser",
  TAT: "Tyr", TAC: "Tyr", TAA: "Stop", TAG: "Stop",
  TGT: "Cys", TGC: "Cys", TGA: "Stop", TGG: "Trp",
  CTT: "Leu", CTC: "Leu", CTA: "Leu", CTG: "Leu",
  CCT: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
  CAT: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
  CGT: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
  ATT: "Ile", ATC: "Ile", ATA: "Ile", ATG: "Met",
  ACT: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
  AAT: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
  AGT: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
  GTT: "Val", GTC: "Val", GTA: "Val", GTG: "Val",
  GCT: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
  GAT: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
  GGT: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
};

type Window = { start: number; end: number };

type VisibleFeature = {
  annotation: Annotation;
  sourceAnnotation?: Annotation;
  start: number;
  end: number;
  lane: number;
};

type TranslationCell = {
  start: number;
  end: number;
  aminoAcid: string;
  codon: string;
  residue: number;
};

function clamp(value: number, low: number, high: number): number {
  return Math.min(high, Math.max(low, value));
}

/** Keep feature labels readable on both pale and saturated imported colours. */
export function readableTextColour(hex: string): "#ffffff" | "#0b1f3b" {
  const match = /^#([0-9a-f]{6})$/i.exec(hex.trim());
  if (!match) return "#0b1f3b";
  const value = Number.parseInt(match[1], 16);
  const red = (value >> 16) & 255;
  const green = (value >> 8) & 255;
  const blue = value & 255;
  const luminance = (0.2126 * red + 0.7152 * green + 0.0722 * blue) / 255;
  return luminance > 0.52 ? "#0b1f3b" : "#ffffff";
}

function reverseComplement(sequence: string): string {
  return sequence
    .toUpperCase()
    .split("")
    .reverse()
    .map((base) => COMPLEMENT[base] ?? "N")
    .join("");
}

/** Pick a useful locus rather than shrinking a 5 kb plasmid to illegibility. */
export function chooseAnnotationWindow(
  annotations: Annotation[],
  sequenceLength: number,
  focus: Annotation | null = null,
  preferredName = "",
  circular = false,
): Window {
  if (sequenceLength <= 0) return { start: 0, end: 0 };
  if (circular && sequenceLength <= MAX_WINDOW) return { start: 0, end: sequenceLength };

  const preferred = preferredName
    ? annotations.find((annotation) => annotation.name === preferredName)
    : undefined;
  const target = focus
    ?? preferred
    ?? annotations.find((annotation) => annotation.name.toLowerCase() === "insert")
    ?? annotations.find((annotation) => annotation.type.toLowerCase() === "cds")
    ?? annotations[0];

  if (!target) return { start: 0, end: Math.min(sequenceLength, 300) };

  let start = circular
    ? target.start - WINDOW_PADDING
    : clamp(target.start - WINDOW_PADDING, 0, sequenceLength);
  let end = circular
    ? target.end + WINDOW_PADDING
    : clamp(target.end + WINDOW_PADDING, 0, sequenceLength);

  // Bring nearby promoter/operator/RBS tracks into the same scientific view.
  const placements = circular
    ? annotations.flatMap((annotation) => [-sequenceLength, 0, sequenceLength].map((offset) => ({
      start: annotation.start + offset,
      end: annotation.end + offset,
      type: annotation.type,
    })))
    : annotations;
  const searchStart = start;
  const searchEnd = end;
  for (const annotation of placements) {
    // Compare with the initial padded locus. Expanding against an already
    // expanded range can chain from promoter to promoter around a plasmid
    // until the expression cassette is a tiny part of a needlessly huge view.
    if (
      annotation.end - annotation.start <= 300
      && annotation.end >= searchStart - 40
      && annotation.start <= searchEnd + 40
    ) {
      start = Math.min(
        start,
        circular ? annotation.start - 10 : clamp(annotation.start - 10, 0, sequenceLength),
      );
      end = Math.max(
        end,
        circular ? annotation.end + 10 : clamp(annotation.end + 10, 0, sequenceLength),
      );
    }
  }

  if (end - start > MAX_WINDOW) {
    const targetMidpoint = (target.start + target.end) / 2;
    if (circular) {
      start = Math.round(targetMidpoint - MAX_WINDOW / 2);
      end = start + MAX_WINDOW;
    } else {
      start = clamp(Math.round(targetMidpoint - MAX_WINDOW / 2), 0, sequenceLength);
      end = Math.min(sequenceLength, start + MAX_WINDOW);
      start = Math.max(0, end - MAX_WINDOW);
    }
  }

  start = Math.floor(start / 10) * 10;
  end = Math.ceil(end / 10) * 10;
  if (!circular) {
    start = Math.max(0, start);
    end = Math.min(sequenceLength, end);
  }
  return { start, end };
}

/** Greedy interval colouring keeps overlapping biological features legible. */
export function visibleFeaturesForRow(
  annotations: Annotation[], rowStart: number, rowEnd: number,
): VisibleFeature[] {
  const visible = annotations
    .filter((annotation) => annotation.end > rowStart && annotation.start < rowEnd)
    .map((annotation) => ({
      annotation,
      start: Math.max(rowStart, annotation.start),
      end: Math.min(rowEnd, annotation.end),
      lane: 0,
    }))
    .sort((left, right) => left.start - right.start || right.end - left.end);

  const laneEnds: number[] = [];
  for (const feature of visible) {
    let lane = laneEnds.findIndex((end) => end <= feature.start);
    if (lane < 0) lane = laneEnds.length;
    feature.lane = lane;
    laneEnds[lane] = feature.end;
  }
  return visible;
}

type PlacedAnnotation = Annotation & { sourceAnnotation: Annotation };

function placedAnnotationsForWindow(
  annotations: Annotation[], sequenceLength: number, window: Window, circular: boolean,
): PlacedAnnotation[] {
  const offsets = circular ? [-sequenceLength, 0, sequenceLength] : [0];
  return annotations.flatMap((annotation) => offsets
    .map((offset) => ({
      ...annotation,
      start: annotation.start + offset,
      end: annotation.end + offset,
      translation_start: annotation.translation_start === undefined
        ? undefined
        : annotation.translation_start + offset,
      translation_end: annotation.translation_end === undefined
        ? undefined
        : annotation.translation_end + offset,
      sourceAnnotation: annotation,
    }))
    .filter((annotation) => annotation.end > window.start && annotation.start < window.end));
}

function virtualSequence(sequence: string, start: number, end: number): string {
  if (!sequence.length || end <= start) return "";
  return Array.from({ length: end - start }, (_, index) => {
    const position = ((start + index) % sequence.length + sequence.length) % sequence.length;
    return sequence[position];
  }).join("");
}

function translationCells(
  sequence: string,
  feature: Annotation,
  rowStart: number,
  rowEnd: number,
): TranslationCell[] {
  const featureStart = feature.translation_start ?? feature.start;
  const featureEnd = feature.translation_end ?? feature.end;
  const cells: TranslationCell[] = [];

  if (feature.direction === -1) {
    let residue = Math.max(1, Math.floor((featureEnd - rowEnd) / 3) + 1);
    for (let codonEnd = featureEnd - (residue - 1) * 3; codonEnd - 3 >= featureStart && codonEnd > rowStart; codonEnd -= 3) {
      const codonStart = codonEnd - 3;
      if (codonEnd > rowStart && codonStart < rowEnd) {
        const codon = reverseComplement(virtualSequence(sequence, codonStart, codonEnd));
        cells.push({
          start: Math.max(rowStart, codonStart),
          end: Math.min(rowEnd, codonEnd),
          aminoAcid: CODON_TO_AA[codon] ?? "Xaa",
          codon,
          residue,
        });
      }
      residue += 1;
    }
    return cells;
  }

  let residue = Math.max(1, Math.floor((rowStart - featureStart) / 3) + 1);
  for (let codonStart = featureStart + (residue - 1) * 3; codonStart + 3 <= featureEnd && codonStart < rowEnd; codonStart += 3) {
    const codonEnd = codonStart + 3;
    if (codonEnd > rowStart && codonStart < rowEnd) {
      const codon = virtualSequence(sequence, codonStart, codonEnd).toUpperCase();
      cells.push({
        start: Math.max(rowStart, codonStart),
        end: Math.min(rowEnd, codonEnd),
        aminoAcid: CODON_TO_AA[codon] ?? "Xaa",
        codon,
        residue,
      });
    }
    residue += 1;
  }
  return cells;
}

type Props = {
  sequence: string;
  annotations: Annotation[];
  selected: Annotation | null;
  preferredName?: string;
  circular?: boolean;
  onSelect: (annotation: Annotation) => void;
  onAnnotateRange?: (range: { start: number; end: number }) => void;
};

export default function AnnotatedSequenceView({
  sequence, annotations, selected, preferredName = "", circular = false, onSelect, onAnnotateRange,
}: Props) {
  const [rangeAnchor, setRangeAnchor] = useState<number | null>(null);
  const [rangeEnd, setRangeEnd] = useState<number | null>(null);
  const dragging = useRef(false);
  useEffect(() => {
    const finish = () => { dragging.current = false; };
    globalThis.addEventListener("pointerup", finish);
    globalThis.addEventListener("pointercancel", finish);
    return () => { globalThis.removeEventListener("pointerup", finish); globalThis.removeEventListener("pointercancel", finish); };
  }, []);
  useEffect(() => { setRangeAnchor(null); setRangeEnd(null); }, [sequence]);
  const range = rangeAnchor === null || rangeEnd === null ? null : {
    start: Math.min(rangeAnchor, rangeEnd),
    end: Math.min(Math.max(rangeAnchor, rangeEnd) + 1, Math.min(rangeAnchor, rangeEnd) + sequence.length),
  };
  function selectBase(coordinate: number, extend: boolean) {
    if (!extend || rangeAnchor === null) setRangeAnchor(coordinate);
    setRangeEnd(coordinate);
  }
  const [scope, setScope] = useState<"locus" | "whole">("locus");
  const [manualWindow, setManualWindow] = useState<Window | null>(null);
  const [position, setPosition] = useState("");
  const [navigationError, setNavigationError] = useState("");
  const [rowBases, setRowBases] = useState(ROW_BASES);
  const [scrollTop, setScrollTop] = useState(0);
  const [viewportHeight, setViewportHeight] = useState(660);
  const scrollRef = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const node = scrollRef.current;
    if (!node || typeof ResizeObserver === "undefined") return;
    const observer = new ResizeObserver(() => setViewportHeight(node.clientHeight || 660));
    observer.observe(node);
    return () => observer.disconnect();
  }, []);
  useEffect(() => { setManualWindow(null); }, [selected, sequence]);
  const locus = useMemo(
    () => chooseAnnotationWindow(
      annotations, sequence.length, selected, preferredName, circular,
    ),
    [annotations, circular, preferredName, selected, sequence.length],
  );

  const window = useMemo(() => scope === "whole"
    ? { start: 0, end: sequence.length } : manualWindow ?? locus,
  [scope, sequence.length, manualWindow, locus]);
  useEffect(() => {
    setScrollTop(0);
    if (scrollRef.current) scrollRef.current.scrollTop = 0;
  }, [window.start, window.end, rowBases]);

  const placedAnnotations = useMemo(
    () => placedAnnotationsForWindow(annotations, sequence.length, window, circular),
    [annotations, circular, sequence.length, window],
  );

  const rows = useMemo(() => {
    const starts: number[] = [];
    for (let start = window.start; start < window.end; start += rowBases) starts.push(start);
    return starts;
  }, [window, rowBases]);

  const codingFeature = useMemo(() => {
    const candidates = placedAnnotations.filter((annotation) =>
      annotation.type.toLowerCase() === "cds"
      && annotation.end > window.start
      && annotation.start < window.end,
    );
    if (selected?.type.toLowerCase() === "cds") {
      return candidates.find((annotation) => annotation.sourceAnnotation === selected)
        ?? candidates[0]
        ?? null;
    }
    return candidates.find((annotation) => annotation.name === preferredName) ?? candidates[0] ?? null;
  }, [placedAnnotations, preferredName, selected, window]);

  const blocks = useMemo(() => {
    let offset = 0;
    return rows.map((rowStart) => {
      const rowEnd = Math.min(window.end, rowStart + rowBases);
      const visible = visibleFeaturesForRow(placedAnnotations, rowStart, rowEnd);
      const laneCount = Math.max(1, ...visible.map((feature) => feature.lane + 1));
      const translations = codingFeature ? translationCells(sequence, codingFeature, rowStart, rowEnd) : [];
      const height = 76 + laneCount * 28 + (translations.length ? 29 : 0);
      const block = { rowStart, rowEnd, visible, laneCount, translations, offset, height };
      offset += height;
      return block;
    });
  }, [rows, window.end, rowBases, placedAnnotations, codingFeature, sequence]);
  const totalHeight = blocks.length ? blocks[blocks.length - 1].offset + blocks[blocks.length - 1].height : 0;
  const visibleBlocks = blocks.filter((block) => block.offset + block.height >= scrollTop - 350
    && block.offset <= scrollTop + viewportHeight + 350);

  function moveRegion(start: number) {
    const size = Math.min(600, sequence.length);
    const nextStart = clamp(start, 0, Math.max(0, sequence.length - size));
    setScope("locus");
    setManualWindow({ start: nextStart, end: nextStart + size });
  }
  function goToPosition(event: React.FormEvent) {
    event.preventDefault();
    const coordinate = Number(position);
    if (!Number.isInteger(coordinate) || coordinate < 1 || coordinate > sequence.length) {
      setNavigationError(`Enter a base from 1 to ${sequence.length.toLocaleString()}.`);
      return;
    }
    setNavigationError("");
    if (scope === "whole") {
      const block = blocks.find((item) => item.rowStart <= coordinate - 1 && item.rowEnd > coordinate - 1);
      if (block && scrollRef.current) { scrollRef.current.scrollTop = block.offset; setScrollTop(block.offset); }
    } else moveRegion(coordinate - 1 - rowBases);
  }

  const displayCoordinate = (position: number) => {
    if (!circular || !sequence.length) return position + 1;
    return ((position % sequence.length + sequence.length) % sequence.length) + 1;
  };

  return (
    <section className="annotated-sequence" aria-label="Coordinate-level annotated sequence">
      <header className="annotated-sequence-head">
        <div>
          <h2>Annotated sequence</h2>
          <p>
            Showing bases {displayCoordinate(window.start).toLocaleString()}–
            {displayCoordinate(window.end - 1).toLocaleString()}
            {circular && (window.end > sequence.length || window.start < 0) ? " across the circular origin" : ""}
            {codingFeature ? ` · translation from ${codingFeature.name}` : ""}
          </p>
        </div>
        <span className="label">1-based coordinates</span>
      </header>

      <div className="sequence-navigation" aria-label="Sequence navigation">
        <div className="seg-toggle" role="group" aria-label="Sequence range">
          <button type="button" className={scope === "locus" ? "on" : ""} aria-pressed={scope === "locus"}
            onClick={() => { setScope("locus"); setManualWindow(null); }}>Feature context</button>
          <button type="button" className={scope === "whole" ? "on" : ""} aria-pressed={scope === "whole"}
            onClick={() => setScope("whole")}>{circular ? "Whole plasmid" : "Whole sequence"}</button>
        </div>
        <label>Bases per row <select value={rowBases} onChange={(event) => setRowBases(Number(event.target.value))}>
          {[30, 60, 90].map((value) => <option key={value} value={value}>{value}</option>)}
        </select></label>
        <form onSubmit={goToPosition}>
          <label>Go to base <input type="text" inputMode="numeric" value={position} onChange={(event) => setPosition(event.target.value)} placeholder={`1–${sequence.length}`} /></label>
          <button type="submit" className="btn btn-outline">Go</button>
        </form>
        {scope === "locus" && <div className="sequence-page-buttons">
          <button type="button" className="btn btn-outline" disabled={window.start <= 0} onClick={() => moveRegion(window.start - 600)}>Previous region</button>
          <button type="button" className="btn btn-outline" disabled={window.end >= sequence.length} onClick={() => moveRegion(window.end)}>Next region</button>
        </div>}
        {navigationError && <p role="alert">{navigationError}</p>}
      </div>
      {onAnnotateRange && <div className="sequence-selection-toolbar" role="group" aria-label="Nucleotide selection">
        <span>{range ? `${displayCoordinate(range.start)}–${displayCoordinate(range.end - 1)} · ${range.end - range.start} nt selected` : "Drag across bases, or click the first base then Shift-click the last."}</span>
        <button type="button" className="btn btn-primary" disabled={!range} onClick={() => {
          if (!range) return;
          const start = ((range.start % sequence.length) + sequence.length) % sequence.length;
          onAnnotateRange({ start, end: start + range.end - range.start });
        }}>Annotate selection</button>
        {range && <button type="button" className="btn btn-ghost" onClick={() => { setRangeAnchor(null); setRangeEnd(null); }}>Clear range</button>}
      </div>}
      <div className="annotated-sequence-scroll" tabIndex={0} ref={scrollRef}
        aria-label="Scrollable annotated DNA" onScroll={(event) => setScrollTop(event.currentTarget.scrollTop)}>
        <div style={{ height: totalHeight, position: "relative", minWidth: rowBases * 13,
          "--sequence-row-bases": rowBases } as CSSProperties}>
        {visibleBlocks.map(({ rowStart, rowEnd, visible, laneCount, translations, offset, height }) => {
          const bases = virtualSequence(sequence, rowStart, rowEnd).toUpperCase().split("");

          return (
            <div className="annotation-block" key={rowStart} style={{ position: "absolute", top: offset, height, width: "100%" }}>
              <div className="annotation-ruler" aria-hidden="true">
                {Array.from({ length: rowEnd - rowStart }, (_, index) => {
                  const coordinate = rowStart + index;
                  const shownCoordinate = displayCoordinate(coordinate);
                  const show = index === 0 || shownCoordinate % 10 === 0 || shownCoordinate === 1;
                  return (
                    <span className={show ? "tick major" : "tick"} key={coordinate}>
                      {show ? shownCoordinate.toLocaleString() : ""}
                    </span>
                  );
                })}
              </div>

              <div
                className="annotation-tracks"
                role="group"
                style={{ height: `${laneCount * 28}px` }}
                aria-label={`Features across bases ${displayCoordinate(rowStart)} to ${displayCoordinate(rowEnd - 1)}`}
              >
                {visible.map((feature, index) => {
                  const width = feature.end - feature.start;
                  const clippedLeft = feature.annotation.start < rowStart;
                  const clippedRight = feature.annotation.end > rowEnd;
                  const direction = feature.annotation.direction === -1 ? "←" : feature.annotation.direction === 1 ? "→" : "";
                  return (
                    <button
                      type="button"
                      className={`annotation-feature ${feature.annotation.inferred ? "inferred" : ""} ${selected === (feature.annotation as PlacedAnnotation).sourceAnnotation ? "selected" : ""}`}
                      key={`${feature.annotation.name}-${feature.annotation.start}-${index}`}
                      style={{
                        left: `${100 * (feature.start - rowStart) / rowBases}%`,
                        width: `${100 * width / rowBases}%`,
                        top: `${feature.lane * 28}px`,
                        backgroundColor: feature.annotation.color,
                        color: readableTextColour(feature.annotation.color),
                      }}
                      onClick={() => onSelect((feature.annotation as PlacedAnnotation).sourceAnnotation)}
                      title={`${feature.annotation.inferred ? "Candidate · function unconfirmed · " : ""}${feature.annotation.name}: ${displayCoordinate(feature.annotation.start)}–${displayCoordinate(feature.annotation.end - 1)} (${feature.annotation.direction === -1 ? "reverse" : feature.annotation.direction === 1 ? "forward" : "unstranded"})`}
                      aria-label={`${feature.annotation.name}, ${feature.annotation.type}, bases ${displayCoordinate(feature.annotation.start)} to ${displayCoordinate(feature.annotation.end - 1)}, ${feature.annotation.direction === -1 ? "reverse" : feature.annotation.direction === 1 ? "forward" : "unstranded"} strand${feature.annotation.inferred ? ", candidate, function unconfirmed" : ""}`}
                    >
                      <span aria-hidden="true">{featureLabel(feature.annotation)} {!clippedLeft && direction}{clippedRight ? "…" : ""}</span>
                    </button>
                  );
                })}
              </div>

              <div className="annotation-bases" aria-label={`Sequence bases ${displayCoordinate(rowStart)} to ${displayCoordinate(rowEnd - 1)}`}>
                {bases.map((base, index) => (
                  <span
                    className={`annotation-base base-${base}${range && rowStart + index >= range.start && rowStart + index < range.end ? " range-selected" : ""}`}
                    key={rowStart + index}
                    title={`${displayCoordinate(rowStart + index).toLocaleString()}: ${base}`}
                    data-base-position={rowStart + index}
                    role={onAnnotateRange ? "button" : undefined}
                    tabIndex={onAnnotateRange ? (index === 0 ? 0 : -1) : undefined}
                    aria-label={onAnnotateRange ? `Base ${displayCoordinate(rowStart + index)}: ${base}` : undefined}
                    aria-pressed={onAnnotateRange ? Boolean(range && rowStart + index >= range.start && rowStart + index < range.end) : undefined}
                    onPointerDown={onAnnotateRange ? (event) => {
                      if (event.button !== 0) return;
                      event.preventDefault(); dragging.current = true;
                      selectBase(rowStart + index, event.shiftKey);
                    } : undefined}
                    onPointerEnter={onAnnotateRange ? () => { if (dragging.current) setRangeEnd(rowStart + index); } : undefined}
                    onKeyDown={onAnnotateRange ? (event) => {
                      if (event.key === "Enter" || event.key === " ") { event.preventDefault(); selectBase(rowStart + index, event.shiftKey); }
                      const step = { ArrowLeft: -1, ArrowRight: 1, ArrowUp: -rowBases, ArrowDown: rowBases }[event.key];
                      if (step !== undefined) {
                        event.preventDefault();
                        const next = rowStart + index + step;
                        const base = scrollRef.current?.querySelector<HTMLElement>(`[data-base-position="${next}"]`);
                        if (base) { base.focus(); if (event.shiftKey) selectBase(next, true); }
                      }
                    } : undefined}
                  >
                    {base}
                  </span>
                ))}
              </div>

              {translations.length > 0 && (
                <div className="annotation-translation" aria-label={`Translation of ${codingFeature?.name ?? "CDS"}`}>
                  {translations.map((cell) => (
                    <span
                      className={cell.aminoAcid === "Stop" ? "stop" : ""}
                      key={`${cell.start}-${cell.residue}`}
                      style={{
                        gridColumn: `${cell.start - rowStart + 1} / span ${Math.max(1, cell.end - cell.start)}`,
                      }}
                      title={`Residue ${cell.residue}: ${cell.aminoAcid} (${cell.codon})`}
                    >
                      {cell.aminoAcid === "Stop" ? "*" : cell.aminoAcid}
                    </span>
                  ))}
                </div>
              )}
            </div>
          );
        })}
        </div>
      </div>

      <footer className="annotation-legend" aria-label="Nucleotide colour legend">
        {(["A", "C", "G", "T"] as const).map((base) => (
          <span key={base}><i className={`base-${base}`} />{base}</span>
        ))}
        <span className="annotation-legend-note">Drag the lower edge to resize. ? and dashed borders mark candidates with unconfirmed function.</span>
      </footer>
    </section>
  );
}
