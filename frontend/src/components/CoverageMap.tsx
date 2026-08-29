import type { VerifyReport } from "../api/client";

type Span = { start: number; end: number };

function pieces(span: Span, length: number): Span[] {
  if (length <= 0 || span.end <= span.start) return [];
  const width = span.end - span.start;
  if (width >= length) return [{ start: 0, end: length }];
  const start = ((span.start % length) + length) % length;
  const end = start + width;
  return end <= length
    ? [{ start, end }]
    : [{ start, end: length }, { start: 0, end: end - length }];
}

function style(span: Span, length: number) {
  return {
    left: `${(span.start / length) * 100}%`,
    width: `${Math.max(0.3, ((span.end - span.start) / length) * 100)}%`,
  };
}

/** Construct-coordinate overview: coverage, gaps and each placed read. */
export default function CoverageMap({ report }: { report: VerifyReport }) {
  const length = report.design_length;
  if (!length) return null;
  const region = {
    start: report.region_start ?? 0,
    end: report.region_end ?? length,
  };

  return (
    <div className="coverage-map" role="img" aria-label={`${report.coverage}% sequencing coverage; ${report.gaps.length} gap${report.gaps.length === 1 ? "" : "s"}`}>
      <div className="coverage-axis">
        {pieces(region, length).map((part, index) => (
          <span className="coverage-region" style={style(part, length)} key={`region-${index}`} />
        ))}
        {report.gaps.flatMap(([start, end]) => pieces({ start, end }, length)).map((part, index) => (
          <span className="coverage-gap" style={style(part, length)} key={`gap-${index}`} />
        ))}
      </div>
      <div className="coverage-reads" aria-hidden="true" style={{ height: `${Math.max(12, report.reads.length * 9)}px` }}>
        {report.reads.map((read, row) => pieces(read, length).map((part, index) => (
          <span
            className={`coverage-read ${read.reverse_complemented ? "reverse" : "forward"}`}
            style={{ ...style(part, length), top: `${row * 9}px` }}
            key={`${read.name}-${index}`}
            title={`${read.name}: ${read.start + 1}–${read.end}`}
          />
        )))}
      </div>
      <div className="coverage-labels"><span>1</span><span>{length.toLocaleString()} bp</span></div>
      <div className="coverage-legend" aria-hidden="true">
        <span><i className="forward" />forward read</span>
        <span><i className="reverse" />reverse read</span>
        <span><i className="gap" />uncovered</span>
      </div>
    </div>
  );
}
