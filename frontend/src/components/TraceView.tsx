import type { TraceWindow } from "../api/client";


const CHANNEL: Record<string, string> = {


  A: "#3f7a52",
  C: "#1d63a8",
  G: "#16202c",
  T: "#a83232",
};

type Props = {
  window: TraceWindow;
  height?: number;
};

export default function TraceView({ window: view, height = 92 }: Props) {
  const [first, last] = view.samples;
  const width = last - first;
  if (width <= 0) return null;

  const peak = Math.max(
    1,
    ...Object.values(view.traces).flatMap((channel) => channel),
  );


  const line = (samples: number[]) =>
    samples
      .map((v, i) => `${(i / width) * 100},${height - (v / peak) * (height - 14)}`)
      .join(" ");

  const centre = view.bases.find((b) => b.index === view.centre);

  return (
    <div className="trace-view">
      <svg
        viewBox={`0 0 100 ${height}`}
        preserveAspectRatio="none"
        className="trace-plot"
        style={{ height }}
        role="img"
        aria-label={
          centre
            ? `Sequencing trace around base ${centre.base}, quality ${centre.quality}`
            : "Sequencing trace"
        }
      >

        {centre && (
          <rect
            x={((centre.at - 5) / width) * 100}
            y={0}
            width={(10 / width) * 100}
            height={height}
            className="trace-focus"
          />
        )}
        {Object.entries(view.traces).map(([base, samples]) =>
          samples.length ? (
            <polyline
              key={base}
              points={line(samples)}
              fill="none"
              stroke={CHANNEL[base] ?? "#78889b"}
              strokeWidth={0.6}
              vectorEffect="non-scaling-stroke"
            />
          ) : null,
        )}
      </svg>


      <div className="trace-bases" style={{ height: 20 }}>
        {view.bases.map((b) => (
          <span
            key={b.index}
            className={
              "trace-base" +
              (b.index === view.centre ? " is-centre" : "") +
              (b.quality < 20 ? " is-poor" : "")
            }
            style={{ left: `${(b.at / width) * 100}%`, color: CHANNEL[b.base] }}
            title={`Base ${b.index + 1} · ${b.base} · Q${b.quality}`}
          >
            {b.base}
          </span>
        ))}
      </div>

      {centre && (
        <p className="note trace-note">
          {centre.quality >= 20 ? (
            <>
              Called <b>{centre.base}</b> at Q{centre.quality} — a clean peak.
              This difference is in the construct, not in the read.
            </>
          ) : (
            <>
              Called <b>{centre.base}</b> at only Q{centre.quality}. At this
              confidence the basecaller is choosing between a peak and its
              neighbour&rsquo;s shoulder, so re-read before treating this as
              a real change.
            </>
          )}
        </p>
      )}
    </div>
  );
}
