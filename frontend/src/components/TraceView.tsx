import type { TraceWindow } from "../api/client";

/**
 * The peaks around one disputed base.
 *
 * The reason this exists: a read differs from the design at one position,
 * and whether that is a mutation or a bad call is not in the letters. It is
 * in the shape underneath them. One sharp peak clear of its neighbours is a
 * real base; a shoulder on the peak beside it is the basecaller guessing,
 * and it arrives spelled exactly the same way.
 *
 * So this draws what the instrument actually recorded — four channels, the
 * called bases at their peak positions, the disputed one marked — and lets
 * the reader make the judgement they would have made in another program.
 *
 * Only the window travels from the server. A full trace is ~15 000 samples
 * across four channels; sending them all to draw ten bases would be
 * megabytes per difference.
 */

const CHANNEL: Record<string, string> = {
  // The convention every sequencing viewer uses. Worth matching exactly:
  // people read these by colour before they read the letters.
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

  // One polyline per channel, scaled to the tallest peak in this window so a
  // weak read still fills the box rather than sitting flat against the axis.
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
        {/* The disputed base's own column, so the eye lands on it first. */}
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

      {/* The base calls, placed at their own peaks rather than evenly: an
          evenly-spaced row drifts out of register with the trace it labels. */}
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
