import type { JunctionView } from "../api/client";

/**
 * One ligation seam, drawn as the two ends that made it.
 *
 * A banner saying "the overhangs match" asks to be believed. Showing the two
 * ends about to anneal, base against base, can be checked — and it is the
 * only rendering in which an overhang one base short, or the right sequence
 * on the wrong strand, looks obviously wrong instead of looking like a
 * passing test.
 */

type Props = {
  view: JunctionView;
  /** Show the pre-ligation ends as well as the joined result. */
  showEnds?: boolean;
};

function Strand({ text, span, className }: {
  text: string;
  span?: [number, number];
  className?: string;
}) {
  return (
    <span className={`dx-seq ${className ?? ""}`}>
      {[...text].map((base, index) => {
        const inOverhang = span && index >= span[0] && index < span[1];
        // A padding column must hold its width. A plain space collapses
        // wherever `white-space: pre` does not reach, and the two pieces
        // stop lining up — which is the one thing this drawing is for.
        return (
          <span
            key={index}
            className={
              base === " " ? "dx-gap" : inOverhang ? "dx-base jx-overhang" : "dx-base"
            }
          >
            {base === " " ? "\u00A0" : base}
          </span>
        );
      })}
    </span>
  );
}

export default function JunctionDuplex({ view, showEnds = true }: Props) {
  const span: [number, number] = [view.overhang_span[0], view.overhang_span[1]];

  return (
    <div className="junction-duplex">
      <div className="junction-head">
        <strong>{view.name}</strong>
        <span className="label">
          {view.enzyme} · {view.kind} {view.overhang || "blunt"}
        </span>
        <span className="grow" />
        <span className={view.compatible ? "pill pill-ok" : "pill pill-bad"}>
          {view.compatible ? "overhangs match" : "does not match"}
        </span>
      </div>

      {!view.compatible && <p className="note vector-note">{view.reason}</p>}

      {showEnds && (
        <>
          <div className="jx-caption">Before ligation — the two ends</div>
          <div className="duplex-scroll jx-block">
            <div className="duplex-row">
              <span className="dx-end">5'</span>
              <Strand text={view.left_top} />
              <span className="jx-space" />
              <Strand text={view.right_top} />
            </div>
            <div className="duplex-row">
              <span className="dx-end">3'</span>
              <Strand text={view.left_bottom} />
              <span className="jx-space" />
              <Strand text={view.right_bottom} />
            </div>
          </div>
          <p className="note jx-note">
            Each piece carries the overhang, on opposite strands — that is what
            lets them anneal. One that carried it on the same strand as its
            partner could not join.
          </p>
        </>
      )}

      <div className="jx-caption">After ligation</div>
      <div className="duplex-scroll jx-block">
        <div className="duplex-row">
          <span className="dx-end">5'</span>
          <Strand text={view.joined_top} span={span} />
        </div>
        <div className="duplex-row">
          <span className="dx-end" />
          <span className="dx-seq dx-ticks">
            {[...view.joined_pairs].map((mark, index) => (
              <span
                key={index}
                className={index >= span[0] && index < span[1] ? "dx-stagger" : ""}
              >
                {mark}
              </span>
            ))}
          </span>
        </div>
        <div className="duplex-row">
          <span className="dx-end">3'</span>
          <Strand text={view.joined_bottom} span={span} />
        </div>
      </div>
    </div>
  );
}
