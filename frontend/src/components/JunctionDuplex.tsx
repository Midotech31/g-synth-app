import type { JunctionView } from "../api/client";


type Props = {
  view: JunctionView;
  showEnds?: boolean;
  ligated?: boolean;
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

export default function JunctionDuplex({ view, showEnds = true, ligated = false }: Props) {
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

      <div className="jx-caption">{ligated ? "Ligated junction" : "Expected junction after ligation"}</div>
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
