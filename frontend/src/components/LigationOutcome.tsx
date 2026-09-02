import type { JunctionView } from "../api/client";
import Icon from "./Icon";

type Props = {
  vectorName: string;
  backboneLength: number;
  insertName: string;
  insertLength: number;
  productName: string;
  productLength: number;
  junctions: JunctionView[];
};

function JoinedSequence({ sequence, span }: { sequence: string; span: [number, number] }) {
  const [start, end] = span;
  return (
    <code>
      {sequence.slice(0, start)}
      <mark>{sequence.slice(start, end)}</mark>
      {sequence.slice(end)}
    </code>
  );
}

export default function LigationOutcome({
  vectorName,
  backboneLength,
  insertName,
  insertLength,
  productName,
  productLength,
  junctions,
}: Props) {
  return (
    <section
      className="ligation-outcome"
      role="status"
      aria-live="polite"
      aria-atomic="true"
      aria-label="Ligation result"
    >
      <div className="ligation-outcome-verdict">
        <span className="ligation-outcome-icon" aria-hidden="true">
          <Icon name="check" size={24} />
        </span>
        <div>
          <span className="eyebrow">Ligation confirmed</span>
          <strong>{productName} · {productLength.toLocaleString()} bp</strong>
          <small>Both insert–vector junctions are closed in the recombinant product.</small>
        </div>
        <span className="pill pill-ok">2 junctions verified</span>
      </div>

      <div className="ligation-product-flow" aria-label={`${vectorName} plus ${insertName} forms ${productName}`}>
        <div>
          <span>Linearized vector</span>
          <strong>{vectorName}</strong>
          <small>{backboneLength.toLocaleString()} bp</small>
        </div>
        <b aria-hidden="true">+</b>
        <div>
          <span>Duplex insert</span>
          <strong>{insertName}</strong>
          <small>{insertLength.toLocaleString()} bp</small>
        </div>
        <Icon name="arrowRight" size={22} />
        <div className="ligation-product">
          <span>Recombinant product</span>
          <strong>{productName}</strong>
          <small>{productLength.toLocaleString()} bp · circular</small>
        </div>
      </div>

      <div className="ligation-junction-closeups">
        {junctions.map((junction) => (
          <div key={junction.name} className="ligation-junction-closeup">
            <div>
              <span className="ligation-junction-check"><Icon name="check" size={14} /> Closed</span>
              <strong>{junction.name}</strong>
              <small>{junction.enzyme} · {junction.kind} {junction.overhang || "blunt"}</small>
            </div>
            <div className="ligation-joined-duplex" aria-label={`${junction.name} ligated duplex`}>
              <span>5′</span><JoinedSequence sequence={junction.joined_top} span={junction.overhang_span} /><span>3′</span>
              <span>3′</span><JoinedSequence sequence={junction.joined_bottom} span={junction.overhang_span} /><span>5′</span>
            </div>
          </div>
        ))}
      </div>
    </section>
  );
}
