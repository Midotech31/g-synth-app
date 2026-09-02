import type { PcrPrimer } from "../api/client";

const COMPLEMENT: Record<string, string> = { A: "T", T: "A", G: "C", C: "G" };

export function annealingRows(primer: PcrPrimer) {
  const unpaired = " ".repeat(primer.tail.length);
  const template = primer.anneals
    .split("")
    .map((base) => COMPLEMENT[base] ?? "N")
    .join("");
  return {
    primer: primer.sequence,
    pairs: unpaired + "|".repeat(primer.anneals.length),
    template: unpaired + template,
    unpairedTailLength: primer.tail.length,
  };
}

function AnnealingDiagram({ primer, label }: { primer: PcrPrimer; label: string }) {
  const rows = annealingRows(primer);
  return (
    <section className="annealing-diagram" aria-label={`${label} primer annealing`}>
      <div className="annealing-title">
        <strong>{label} primer</strong>
        <span>{primer.start + 1}–{primer.end} on template</span>
      </div>
      <div className="annealing-scroll" tabIndex={0} aria-label={`${label} primer aligned to template`}>
        <div className="annealing-row">
          <span className="strand-label">Primer</span><b>5′</b>
          <code>
            {primer.tail && <span className="annealing-unpaired">{primer.tail}</span>}
            <span className="annealing-paired">{primer.anneals}</span>
          </code>
          <b>3′</b>
        </div>
        <div className="annealing-row annealing-pairs" aria-hidden="true">
          <span className="strand-label" /><b />
          <code>{rows.pairs}</code><b />
        </div>
        <div className="annealing-row">
          <span className="strand-label">Template</span><b>3′</b>
          <code>{rows.template}</code>
          <b>5′</b>
        </div>
      </div>
      {primer.tail ? (
        <p className="annealing-callout">
          <strong>{primer.tail.length} nt at the primer’s 5′ end are intentionally unpaired in cycle 1.</strong>
          {" "}Only the {primer.anneal_length}-nt 3′ region hybridizes and can be extended.
        </p>
      ) : (
        <p className="annealing-callout">
          All {primer.anneal_length} primer bases hybridize to the template in cycle 1.
        </p>
      )}
    </section>
  );
}

export default function PrimerAnnealingView({ forward, reverse }: {
  forward: PcrPrimer;
  reverse: PcrPrimer;
}) {
  const hasTail = Boolean(forward.tail || reverse.tail);
  return (
    <div className="primer-annealing-view">
      <div className="annealing-summary">
        <div className="annealing-step current">
          <span>Cycle 1</span>
          <strong>{hasTail ? "3′ annealing regions bind" : "Primers bind"}</strong>
          <small>{hasTail ? "5′ cloning tails remain single-stranded" : "The complete primers are complementary"}</small>
        </div>
        <div className="annealing-arrow" aria-hidden="true">→</div>
        <div className="annealing-step">
          <span>Extension</span>
          <strong>Polymerase copies through</strong>
          <small>The 5′ tails become part of the new product</small>
        </div>
        <div className="annealing-arrow" aria-hidden="true">→</div>
        <div className="annealing-step">
          <span>Cycles 2+</span>
          <strong>Full-length products amplify</strong>
          <small>Tail-derived restriction sites now have complements</small>
        </div>
      </div>
      <AnnealingDiagram primer={forward} label="Forward" />
      <AnnealingDiagram primer={reverse} label="Reverse" />
      <div className="duplex keys annealing-key">
        <span className="key"><i style={{ background: "var(--amber)" }} /> 5′ tail — not hybridized in cycle 1</span>
        <span className="key"><i style={{ background: "var(--accent)" }} /> 3′ annealing region — base-paired</span>
        <span className="key"><code>|</code> Watson–Crick base pair</span>
      </div>
    </div>
  );
}
