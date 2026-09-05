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

export function PrimerTail({ primer }: { primer: PcrPrimer }) {
  if (!primer.tail) return null;
  const siteAt = primer.restriction_site
    ? primer.tail.indexOf(primer.restriction_site)
    : -1;
  if (siteAt < 0) {
    return <span className="annealing-clamp">{primer.tail}</span>;
  }
  const siteEnd = siteAt + primer.restriction_site.length;
  return (
    <>
      {siteAt > 0 && <span className="annealing-clamp">{primer.tail.slice(0, siteAt)}</span>}
      <span className="annealing-site">{primer.tail.slice(siteAt, siteEnd)}</span>
      {siteEnd < primer.tail.length && (
        <span className="annealing-spacer">{primer.tail.slice(siteEnd)}</span>
      )}
    </>
  );
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
            {primer.tail && <span className="annealing-unpaired"><PrimerTail primer={primer} /></span>}
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

export default function PrimerAnnealingView({ forward, reverse, templateLength }: {
  forward: PcrPrimer;
  reverse: PcrPrimer;
  templateLength?: number;
}) {
  const hasTail = Boolean(forward.tail || reverse.tail);
  const span = Math.max(templateLength ?? reverse.end, reverse.end, forward.end, 1);
  const forwardLeft = Math.max(0, Math.min(100, (forward.start / span) * 100));
  const forwardWidth = Math.max(2, ((forward.end - forward.start) / span) * 100);
  const reverseLeft = Math.max(0, Math.min(100, (reverse.start / span) * 100));
  const reverseWidth = Math.max(2, ((reverse.end - reverse.start) / span) * 100);
  return (
    <div className="primer-annealing-view">
      <div className="primer-target-map" aria-label={`Primer positions on ${span}-base target`}>
        <div className="primer-target-labels">
          <strong>Target sequence</strong>
          <span>{span} bp · 5′→3′</span>
        </div>
        <div className="primer-target-track">
          <span
            className="primer-target-hit forward"
            style={{ left: `${forwardLeft}%`, width: `${forwardWidth}%` }}
            title={`Forward primer: ${forward.start + 1}–${forward.end}`}
          >→</span>
          <span
            className="primer-target-hit reverse"
            style={{ left: `${reverseLeft}%`, width: `${reverseWidth}%` }}
            title={`Reverse primer: ${reverse.start + 1}–${reverse.end}`}
          >←</span>
        </div>
        <div className="primer-target-coordinates">
          <span>Forward {forward.start + 1}–{forward.end} →</span>
          <span>← Reverse {reverse.start + 1}–{reverse.end}</span>
        </div>
      </div>
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
        {hasTail && <span className="key"><i className="key-clamp" /> terminal clamp — unpaired</span>}
        {hasTail && <span className="key"><i className="key-site" /> restriction site — unpaired</span>}
        <span className="key"><i style={{ background: "var(--accent)" }} /> 3′ annealing region — base-paired</span>
        <span className="key"><code>|</code> Watson–Crick base pair</span>
      </div>
    </div>
  );
}
