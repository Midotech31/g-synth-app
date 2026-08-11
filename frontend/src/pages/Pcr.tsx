import { useEffect, useMemo, useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";

import { ApiError, api, type Catalogue, type PcrResult } from "../api/client";
import Icon from "../components/Icon";

/**
 * PCR, and what cutting the product leaves.
 *
 * The page is deliberately one column of stages rather than a form and a
 * result: the whole point is that each step feeds the next, and a reader who
 * cannot see the tail becoming a site, and the site becoming an overhang, has
 * to take the insert on trust.
 */

/** A primer drawn so the tail and the annealing part are visibly different
 *  things — they are ordered as one oligo but behave as two, and the Tm that
 *  matters belongs to only one of them. */
function PrimerRow({ primer, label }: { primer: PcrResult["forward"]; label: string }) {
  return (
    <div className="primer-card">
      <div className="primer-head">
        <strong>{primer.name}</strong>
        <span className="label">{label}</span>
        <span className="grow" />
        <span className="note nums">{primer.length} nt</span>
      </div>

      <div className="primer-seq mono">
        {primer.tail && (
          <span className="seq-tail" title={`5' tail — ${primer.enzyme ?? "addition"}`}>
            {primer.tail}
          </span>
        )}
        <span className="seq-anneal" title="Anneals to the template">
          {primer.anneals}
        </span>
      </div>

      <div className="primer-facts">
        <span>
          Tm <b className="nums">{primer.tm.toFixed(1)} °C</b>
          <span className="note"> annealing part</span>
        </span>
        {primer.tail && (
          <span>
            Tm <b className="nums">{primer.tm_full.toFixed(1)} °C</b>
            <span className="note"> whole oligo</span>
          </span>
        )}
        <span>GC <b className="nums">{primer.gc.toFixed(0)}%</b></span>
        {primer.enzyme && <span className="pill pill-enzyme">{primer.enzyme}</span>}
      </div>

      {primer.warnings.length > 0 && (
        <ul className="primer-warnings">
          {primer.warnings.map((w) => <li key={w}>{w}</li>)}
        </ul>
      )}
    </div>
  );
}

/** The product with its parts coloured: tail, gene, tail. */
function ProductView({ result }: { result: PcrResult }) {
  const { product, forward, reverse } = result;
  const leftTail = forward.tail.length;
  const rightTail = reverse.tail.length;
  const body = product.slice(leftTail, product.length - rightTail);

  return (
    <div className="seq-block product-view mono">
      {leftTail > 0 && <span className="seq-tail">{product.slice(0, leftTail)}</span>}
      <span className="seq-body">{body}</span>
      {rightTail > 0 && <span className="seq-tail">{product.slice(product.length - rightTail)}</span>}
    </div>
  );
}

export default function Pcr() {
  const navigate = useNavigate();
  const location = useLocation() as { state?: { sequence?: string } | null };

  const [template, setTemplate] = useState(location.state?.sequence ?? "");
  const [mode, setMode] = useState<"conventional" | "cloning">("cloning");
  const [leftEnzyme, setLeftEnzyme] = useState("NdeI");
  const [rightEnzyme, setRightEnzyme] = useState("XhoI");
  const [clamp, setClamp] = useState(6);
  const [keepFrame, setKeepFrame] = useState(true);
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);

  const [result, setResult] = useState<PcrResult | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");

  useEffect(() => {
    api.catalogue().then(setCatalogue).catch(() => setCatalogue(null));
  }, []);

  const cleanLength = useMemo(
    () => template.replace(/[^A-Za-z]/g, "").length, [template],
  );

  async function run() {
    setBusy(true);
    setError("");
    setResult(null);
    try {
      setResult(await api.pcr({
        template,
        left_enzyme: mode === "cloning" ? leftEnzyme : null,
        right_enzyme: mode === "cloning" ? rightEnzyme : null,
        clamp,
        keep_frame: mode === "cloning" && keepFrame,
      }));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "The design could not be run.");
    } finally {
      setBusy(false);
    }
  }

  /** Hand the cut insert to Clone. Both strands travel: the stagger between
   *  them is the overhang, and Clone measures the ends off them rather than
   *  trusting what it is told. */
  function sendToClone() {
    if (!result?.digest) return;
    navigate("/clone", {
      state: {
        preDigested: {
          top: result.digest.top,
          bottom: result.digest.bottom,
          leftEnzyme: result.left_enzyme,
          rightEnzyme: result.right_enzyme,
        },
      },
    });
  }

  const enzymes = catalogue?.enzymes ?? [];

  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>PCR</h1>
          <p className="sub">
            Amplify a region, or add restriction sites to it and cut the product
            into an insert.
          </p>
        </div>
        <button className="btn btn-primary" onClick={() => void run()} disabled={busy || !cleanLength}>
          {busy && <span className="spinner" />}
          {busy ? "Designing…" : "Design PCR"}
        </button>
      </div>

      <div className="content design-layout">
        <div className="card">
          <div className="card-head"><h2>Template</h2></div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.9rem" }}>
            <div className="field">
              <label htmlFor="template">Sequence to amplify (A/C/G/T)</label>
              <textarea
                id="template"
                value={template}
                onChange={(e) => setTemplate(e.target.value)}
                rows={7}
                placeholder="Paste the gene or region you want to amplify"
                aria-describedby="template-count"
              />
              <span id="template-count" className="label">{cleanLength} nt entered</span>
            </div>

            <div className="field">
              <span className="field-label">Reaction</span>
              <div className="mode-list">
                <button
                  type="button"
                  className={`mode ${mode === "conventional" ? "on" : ""}`}
                  onClick={() => setMode("conventional")}
                  aria-pressed={mode === "conventional"}
                >
                  <strong>Conventional</strong>
                  <span>Copy the region. The product is what was there.</span>
                </button>
                <button
                  type="button"
                  className={`mode ${mode === "cloning" ? "on" : ""}`}
                  onClick={() => setMode("cloning")}
                  aria-pressed={mode === "cloning"}
                >
                  <strong>Cloning</strong>
                  <span>Add a site to each end, then cut the product into an insert.</span>
                </button>
              </div>
            </div>

            {mode === "cloning" && (
              <>
                <div className="row-2">
                  <div className="field">
                    <label htmlFor="left-enzyme">5&prime; enzyme</label>
                    <select id="left-enzyme" value={leftEnzyme} onChange={(e) => setLeftEnzyme(e.target.value)}>
                      {enzymes.map((e) => (
                        <option key={e.name} value={e.name}>{e.name} &middot; {e.recognition}</option>
                      ))}
                    </select>
                  </div>
                  <div className="field">
                    <label htmlFor="right-enzyme">3&prime; enzyme</label>
                    <select id="right-enzyme" value={rightEnzyme} onChange={(e) => setRightEnzyme(e.target.value)}>
                      {enzymes.map((e) => (
                        <option key={e.name} value={e.name}>{e.name} &middot; {e.recognition}</option>
                      ))}
                    </select>
                  </div>
                </div>

                <div className="field">
                  <label htmlFor="clamp">Clamp bases outside each site</label>
                  <input
                    id="clamp" type="number" min={0} max={20} value={clamp}
                    onChange={(e) => setClamp(Number(e.target.value))}
                    aria-describedby="clamp-note"
                  />
                  <span id="clamp-note" className="note">
                    Enzymes cut poorly at a fragment&rsquo;s end. Six is the usual minimum.
                  </span>
                </div>

                <div className="checks">
                  <label>
                    <input type="checkbox" checked={keepFrame} onChange={(e) => setKeepFrame(e.target.checked)} />
                    Keep the vector&rsquo;s reading frame
                  </label>
                </div>
              </>
            )}
          </div>
        </div>

        <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem", minWidth: 0 }}>
          {error && <div className="notice notice-error" role="alert">{error}</div>}

          {!result && !error && (
            <div className="card">
              <div className="empty">
                <Icon name="helix" size={34} className="glyph" />
                <strong>No reaction yet</strong>
                <span>Paste a template and press Design PCR.</span>
              </div>
            </div>
          )}

          {result && (
            <>
              {result.problems.length > 0 && (
                <div className="notice notice-error" role="alert">
                  <strong>This will not give you an insert.</strong>
                  <ul style={{ margin: "0.45rem 0 0", paddingLeft: "1.1rem" }}>
                    {result.problems.map((p) => <li key={p}>{p}</li>)}
                  </ul>
                </div>
              )}

              <div className="card">
                <div className="card-head">
                  <h2 style={{ flex: 1 }}>Primers to order</h2>
                  <span className="terminal-end label">
                    anneal at <strong>{result.annealing_temperature.toFixed(1)} °C</strong>
                  </span>
                </div>
                <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.8rem" }}>
                  <PrimerRow primer={result.forward} label="forward" />
                  <PrimerRow primer={result.reverse} label="reverse" />
                  {result.forward.tail && (
                    <p className="note">
                      The annealing temperature comes from the annealing part alone. In
                      the first cycle the tail has nothing to pair with, so setting it
                      from the whole oligo&rsquo;s Tm runs far too hot and nothing
                      amplifies.
                    </p>
                  )}
                </div>
              </div>

              <div className="card">
                <div className="card-head">
                  <h2 style={{ flex: 1 }}>Product</h2>
                  <span className="note nums">{result.product_length} bp</span>
                </div>
                <div className="card-body">
                  <ProductView result={result} />
                  {result.forward.tail && (
                    <div className="duplex .keys keys" style={{ marginTop: "0.6rem" }}>
                      <span className="key"><i style={{ background: "var(--amber)" }} /> added by the primers</span>
                      <span className="key"><i style={{ background: "var(--accent)" }} /> copied from the template</span>
                    </div>
                  )}
                </div>
              </div>

              {result.digest && (
                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>
                      Cut with {result.left_enzyme} and {result.right_enzyme}
                    </h2>
                    <span className="note nums">{result.digest.length} bp</span>
                  </div>
                  <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.8rem" }}>
                    <div className="stat-row">
                      <div className="stat">
                        <div className="k">left end</div>
                        <div className="v mono">
                          {result.digest.left_end.sequence || "blunt"}
                          <small> {result.digest.left_end.kind}</small>
                        </div>
                      </div>
                      <div className="stat">
                        <div className="k">right end</div>
                        <div className="v mono">
                          {result.digest.right_end.sequence || "blunt"}
                          <small> {result.digest.right_end.kind}</small>
                        </div>
                      </div>
                      <div className="stat">
                        <div className="k">trimmed away</div>
                        <div className="v nums">
                          {result.digest.trimmed_left + result.digest.trimmed_right}
                          <small> bp</small>
                        </div>
                      </div>
                    </div>

                    <p className="note">
                      These ends were read off the cut molecule, not looked up from the
                      enzyme table &mdash; a value copied from the table agrees with the
                      table whatever the bases actually spell.
                    </p>

                    <div className="seq-block">{result.digest.top}</div>

                    <div>
                      <button className="btn btn-primary" onClick={sendToClone}>
                        Clone into a vector <Icon name="arrowRight" size={16} />
                      </button>
                    </div>
                  </div>
                </div>
              )}

              {result.warnings.length > 0 && (
                <div className="card">
                  <div className="card-head"><h2>Notes</h2></div>
                  <div className="card-body">
                    <ul style={{ margin: 0, paddingLeft: "1.1rem", color: "var(--ink-soft)" }}>
                      {[...new Set(result.warnings)].map((w) => (
                        <li key={w} style={{ marginBottom: "0.35rem", lineHeight: 1.5 }}>{w}</li>
                      ))}
                    </ul>
                  </div>
                </div>
              )}
            </>
          )}
        </div>
      </div>
    </>
  );
}
