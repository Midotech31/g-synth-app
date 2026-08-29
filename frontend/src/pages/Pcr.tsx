import { useEffect, useMemo, useRef, useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";

import { ApiError, api, type Catalogue, type PcrResult } from "../api/client";
import Icon from "../components/Icon";
import PreflightPanel from "../components/PreflightPanel";
import { useWorkspaceState } from "../state/WorkspaceStateContext";

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

  const [template, setTemplate, clearTemplate] = useWorkspaceState("pcr.template", "");
  const [mode, setMode, clearMode] = useWorkspaceState<"conventional" | "cloning">("pcr.mode", "cloning");
  const [leftEnzyme, setLeftEnzyme, clearLeftEnzyme] = useWorkspaceState("pcr.leftEnzyme", "NdeI");
  const [rightEnzyme, setRightEnzyme, clearRightEnzyme] = useWorkspaceState("pcr.rightEnzyme", "XhoI");
  const [clamp, setClamp, clearClamp] = useWorkspaceState("pcr.clamp", 6);
  const [keepFrame, setKeepFrame, clearKeepFrame] = useWorkspaceState("pcr.keepFrame", true);
  const [startCodonMode, setStartCodonMode, clearStartCodonMode] = useWorkspaceState<"use_site" | "keep_both">("pcr.startCodonMode", "use_site");
  const [primerSizing, setPrimerSizing, clearPrimerSizing] = useWorkspaceState<"automatic" | "manual">("pcr.primerSizing", "automatic");
  const [forwardLength, setForwardLength, clearForwardLength] = useWorkspaceState("pcr.forwardLength", 20);
  const [reverseLength, setReverseLength, clearReverseLength] = useWorkspaceState("pcr.reverseLength", 20);
  const [experience, setExperience] = useWorkspaceState<"guided" | "expert">("pcr.experience", "guided");
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);

  const [result, setResult, clearResult] = useWorkspaceState<PcrResult | null>("pcr.result", null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");
  // Incremented whenever the inputs change or a new request starts. If an
  // older request finishes after an edit, its primers belong to the previous
  // form state and must not be put back on screen.
  const requestVersion = useRef(0);

  useEffect(() => {
    api.catalogue().then(setCatalogue).catch(() => setCatalogue(null));
  }, []);

  useEffect(() => {
    if (!location.state?.sequence) return;
    setTemplate(location.state.sequence);
    setResult(null);
  }, [location.state?.sequence, setResult, setTemplate]);

  const cleanLength = useMemo(
    () => template.replace(/[^A-Za-z]/g, "").length, [template],
  );

  function inputsChanged() {
    requestVersion.current += 1;
    setResult(null);
    setError("");
  }

  function clearWorkspace() {
    requestVersion.current += 1;
    clearTemplate();
    clearMode();
    clearLeftEnzyme();
    clearRightEnzyme();
    clearClamp();
    clearKeepFrame();
    clearStartCodonMode();
    clearPrimerSizing();
    clearForwardLength();
    clearReverseLength();
    clearResult();
    setError("");
  }

  async function run() {
    const version = ++requestVersion.current;
    setBusy(true);
    setError("");
    setResult(null);
    try {
      const next = await api.pcr({
        template,
        left_enzyme: mode === "cloning" ? leftEnzyme : null,
        right_enzyme: mode === "cloning" ? rightEnzyme : null,
        clamp,
        keep_frame: mode === "cloning" && keepFrame,
        start_codon_mode: startCodonMode,
        forward_anneal_length: primerSizing === "manual" ? forwardLength : null,
        reverse_anneal_length: primerSizing === "manual" ? reverseLength : null,
      });
      if (requestVersion.current === version) setResult(next);
    } catch (err) {
      if (requestVersion.current === version) {
        setError(err instanceof ApiError ? err.message : "The design could not be run.");
      }
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
  const leftSuppliesStart = enzymes.find((enzyme) => enzyme.name === leftEnzyme)
    ?.supplies_start_codon ?? leftEnzyme === "NdeI";

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
        <button className="btn btn-outline" onClick={clearWorkspace} disabled={busy}>
          Clear
        </button>
        <button className="btn btn-primary" onClick={() => void run()} disabled={busy || !cleanLength}>
          {busy && <span className="spinner" />}
          {busy ? "Designing…" : "Design PCR"}
        </button>
      </div>

      <div className="content design-layout">
        <div className="card">
          <div className="card-head">
            <h2 style={{ flex: 1 }}>Template</h2>
            <div className="mode-switch" role="group" aria-label="PCR detail level">
              <button className={experience === "guided" ? "active" : ""} onClick={() => {
                setExperience("guided");
                setClamp(6);
                setKeepFrame(true);
                setStartCodonMode("use_site");
                setPrimerSizing("automatic");
                inputsChanged();
              }}>Guided</button>
              <button className={experience === "expert" ? "active" : ""} onClick={() => setExperience("expert")}>Expert</button>
            </div>
          </div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.9rem" }}>
            <div className="field">
              <label htmlFor="template">Sequence to amplify (A/C/G/T)</label>
              <textarea
                id="template"
                value={template}
                onChange={(e) => {
                  inputsChanged();
                  setTemplate(e.target.value);
                }}
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
                  onClick={() => {
                    inputsChanged();
                    setMode("conventional");
                  }}
                  aria-pressed={mode === "conventional"}
                >
                  <strong>Conventional</strong>
                  <span>Copy the region. The product is what was there.</span>
                </button>
                <button
                  type="button"
                  className={`mode ${mode === "cloning" ? "on" : ""}`}
                  onClick={() => {
                    inputsChanged();
                    setMode("cloning");
                  }}
                  aria-pressed={mode === "cloning"}
                >
                  <strong>Cloning</strong>
                  <span>Add a site to each end, then cut the product into an insert.</span>
                </button>
              </div>
            </div>

            {mode === "cloning" && (
              <>
                {experience === "guided" && (
                  <div className="notice notice-info compact">
                    Uses a six-base terminal clamp, preserves the reading frame and avoids a duplicated start codon when the enzyme supplies ATG.
                  </div>
                )}
                <div className="row-2">
                  <div className="field">
                    <label htmlFor="left-enzyme">5&prime; enzyme</label>
                    <select
                      id="left-enzyme"
                      value={leftEnzyme}
                      onChange={(e) => {
                        inputsChanged();
                        setLeftEnzyme(e.target.value);
                      }}
                    >
                      {enzymes.map((e) => (
                        <option key={e.name} value={e.name}>{e.name} &middot; {e.recognition}</option>
                      ))}
                    </select>
                  </div>
                  <div className="field">
                    <label htmlFor="right-enzyme">3&prime; enzyme</label>
                    <select
                      id="right-enzyme"
                      value={rightEnzyme}
                      onChange={(e) => {
                        inputsChanged();
                        setRightEnzyme(e.target.value);
                      }}
                    >
                      {enzymes.map((e) => (
                        <option key={e.name} value={e.name}>{e.name} &middot; {e.recognition}</option>
                      ))}
                    </select>
                  </div>
                </div>

                {experience === "expert" && <div className="field">
                  <label htmlFor="clamp">Clamp bases outside each site</label>
                  <input
                    id="clamp" type="number" min={0} max={20} value={clamp}
                    onChange={(e) => {
                      inputsChanged();
                      setClamp(Number(e.target.value));
                    }}
                    aria-describedby="clamp-note"
                  />
                  <span id="clamp-note" className="note">
                    Enzymes cut poorly at a fragment&rsquo;s end. Six is the usual minimum.
                  </span>
                </div>}

                {experience === "expert" && leftSuppliesStart && (
                  <div className="field">
                    <label htmlFor="start-codon-mode">Start codon</label>
                    <select
                      id="start-codon-mode"
                      value={startCodonMode}
                      onChange={(e) => {
                        inputsChanged();
                        setStartCodonMode(e.target.value as "use_site" | "keep_both");
                      }}
                      aria-describedby="start-codon-note"
                    >
                      <option value="use_site">Use {leftEnzyme}&rsquo;s ATG (recommended)</option>
                      <option value="keep_both">Keep both ATGs (adds N-terminal Met)</option>
                    </select>
                    <span id="start-codon-note" className="note">
                      {leftEnzyme}&rsquo;s recognition site supplies ATG. Using it alone
                      matches the legacy G-Synth logic and avoids a Met-Met start.
                    </span>
                  </div>
                )}

                {experience === "expert" && <div className="checks">
                  <label>
                    <input
                      type="checkbox"
                      checked={keepFrame}
                      onChange={(e) => {
                        inputsChanged();
                        setKeepFrame(e.target.checked);
                      }}
                    />
                    Keep the vector&rsquo;s reading frame
                  </label>
                </div>}
              </>
            )}

            {experience === "expert" && <div className="field">
              <span className="field-label" id="primer-sizing-label">Annealing footprints</span>
              <div className="mode-list" role="group" aria-labelledby="primer-sizing-label">
                <button type="button" className={primerSizing === "automatic" ? "mode on" : "mode"}
                  aria-pressed={primerSizing === "automatic"}
                  onClick={() => { inputsChanged(); setPrimerSizing("automatic"); }}>
                  <strong>Automatic</strong><span>Choose 18–30 nt toward 60 °C.</span>
                </button>
                <button type="button" className={primerSizing === "manual" ? "mode on" : "mode"}
                  aria-pressed={primerSizing === "manual"}
                  onClick={() => { inputsChanged(); setPrimerSizing("manual"); }}>
                  <strong>Manual lengths</strong><span>Use exact 15–60 nt footprints and retain quality warnings.</span>
                </button>
              </div>
            </div>}

            {experience === "expert" && primerSizing === "manual" && <div className="row-2">
              <div className="field">
                <label htmlFor="forward-length">Forward annealing length</label>
                <input id="forward-length" type="number" min={15} max={60} value={forwardLength}
                  onChange={(event) => { inputsChanged(); setForwardLength(Number(event.target.value)); }} />
              </div>
              <div className="field">
                <label htmlFor="reverse-length">Reverse annealing length</label>
                <input id="reverse-length" type="number" min={15} max={60} value={reverseLength}
                  onChange={(event) => { inputsChanged(); setReverseLength(Number(event.target.value)); }} />
              </div>
            </div>}
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

              {result.alternative_pairs.length > 0 && (
                <div className="card">
                  <div className="card-head"><h2>Clean alternative enzyme pairs</h2></div>
                  <div className="card-body">
                    <p className="note" style={{ marginBottom: "0.7rem" }}>
                      Each pair passed the same complete PCR-product and digest simulation.
                      Confirm that both sites are unique in the vector you will actually use.
                    </p>
                    <div className="checks">
                      {result.alternative_pairs.map((pair) => (
                        <button type="button" className="btn btn-outline"
                          key={`${pair.left_enzyme}/${pair.right_enzyme}`}
                          onClick={() => {
                            setLeftEnzyme(pair.left_enzyme);
                            setRightEnzyme(pair.right_enzyme);
                            inputsChanged();
                          }}>
                          {pair.left_enzyme} ({pair.left_overhang || "blunt"}) / {pair.right_enzyme} ({pair.right_overhang || "blunt"})
                        </button>
                      ))}
                    </div>
                  </div>
                </div>
              )}

              <PreflightPanel report={result.preflight} />

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
                      <button className="btn btn-primary" onClick={sendToClone} disabled={result.preflight?.can_export === false}>
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
