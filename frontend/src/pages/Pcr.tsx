import { useEffect, useMemo, useRef, useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";

import { ApiError, api, type Catalogue, type PcrResult } from "../api/client";
import Icon from "../components/Icon";
import GelSimulation from "../components/GelSimulation";
import EnzymePicker from "../components/EnzymePicker";
import PrimerAnnealingView, { PrimerTail } from "../components/PrimerAnnealingView";
import PreflightPanel from "../components/PreflightPanel";
import { useWorkspaceState } from "../state/WorkspaceStateContext";

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
        {primer.tail && <span className="seq-tail"><PrimerTail primer={primer} /></span>}
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

function DigestView({ result }: { result: PcrResult }) {
  if (!result.digest) return null;
  const { digest } = result;
  const leftLength = digest.left_end.sequence.length;
  const rightLength = digest.right_end.sequence.length;
  const top = `${digest.left_end.strand === "bottom" ? " ".repeat(leftLength) : ""}${digest.top}${digest.right_end.strand === "bottom" ? " ".repeat(rightLength) : ""}`;
  const bottomCore = [...digest.bottom].reverse().join("");
  const bottom = `${digest.left_end.strand === "top" ? " ".repeat(leftLength) : ""}${bottomCore}${digest.right_end.strand === "top" ? " ".repeat(rightLength) : ""}`;
  const width = Math.max(top.length, bottom.length);
  const columns = width <= 54
    ? [...Array(width).keys()]
    : [...Array(24).keys(), -1, ...Array(24).keys()].map((value, index) => (
      value === -1 ? -1 : index < 24 ? value : width - 24 + value
    ));
  const line = (sequence: string) => columns.map((position, index) => {
    if (position === -1) return <span key={`ellipsis-${index}`} className="digest-ellipsis">…</span>;
    const base = sequence[position] ?? " ";
    const partner = sequence === top ? bottom[position] : top[position];
    return (
      <span
        key={position}
        className={base !== " " && partner === " " ? "hybrid-overhang-base" : base === " " ? "dx-gap" : "dx-base"}
      >
        {base === " " ? "\u00a0" : base}
      </span>
    );
  });
  const pairs = columns.map((position, index) => (
    <span key={position === -1 ? `ellipsis-${index}` : position}>
      {position === -1 ? "…" : top[position] !== " " && bottom[position] !== " " ? "|" : " "}
    </span>
  ));
  return (
    <div
      className="pcr-digest-view"
      aria-label="Nucleotide-level digested insert and cohesive ends"
    >
      <div className="pcr-digest-head">
        <div>
          <span>Left · {result.left_enzyme}</span>
          <strong className="mono">{digest.left_end.sequence ? `${digest.left_end.kind.replace("'", "′")}-${digest.left_end.sequence}` : "blunt"}</strong>
        </div>
        <span className="pill pill-ok">Digestion verified</span>
        <div>
          <span>Right · {result.right_enzyme}</span>
          <strong className="mono">{digest.right_end.sequence ? `${digest.right_end.kind.replace("'", "′")}-${digest.right_end.sequence}` : "blunt"}</strong>
        </div>
      </div>
      <div className="duplex-scroll pcr-digest-duplex" tabIndex={0}>
        <div className="duplex-row">
          <span className="dx-end">5′</span>
          <span className="dx-seq">{line(top)}</span>
          <span className="dx-end">3′</span>
        </div>
        <div className="duplex-row">
          <span className="dx-end" />
          <span className="dx-seq dx-ticks">{pairs}</span>
          <span className="dx-end" />
        </div>
        <div className="duplex-row">
          <span className="dx-end">3′</span>
          <span className="dx-seq">{line(bottom)}</span>
          <span className="dx-end">5′</span>
        </div>
      </div>
      <div className="duplex keys">
        <span className="key"><i className="hybrid-key-overhang" /> exposed cohesive end</span>
        <span className="key"><code>|</code> Watson–Crick base pair</span>
      </div>
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
  const [experience, setExperience] = useWorkspaceState<"guided" | "expert">("pcr.experience", "guided");
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [customPrimers, setCustomPrimers, clearCustomPrimers] = useWorkspaceState(
    "pcr.customPrimers", false,
  );
  const [forwardPrimer, setForwardPrimer, clearForwardPrimer] = useWorkspaceState(
    "pcr.forwardPrimer", "",
  );
  const [reversePrimer, setReversePrimer, clearReversePrimer] = useWorkspaceState(
    "pcr.reversePrimer", "",
  );

  const [result, setResult, clearResult] = useWorkspaceState<PcrResult | null>("pcr.result", null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");
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
    clearCustomPrimers();
    clearForwardPrimer();
    clearReversePrimer();
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
        forward_primer: customPrimers ? forwardPrimer : null,
        reverse_primer: customPrimers ? reversePrimer : null,
      });
      if (requestVersion.current === version) {
        setResult(next);
        if (customPrimers) {
          setForwardPrimer(next.forward.sequence);
          setReversePrimer(next.reverse.sequence);
        }
      }
    } catch (err) {
      if (requestVersion.current === version) {
        setError(err instanceof ApiError ? err.message : "The design could not be run.");
      }
    } finally {
      setBusy(false);
    }
  }

  function sendToClone() {
    if (!result?.digest) return;
    navigate("/clone", {
      state: {
        preDigested: {
          top: result.digest.top,
          bottom: result.digest.bottom,
          leftEnzyme: result.left_enzyme,
          rightEnzyme: result.right_enzyme,
          orfStart: result.insert_orf_start,
          origin: "pcr",
        },
      },
    });
  }

  const enzymes = catalogue?.enzymes ?? [];
  const leftSuppliesStart = enzymes.find((enzyme) => enzyme.name === leftEnzyme)
    ?.supplies_start_codon ?? leftEnzyme === "NdeI";
  const customPrimerPairReady = forwardPrimer.replace(/[^A-Za-z]/g, "").length > 0
    && reversePrimer.replace(/[^A-Za-z]/g, "").length > 0;

  function editPrimers() {
    if (result) {
      setForwardPrimer(result.forward.sequence);
      setReversePrimer(result.reverse.sequence);
    }
    setCustomPrimers(true);
  }

  function useAutomaticPrimers() {
    requestVersion.current += 1;
    setCustomPrimers(false);
    setForwardPrimer("");
    setReversePrimer("");
    setResult(null);
    setError("");
  }

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
        <button
          className="btn btn-primary"
          onClick={() => void run()}
          disabled={busy || !cleanLength || (customPrimers && !customPrimerPairReady)}
        >
          {busy && <span className="spinner" />}
          {busy ? "Analysing…" : customPrimers ? "Revalidate primers" : "Design PCR"}
        </button>
      </div>

      <div className="content design-layout pcr-layout">
        <div className="card">
          <div className="card-head">
            <h2 style={{ flex: 1 }}>Template</h2>
            <div className="mode-switch" role="group" aria-label="PCR detail level">
              <button className={experience === "guided" ? "active" : ""} onClick={() => {
                setExperience("guided");
                setClamp(6);
                setKeepFrame(true);
                setStartCodonMode("use_site");
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
                    setCustomPrimers(false);
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
                    setCustomPrimers(false);
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
                  <EnzymePicker
                    id="left-enzyme"
                    label="5′ enzyme"
                    enzymes={enzymes}
                    value={leftEnzyme}
                    onChange={(value) => {
                      inputsChanged();
                      setLeftEnzyme(value);
                    }}
                  />
                  <EnzymePicker
                    id="right-enzyme"
                    label="3′ enzyme"
                    enzymes={enzymes}
                    value={rightEnzyme}
                    onChange={(value) => {
                      inputsChanged();
                      setRightEnzyme(value);
                    }}
                  />
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
                      preserves the validated reading frame and avoids a Met-Met start.
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

            {customPrimers && (
              <div className="manual-primer-editor" aria-label="Edit primer sequences">
                <div className="manual-primer-head">
                  <div>
                    <strong>Custom primers</strong>
                    <span>Enter complete oligos in 5′→3′ orientation.</span>
                  </div>
                  <button className="btn btn-outline" type="button" onClick={useAutomaticPrimers}>
                    Restore automatic design
                  </button>
                </div>
                <div className="manual-primer-grid">
                  <div className="field">
                    <label htmlFor="forward-primer">Forward primer (5′→3′)</label>
                    <textarea
                      id="forward-primer"
                      className="mono"
                      rows={3}
                      value={forwardPrimer}
                      onChange={(event) => {
                        inputsChanged();
                        setForwardPrimer(event.target.value);
                      }}
                    />
                  </div>
                  <div className="field">
                    <label htmlFor="reverse-primer">Reverse primer (5′→3′)</label>
                    <textarea
                      id="reverse-primer"
                      className="mono"
                      rows={3}
                      value={reversePrimer}
                      onChange={(event) => {
                        inputsChanged();
                        setReversePrimer(event.target.value);
                      }}
                    />
                  </div>
                </div>
                <span className="field-hint">
                  G-Synth will locate each 3′ annealing region, separate any 5′ addition,
                  and recompute Tm, specificity checks, product geometry and digestion.
                </span>
              </div>
            )}
          </div>
        </div>

        <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem", minWidth: 0 }}>
          {error && <div className="notice notice-error" role="alert">{error}</div>}

          {!result && !error && (
            <div className="card">
              <div className="empty">
                <Icon name="helix" size={34} className="glyph" />
                <strong>{customPrimers ? "Edited primers need validation" : "No reaction yet"}</strong>
                <span>{customPrimers ? "Press Revalidate primers before ordering or cloning." : "Paste a template and press Design PCR."}</span>
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

              <PreflightPanel report={result.preflight} />

              <div className="card">
                <div className="card-head">
                  <h2 style={{ flex: 1 }}>Primers to order</h2>
                  <span className={`pill ${result.primer_source === "custom" ? "pill-ok" : ""}`}>
                    {result.primer_source === "custom" ? "Edited · validated" : "G-Synth design"}
                  </span>
                  <button className="btn btn-outline" type="button" onClick={editPrimers}>
                    Edit primers
                  </button>
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
                  <h2 style={{ flex: 1 }}>Primer–template hybridization</h2>
                  <span className="label">cycle 1 geometry</span>
                </div>
                <div className="card-body">
                  <PrimerAnnealingView
                    forward={result.forward}
                    reverse={result.reverse}
                    templateLength={cleanLength}
                  />
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

              {result.gel && (
                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>{result.gel.title}</h2>
                    <span className="label">agarose gel simulation</span>
                  </div>
                  <div className="card-body">
                    <GelSimulation simulation={result.gel} />
                  </div>
                </div>
              )}

              {result.digest && (
                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>
                      Cut with {result.left_enzyme} and {result.right_enzyme}
                    </h2>
                    <span className="note nums">{result.digest.length} bp</span>
                  </div>
                  <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.8rem" }}>
                    <DigestView result={result} />
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
                      The 5′ primer additions are unpaired only in cycle 1. These cohesive
                      ends appear later, after the completed PCR product is digested.
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
