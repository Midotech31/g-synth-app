import { useEffect, useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";

import {
  ApiError,
  api,
  type CodonHost,
  type OptimiseParams,
  type OptimiseResult,
} from "../api/client";
import Icon from "../components/Icon";
import LiveStatus from "../components/LiveStatus";
import PreflightPanel from "../components/PreflightPanel";
import { useWorkspaceState } from "../state/WorkspaceStateContext";

const SAMPLE =
  "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGAAATAATGGAGCACATATGGGA" +
  "TTAAATTTAGCATTATTAGGATTAGCAAGTTTATTAGGTAAAGGTATTAGTAAATTAGGA";

const DEFAULTS: OptimiseParams = {
  sequence: SAMPLE,
  host: "ecoli",
  is_protein: false,
  protein_context: "auto",
  keep_stop: false,
  avoid_enzymes: ["NdeI", "XhoI"],
  avoid_motifs: [],
  max_homopolymer: 5,
  gc_min: 40,
  gc_max: 60,
  gc_window: 50,
  max_repeat: 15,
  avoid_rare: true,
};


const COMMON = ["NdeI", "XhoI", "BamHI", "EcoRI", "HindIII", "NotI", "SalI", "SacI", "XbaI", "NcoI"];

export default function Optimise() {
  const [params, setParams, clearParams] = useWorkspaceState<OptimiseParams>("optimise.params", DEFAULTS);
  const [result, setResult, clearResult] = useWorkspaceState<OptimiseResult | null>("optimise.result", null);
  const [error, setError] = useState("");
  const [catalogueError, setCatalogueError] = useState("");
  const [busy, setBusy] = useState(false);
  const [hosts, setHosts] = useState<CodonHost[]>([]);
  const navigate = useNavigate();

  useEffect(() => {
    let cancelled = false;
    api.codonHosts()
      .then((catalogue) => {
        if (!cancelled) {
          setHosts(catalogue.hosts);
          setCatalogueError("");
        }
      })
      .catch(() => {
        if (!cancelled) setCatalogueError("Host profiles unavailable. Check the API and reload.");
      });
    return () => { cancelled = true; };
  }, []);

  const selectedHost = useMemo(
    () => hosts.find((host) => host.key === (params.host ?? "ecoli")),
    [hosts, params.host],
  );
  const hostGroups = useMemo(() => {
    const groups = new Map<string, CodonHost[]>();
    hosts.forEach((host) => {
      const group = groups.get(host.category) ?? [];
      group.push(host);
      groups.set(host.category, group);
    });
    return Array.from(groups.entries());
  }, [hosts]);

  function set<K extends keyof OptimiseParams>(key: K, value: OptimiseParams[K]) {
    setParams((current) => ({ ...current, [key]: value }));
    setResult(null);
  }

  function toggleEnzyme(name: string) {
    const current = params.avoid_enzymes ?? [];
    set(
      "avoid_enzymes",
      current.includes(name) ? current.filter((e) => e !== name) : [...current, name],
    );
  }

  async function run() {
    setBusy(true);
    setError("");
    try {
      setResult(await api.optimise(params));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "The optimisation failed.");
      setResult(null);
    } finally {
      setBusy(false);
    }
  }

  function sendToDesign() {
    if (!result) return;
    navigate("/design", {
      state: {
        sequence: result.sequence,
        isCoding: result.recommended_design_is_coding,
      },
    });
  }

  const inputLength = params.sequence.replace(/\s/g, "").length;
  const peptideStartsWithMet = params.sequence.replace(/\s/g, "").toUpperCase().startsWith("M");
  const peptideContext = params.protein_context ?? "auto";
  const peptideDecision = peptideContext === "complete_orf"
    ? peptideStartsWithMet
      ? "One N-terminal Met will be preserved."
      : "One initiator Met (ATG) will be added."
    : peptideContext === "mature_peptide"
      ? "The peptide will be preserved exactly; Design will supply translation initiation upstream."
      : peptideStartsWithMet
        ? "Detected N-terminal Met: treat as a complete ORF."
        : "No N-terminal Met: preserve as a mature peptide."

  function clearWorkspace() {
    clearParams();
    clearResult();
    setError("");
  }

  const status = busy
    ? "Optimising…"
    : result === null
      ? ""
      : result.is_clean
        ? `Optimised: ${result.changed_codons} codons changed, ${result.metric_label} ${result.cai_after}, GC ${result.gc_after}%. Translation verified.`
        : "Optimised, but the sequence cannot be made clean. Read the problems above the result.";

  return (
    <>
      <LiveStatus message={status} />

      <div className="topbar">
        <div className="grow">
          <h1>Host optimisation</h1>
          <p className="sub">
            Optimise coding DNA or back-translate a peptide for its expression host.
          </p>
        </div>
        <button className="btn btn-outline" onClick={clearWorkspace} disabled={busy}>
          Clear
        </button>
        <button className="btn btn-primary" onClick={() => void run()} disabled={busy || !inputLength}>
          {busy && <span className="spinner" />}
          {busy
            ? "Optimising…"
            : params.is_protein
              ? "Back-translate & optimise"
              : "Optimise DNA"}
        </button>
      </div>

      <div
        className="content"
        style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}
        aria-busy={busy}
      >
        {error && <div className="notice notice-error" role="alert">{error}</div>}

        <div className="design-layout">
          <div className="card">
            <div className="card-head"><h2>Input</h2></div>
            <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.85rem" }}>
              <div className="field">
                <label htmlFor="opt-host">Expression host</label>
                <select
                  id="opt-host"
                  value={params.host ?? "ecoli"}
                  onChange={(event) => set("host", event.target.value)}
                  disabled={hosts.length === 0}
                >
                  {hosts.length === 0 && <option value="ecoli">Loading evidence-backed profiles…</option>}
                  {hostGroups.map(([category, entries]) => (
                    <optgroup key={category} label={category}>
                      {entries.map((host) => (
                        <option key={host.key} value={host.key}>{host.name}</option>
                      ))}
                    </optgroup>
                  ))}
                </select>
                {catalogueError && <div className="notice notice-error compact" role="alert">{catalogueError}</div>}
                {selectedHost && (
                  <>
                    <p className="note" style={{ marginTop: "0.4rem" }}>
                      {selectedHost.dataset} · {selectedHost.coding_sequences.toLocaleString()} CDS · {selectedHost.gc_percent.toFixed(2)}% GC
                    </p>
                    <details className="advanced-control" style={{ marginTop: "0.45rem" }}>
                      <summary>Scientific provenance</summary>
                      <div style={{ padding: "0 0.75rem 0.75rem" }}>
                        <p className="note">
                          HIVE-CUTs {selectedHost.dataset_release} · NCBI taxon {selectedHost.taxon_id} · {selectedHost.codon_count.toLocaleString()} codons. {selectedHost.data_scope}.
                        </p>
                        <p className="note" style={{ marginTop: "0.35rem" }}>
                          The score is profile-relative and does not predict expression yield. <a href={selectedHost.source_url} target="_blank" rel="noreferrer">Primary data</a>
                        </p>
                      </div>
                    </details>
                  </>
                )}
              </div>

              <details>
                <summary>Use a strain-, cell- or tissue-specific reference set</summary>
                <div className="field" style={{ marginTop: "0.7rem" }}>
                  <label htmlFor="opt-reference-genes">Highly expressed reference CDSs</label>
                  <textarea
                    id="opt-reference-genes"
                    className="mono"
                    rows={4}
                    value={(params.reference_genes ?? []).join("\n")}
                    onChange={(event) => set(
                      "reference_genes",
                      event.target.value.split(/\r?\n/).map((entry) => entry.trim()).filter(Boolean),
                    )}
                    placeholder="One complete coding sequence per line"
                  />
                  <p className="note" style={{ marginTop: "0.4rem" }}>
                    Overrides the species profile for context-specific CAI.
                  </p>
                </div>
              </details>

              <div className="field">
                <span className="field-label" id="input-molecule-label">Input molecule</span>
                <div className="mode-switch" role="group" aria-labelledby="input-molecule-label">
                  <button
                    type="button"
                    className={!params.is_protein ? "active" : ""}
                    aria-pressed={!params.is_protein}
                    onClick={() => set("is_protein", false)}
                  >DNA</button>
                  <button
                    type="button"
                    className={params.is_protein ? "active" : ""}
                    aria-pressed={params.is_protein}
                    onClick={() => set("is_protein", true)}
                  >Peptide</button>
                </div>
              </div>

              {params.is_protein && (
                <div className="field">
                  <label htmlFor="protein-context">Peptide role</label>
                  <select
                    id="protein-context"
                    value={peptideContext}
                    onChange={(event) => set(
                      "protein_context",
                      event.target.value as NonNullable<OptimiseParams["protein_context"]>,
                    )}
                  >
                    <option value="auto">Automatic from N-terminal Met</option>
                    <option value="mature_peptide">Mature peptide / fusion insert</option>
                    <option value="complete_orf">Complete expression protein</option>
                  </select>
                  <p className="note" style={{ marginTop: "0.4rem" }}>{peptideDecision}</p>
                </div>
              )}

              <div className="field">
                <label htmlFor="opt-seq">
                  {params.is_protein ? "Peptide (one-letter code)" : "Coding DNA (A/C/G/T)"}
                </label>
                <textarea
                  id="opt-seq"
                  value={params.sequence}
                  onChange={(e) => set("sequence", e.target.value)}
                  rows={7}
                  className="mono"
                  style={{ fontSize: "0.78rem" }}
                  aria-describedby="opt-seq-count"
                />
                <span className="label" id="opt-seq-count">
                  {inputLength} {params.is_protein ? "residues" : "nt"} entered
                </span>
              </div>

              <div className="checks">
                <label>
                  <input type="checkbox" checked={params.keep_stop}
                         onChange={(e) => set("keep_stop", e.target.checked)} />
                  Add a stop codon
                </label>
                <label>
                  <input type="checkbox" checked={params.avoid_rare}
                         onChange={(e) => set("avoid_rare", e.target.checked)} />
                  Avoid low-frequency codons in this profile
                </label>
              </div>

              <div className="field">
                <span className="field-label" id="avoid-label">Keep these sites out</span>

                <div className="enzyme-chips" role="group" aria-labelledby="avoid-label">
                  {COMMON.map((name) => (
                    <button
                      key={name}
                      type="button"
                      className={params.avoid_enzymes?.includes(name) ? "chip on" : "chip"}
                      aria-pressed={params.avoid_enzymes?.includes(name) ?? false}
                      onClick={() => toggleEnzyme(name)}
                    >
                      {name}
                    </button>
                  ))}
                </div>

                <p className="note" style={{ marginTop: "0.4rem" }}>
                  Select restriction sites that must be absent from the gene.
                </p>
              </div>

              <div className="row-2">
                <div className="field">
                  <label htmlFor="gc-min">GC floor (%)</label>
                  <input id="gc-min" type="number" min={10} max={90}
                         value={params.gc_min}
                         onChange={(e) => set("gc_min", Number(e.target.value))} />
                </div>
                <div className="field">
                  <label htmlFor="gc-max">GC ceiling (%)</label>
                  <input id="gc-max" type="number" min={10} max={90}
                         value={params.gc_max}
                         onChange={(e) => set("gc_max", Number(e.target.value))} />
                </div>
              </div>

              <div className="row-2">
                <div className="field">
                  <label htmlFor="homo">Longest single-base run</label>
                  <input id="homo" type="number" min={3} max={12}
                         value={params.max_homopolymer}
                         onChange={(e) => set("max_homopolymer", Number(e.target.value))} />
                </div>
                <div className="field">
                  <label htmlFor="rep">Longest repeat (nt)</label>
                  <input id="rep" type="number" min={8} max={40}
                         value={params.max_repeat}
                         onChange={(e) => set("max_repeat", Number(e.target.value))} />
                </div>
              </div>
            </div>
          </div>

          <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
            {!result ? (
              <div className="card">
                <div className="empty">
                  <Icon name="helix" size={38} className="glyph" />
                  <strong>No result yet</strong>
                  <span>Enter DNA or peptide, then run optimisation.</span>
                </div>
              </div>
            ) : (
              <>
                <div className={`notice ${result.is_clean ? "notice-ok" : "notice-error"}`}>
                  {result.is_clean ? (
                    <>
                      <strong>Clean.</strong> Translation verified; excluded sites absent.
                    </>
                  ) : (
                    <>
                      <strong>Cannot be made clean.</strong>{" "}
                      {result.problems.join(" ")}
                    </>
                  )}
                </div>

                {result.initiator_methionine_added && (
                  <div className="notice notice-info compact">
                    Initiator Met added for expression; the supplied mature peptide remains recorded unchanged.
                  </div>
                )}

                <PreflightPanel report={result.preflight} />

                <div className="card">
                  <div className="card-body stat-row">
                    <div className="stat">
                      <div className="k">{result.metric_label}</div>
                      <div className="v">
                        {result.cai_after}
                        {result.cai_before !== null && <small>from {result.cai_before}</small>}
                      </div>
                    </div>
                    <div className="stat">
                      <div className="k">GC</div>
                      <div className="v">
                        {result.gc_after}<small>%</small>
                        {result.gc_before !== null && <small>from {result.gc_before}%</small>}
                      </div>
                    </div>
                    <div className="stat">
                      <div className="k">Low-frequency codons</div>
                      <div className="v">
                        {result.rare_codons_after}
                        <small>from {result.rare_codons_before}</small>
                      </div>
                    </div>
                    <div className="stat">
                      <div className="k">Codons changed</div>
                      <div className="v">{result.changed_codons}</div>
                    </div>
                    <div className="stat">
                      <div className="k">Length</div>
                      <div className="v">{result.length}<small>nt</small></div>
                    </div>
                  </div>
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Optimised gene</h2>
                    <button className="btn btn-outline"
                            onClick={() => void navigator.clipboard?.writeText(result.sequence)}>
                      Copy
                    </button>
                    <button className="btn btn-primary" onClick={sendToDesign} disabled={result.preflight?.can_export === false}>
                      Design oligos →
                    </button>
                  </div>
                  <div className="card-body">
                    <div className="seq-block">{result.sequence}</div>
                    {result.sites_removed.length > 0 && (
                      <p className="note" style={{ marginTop: "0.6rem" }}>
                        Removed from the gene: {result.sites_removed.join(", ")}.
                      </p>
                    )}
                    <details style={{ marginTop: "0.55rem" }}>
                      <summary>Optimisation provenance</summary>
                      <p className="note" style={{ marginTop: "0.4rem" }}>
                        {result.table} — {result.table_source}.
                      </p>
                    </details>
                  </div>
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Protein</h2>
                      <span className="label">
                        {result.protein.length} residues
                        {result.initiator_methionine_added ? " · includes initiator Met" : " · unchanged"}
                      </span>
                  </div>
                  <div className="card-body">
                    <div className="seq-block">{result.protein}</div>
                  </div>
                </div>

                {result.warnings.length > 0 && (
                  <div className="card">
                    <div className="card-head"><h2>Compromises</h2></div>
                    <div className="card-body">
                      <ul style={{ margin: 0, paddingLeft: "1.1rem", color: "var(--ink-soft)" }}>
                        {result.warnings.map((note) => (
                          <li key={note} style={{ marginBottom: "0.3rem", lineHeight: 1.5 }}>
                            {note}
                          </li>
                        ))}
                      </ul>
                      <p className="note" style={{ marginTop: "0.6rem" }}>
                        Non-blocking trade-offs; blocking issues appear above.
                      </p>
                    </div>
                  </div>
                )}
              </>
            )}
          </div>
        </div>
      </div>
    </>
  );
}
