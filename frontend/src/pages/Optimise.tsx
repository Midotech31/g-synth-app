import { useState } from "react";
import { useNavigate } from "react-router-dom";

import {
  ApiError,
  api,
  type OptimiseParams,
  type OptimiseResult,
} from "../api/client";
import Icon from "../components/Icon";

const SAMPLE =
  "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGAAATAATGGAGCACATATGGGA" +
  "TTAAATTTAGCATTATTAGGATTAGCAAGTTTATTAGGTAAAGGTATTAGTAAATTAGGA";

const DEFAULTS: OptimiseParams = {
  sequence: SAMPLE,
  is_protein: false,
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

/** Enzymes worth offering here: the ones a construct is usually cut with. */
const COMMON = ["NdeI", "XhoI", "BamHI", "EcoRI", "HindIII", "NotI", "SalI", "SacI", "XbaI", "NcoI"];

export default function Optimise() {
  const [params, setParams] = useState<OptimiseParams>(DEFAULTS);
  const [result, setResult] = useState<OptimiseResult | null>(null);
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const navigate = useNavigate();

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

  /** Hand the optimised gene to the design page, so the workflow continues. */
  function sendToDesign() {
    if (!result) return;
    navigate("/design", { state: { sequence: result.sequence } });
  }

  const inputLength = params.sequence.replace(/\s/g, "").length;

  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>Optimise for the host</h1>
          <p className="sub">
            Rewrite a gene for the organism that will express it. The protein
            never changes — everything else does.
          </p>
        </div>
        <button className="btn btn-primary" onClick={() => void run()} disabled={busy || !inputLength}>
          {busy && <span className="spinner" />}
          {busy ? "Optimising…" : "Optimise"}
        </button>
      </div>

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
        {error && <div className="notice notice-error">{error}</div>}

        <div className="design-layout">
          {/* ── Inputs ─────────────────────────────────────────────────── */}
          <div className="card">
            <div className="card-head"><h2>Gene</h2></div>
            <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.85rem" }}>
              <div className="field">
                <label htmlFor="opt-seq">
                  {params.is_protein ? "Protein (one letter)" : "Coding sequence (A/C/G/T)"}
                </label>
                <textarea
                  id="opt-seq"
                  value={params.sequence}
                  onChange={(e) => set("sequence", e.target.value)}
                  rows={7}
                  className="mono"
                  style={{ fontSize: "0.78rem" }}
                />
                <span className="label">
                  {inputLength} {params.is_protein ? "residues" : "nt"} entered
                </span>
              </div>

              <div className="checks">
                <label>
                  <input type="checkbox" checked={params.is_protein}
                         onChange={(e) => set("is_protein", e.target.checked)} />
                  This is a protein, not DNA
                </label>
                <label>
                  <input type="checkbox" checked={params.keep_stop}
                         onChange={(e) => set("keep_stop", e.target.checked)} />
                  Add a stop codon
                </label>
                <label>
                  <input type="checkbox" checked={params.avoid_rare}
                         onChange={(e) => set("avoid_rare", e.target.checked)} />
                  Avoid codons the host reads slowly
                </label>
              </div>

              <div className="field">
                <label>Keep these sites out</label>
                <div className="enzyme-chips">
                  {COMMON.map((name) => (
                    <button
                      key={name}
                      type="button"
                      className={params.avoid_enzymes?.includes(name) ? "chip on" : "chip"}
                      onClick={() => toggleEnzyme(name)}
                    >
                      {name}
                    </button>
                  ))}
                </div>
                <p className="note" style={{ marginTop: "0.4rem" }}>
                  Pick the pair you will clone with. A gene carrying an internal
                  NdeI site cannot be cloned NdeI/XhoI, however well it
                  translates.
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

          {/* ── Results ────────────────────────────────────────────────── */}
          <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
            {!result ? (
              <div className="card">
                <div className="empty">
                  <Icon name="helix" size={38} className="glyph" />
                  <strong>Nothing optimised yet</strong>
                  <span>Paste a gene or a protein, then press Optimise.</span>
                </div>
              </div>
            ) : (
              <>
                <div className={`notice ${result.is_clean ? "notice-ok" : "notice-error"}`}>
                  {result.is_clean ? (
                    <>
                      <strong>Clean.</strong> The protein is unchanged and the
                      sequence carries none of the sites you excluded.
                    </>
                  ) : (
                    <>
                      <strong>Cannot be made clean.</strong>{" "}
                      {result.problems.join(" ")}
                    </>
                  )}
                </div>

                <div className="card">
                  <div className="card-body stat-row">
                    <div className="stat">
                      <div className="k">CAI</div>
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
                      <div className="k">Slow codons</div>
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
                    <button className="btn btn-primary" onClick={sendToDesign}>
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
                    <p className="note" style={{ marginTop: "0.4rem" }}>
                      Codon usage: {result.table} — {result.table_source}.
                    </p>
                  </div>
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Protein</h2>
                    <span className="label">{result.protein.length} residues · unchanged</span>
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
                        These cost quality, not viability. A constraint that
                        would stop the construct working appears above instead.
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
