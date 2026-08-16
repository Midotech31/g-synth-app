import { useCallback, useEffect, useState } from "react";
import { useLocation } from "react-router-dom";

import {
  ApiError,
  api,
  type AssemblyResult,
  type Catalogue,
  type DesignParams,
} from "../api/client";
import DuplexView from "../components/DuplexView";
import InsertForm from "../components/InsertForm";
import { segmentColour } from "../components/segmentColour";
import Icon from "../components/Icon";
import LiveStatus from "../components/LiveStatus";
import { useWorkspaceState } from "../state/WorkspaceStateContext";

const SAMPLE = "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAA";

const DEFAULTS: DesignParams = {
  sequence: SAMPLE,
  name: "construct",
  left_enzyme: "NdeI",
  right_enzyme: "XhoI",
  is_coding: false,
  remove_stop: false,
  cleavage_site: "Thrombin",
  include_his_tag: true,
  include_linkers: true,
  target_oligo_length: 90,
  overhang_length: 4,
};

export default function Design() {
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [params, setParams, clearParams] = useWorkspaceState<DesignParams>("design.params", DEFAULTS);
  const [result, setResult, clearResult] = useWorkspaceState<AssemblyResult | null>("design.result", null);
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const [saved, setSaved, clearSaved] = useWorkspaceState("design.saved", "");
  const [exportOpen, setExportOpen] = useState(false);
  const [copied, setCopied] = useState(false);
  const location = useLocation();

  // The optimiser hands its gene over rather than making the user copy it.
  useEffect(() => {
    const handed = (location.state as { sequence?: string } | null)?.sequence;
    if (handed) {
      setParams((current) => ({ ...current, sequence: handed }));
      setResult(null);
    }
  }, [location.state]);

  useEffect(() => {
    api.catalogue().then(setCatalogue).catch(() => {
      setError("Could not load the enzyme catalogue.");
    });
  }, []);

  const set = useCallback(
    <K extends keyof DesignParams>(key: K, value: DesignParams[K]) => {
      setParams((current) => ({ ...current, [key]: value }));
      setResult(null);
      setSaved("");
    },
    [],
  );

  async function design(saveAsProject = false) {
    setBusy(true);
    setError("");
    setSaved("");
    try {
      const data = await api.designAssembly({ ...params, save_as_project: saveAsProject });
      setResult(data);
      if (data.project_id) setSaved(`Saved to your projects (#${data.project_id}).`);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "The design failed.");
      setResult(null);
    } finally {
      setBusy(false);
    }
  }

  async function download(kind: "order-sheet" | "protocol") {
    const suffix = kind === "order-sheet" ? "oligos.csv" : "protocol.txt";
    try {
      await api.download(
        `/api/design/assembly/${kind}/`,
        params,
        `${(params.name || "construct").replace(/\s+/g, "_")}_${suffix}`,
      );
    } catch {
      setError("The download failed. Try designing again first.");
    }
  }

  /** Take the construct out: GenBank for a viewer, FASTA for a supplier. */
  async function exportConstruct(filetype: "genbank" | "fasta" | "oligos") {
    const safe = (params.name || "construct").replace(/\s+/g, "_");
    const names = { genbank: `${safe}.gb`, fasta: `${safe}.fasta`,
                    oligos: `${safe}_oligos.fasta` };
    try {
      await api.download(
        `/api/design/assembly/export/?filetype=${filetype}`, params, names[filetype],
      );
    } catch {
      setError("The download failed. Try designing again first.");
    }
  }

  const verified = result !== null && result.verification.length === 0;

  function clearWorkspace() {
    clearParams();
    clearResult();
    clearSaved();
    setError("");
    setExportOpen(false);
    setCopied(false);
  }

  async function copyConstruct() {
    if (!result) return;
    await navigator.clipboard.writeText(result.construct_forward);
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1800);
  }

  // The verdict is the reason for the wait, so it is what gets said. The
  // error has its own alert; repeating it here would announce it twice.
  const status = busy
    ? "Designing…"
    : result === null
      ? ""
      : verified
        ? `Design verified: ${result.fragment_count} fragments, ${result.oligo_count} oligos to order.`
        : "Design failed verification. Do not order these oligos.";

  return (
    <>
      <LiveStatus message={status} />

      <div className="topbar design-topbar">
        <div className="grow design-heading">
          <h1>Design a construct</h1>
          <div className="design-steps" aria-label="Design progress">
            <div className="design-step active"><span>1</span>Insert</div>
            <i aria-hidden="true" />
            <div className={`design-step ${result ? "complete" : ""}`}><span>2</span>Configure</div>
            <i aria-hidden="true" />
            <div className={`design-step ${result ? "complete" : ""}`}><span>3</span>Review</div>
          </div>
        </div>
        <button
          className="btn btn-primary design-save"
          onClick={() => design(true)}
          disabled={busy || !verified}
        >
          <Icon name="check" size={18} />
          Save project
        </button>
        <div className="design-export">
          <button
            className="btn btn-outline"
            onClick={() => setExportOpen((open) => !open)}
            disabled={!verified}
            aria-expanded={exportOpen}
            aria-controls="design-export-menu"
          >
            <Icon name="arrowRight" size={17} />
            Export
          </button>
          {exportOpen && (
            <div
              className="design-export-menu"
              id="design-export-menu"
              role="menu"
              aria-label="Export formats"
            >
              <button role="menuitem" onClick={() => void download("order-sheet")}>Oligo CSV</button>
              <button role="menuitem" onClick={() => void exportConstruct("oligos")}>Oligo FASTA</button>
              <button role="menuitem" onClick={() => void download("protocol")}>Protocol</button>
              <button role="menuitem" onClick={() => void exportConstruct("genbank")}>GenBank</button>
            </div>
          )}
        </div>
        <button className="btn btn-outline" onClick={() => void copyConstruct()} disabled={!result}>
          <Icon name="book" size={17} />
          {copied ? "Copied" : "Copy"}
        </button>
      </div>

      <div
        className="content design-content"
        style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}
        aria-busy={busy}
      >
        {error && <div className="notice notice-error" role="alert">{error}</div>}
        {saved && <div className="notice notice-info" role="status">{saved}</div>}

        <div className="design-layout">
          {/* ── Inputs ─────────────────────────────────────────────────── */}
          <div className="card design-input-card">
            <div className="card-head"><h2>Insert</h2></div>
            <div className="card-body">
              <InsertForm params={params} catalogue={catalogue} onChange={set} />
            </div>
            <div className="design-form-actions">
              <button className="btn btn-outline" onClick={clearWorkspace} disabled={busy}>
                Clear
              </button>
              <button className="btn btn-primary" onClick={() => design(false)} disabled={busy}>
                {busy && <span className="spinner" />}
                {busy ? "Designing…" : result ? "Update design" : "Design"}
              </button>
            </div>
          </div>

          {/* ── Results ────────────────────────────────────────────────── */}
          <div className="design-results">
            {!result ? (
              <div className="card">
                <div className="empty">
                  <Icon name="helix" size={38} className="glyph" />
                  <strong>No design yet</strong>
                  <span>Set the insert and its ends, then press Design.</span>
                </div>
              </div>
            ) : (
              <>
                <div className={`design-verdict ${verified ? "verified" : "failed"}`}>
                  <div className="design-verdict-icon" aria-hidden="true">
                    <Icon name={verified ? "check" : "cross"} size={28} />
                  </div>
                  {verified ? (
                    <div className="design-verdict-copy">
                      <strong>Verified</strong>
                      <span>Annealing these oligo pairs and ligating them in order reproduces the construct exactly, on both strands.</span>
                    </div>
                  ) : (
                    <div className="design-verdict-copy">
                      <strong>Do not order</strong>
                      <span>{result.verification.join(" ")}</span>
                    </div>
                  )}
                  {verified && (
                    <div className="design-verdict-date">
                      <span>Verified on</span>
                      <strong>{new Intl.DateTimeFormat("en-GB", {
                        day: "2-digit", month: "short", year: "numeric",
                      }).format(new Date())}</strong>
                    </div>
                  )}
                </div>

                <div className="card design-stats">
                  <div className="card-body stat-row">
                    <div className="stat">
                      <Icon name="helix" size={19} />
                      <div className="k">Construct</div>
                      <div className="v">{result.construct_length}<small>bp</small></div>
                    </div>
                    <div className="stat">
                      <Icon name="target" size={19} />
                      <div className="k">GC</div>
                      <div className="v">{result.construct_gc}<small>%</small></div>
                    </div>
                    <div className="stat">
                      <Icon name="plate" size={19} />
                      <div className="k">Fragments</div>
                      <div className="v">{result.fragment_count}</div>
                    </div>
                    <div className="stat">
                      <Icon name="book" size={19} />
                      <div className="k">Oligos</div>
                      <div className="v">{result.oligo_count}</div>
                    </div>
                    <div className="stat">
                      <Icon name="scales" size={19} />
                      <div className="k">Longest oligo</div>
                      <div className="v">{result.longest_oligo}<small>nt</small></div>
                    </div>
                    {/* A long gene needs more distinct junctions than 4 nt can
                        supply, so the design widens them. The form still shows
                        what was asked for; this shows what was built. */}
                    <div className="stat">
                      <Icon name="target" size={19} />
                      <div className="k">Overhang</div>
                      <div className="v">
                        {result.overhang_length}<small>nt</small>
                        {result.overhang_length !== params.overhang_length && (
                          <small className="widened">widened</small>
                        )}
                      </div>
                    </div>
                  </div>
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Construct map</h2>
                    {/* Measured off the assembled fragments. The design's own
                        label cannot disagree with itself, so showing that
                        would confirm nothing. */}
                    <span className="label">
                      {result.terminal_ends.map((end) => (
                        <span key={end.side} className="terminal-end">
                          {end.enzyme}{" "}
                          <strong>{end.overhang || "blunt"}</strong>
                          {end.overhang ? ` ${end.kind}` : ""}
                        </span>
                      ))}
                    </span>
                  </div>
                  <div className="card-body">
                    <div className="track">
                      {result.ssd.segments.map((segment, index) => (
                        <div
                          key={`${segment.name}-${index}`}
                          className="track-part"
                          style={{
                            flexGrow: segment.sequence.length,
                            background: segmentColour(segment.name),
                          }}
                          title={`${segment.name} · ${segment.start + 1}–${segment.end} (${segment.sequence.length} nt)`}
                        >
                          {segment.sequence.length > 12 ? segment.name : ""}
                        </div>
                      ))}
                    </div>
                    <div className="track junctions" aria-label="fragment boundaries">
                      {result.fragments.map((fragment) => (
                        <div
                          key={fragment.index}
                          className="track-part frag"
                          style={{ flexGrow: fragment.forward_length }}
                          title={`${fragment.name}: ${fragment.top_start + 1}–${fragment.top_end}`}
                        >
                          {fragment.name}
                        </div>
                      ))}
                    </div>
                    {result.junction_overhangs.length > 0 && (
                      <p className="label" style={{ marginTop: "0.6rem" }}>
                        Junctions: {result.junction_overhangs.map((o) => `5'-${o}`).join(" · ")}
                      </p>
                    )}
                  </div>
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Hybridisation</h2>
                    <span className="label">
                      {result.duplex.mismatches.length === 0
                        ? "no mismatches"
                        : `${result.duplex.mismatches.length} mismatches`}
                    </span>
                  </div>
                  <div className="card-body">
                    <DuplexView duplex={result.duplex} />
                  </div>
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Oligos to order</h2>
                    <span className="label">Use Export for supplier files and protocol</span>
                  </div>
                  <div className="table-scroll">
                    <table className="data">
                      <thead>
                        <tr>
                          <th>Name</th><th>Sequence (5'→3')</th>
                          <th>Length</th>
                          <th title={`${result.tm_conditions.model} · ${result.tm_conditions.summary}`}>
                            Tm
                          </th>
                          <th>Scale</th><th>Purification</th>
                        </tr>
                      </thead>
                      <tbody>
                        {result.oligos.map((oligo) => (
                          <tr key={String(oligo.Name)}>
                            <td className="mono">{oligo.Name}</td>
                            <td className="mono seq-cell">{oligo["Sequence (5'->3')"]}</td>
                            <td className="num">{oligo["Length (nt)"]}</td>
                            <td className="num">{oligo["Tm (°C)"]}</td>
                            <td>{oligo.Scale}</td>
                            <td>{oligo.Purification}</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                  <div className="card-body" style={{ paddingTop: 0 }}>
                    <p className="note" style={{ margin: 0 }}>
                      Tm from the {result.tm_conditions.model} model, under the
                      conditions of the annealing step — {result.tm_conditions.summary}.
                    </p>
                  </div>
                </div>

                {(result.warnings.length > 0 || result.ssd.warnings.length > 0) && (
                  <div className="card">
                    <div className="card-head"><h2>Notes</h2></div>
                    <div className="card-body">
                      <ul style={{ margin: 0, paddingLeft: "1.1rem", color: "var(--ink-soft)" }}>
                        {[...new Set([...result.ssd.warnings, ...result.warnings])].map((note) => (
                          <li key={note} style={{ marginBottom: "0.3rem" }}>{note}</li>
                        ))}
                      </ul>
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

