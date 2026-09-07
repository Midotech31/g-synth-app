import { useCallback, useEffect, useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";

import {
  ApiError,
  api,
  type AssemblyResult,
  type Catalogue,
  type DesignParams,
  type TerminalEnd,
} from "../api/client";
import CoreWorkflowTrail from "../components/CoreWorkflowTrail";
import EnzymePicker from "../components/EnzymePicker";
import DuplexView from "../components/DuplexView";
import InsertForm from "../components/InsertForm";
import { segmentColour } from "../components/segmentColour";
import Icon from "../components/Icon";
import LiveStatus from "../components/LiveStatus";
import PreflightPanel from "../components/PreflightPanel";
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

function assemblyTokenFor(result: AssemblyResult | null) {
  if (!result) return "";
  return result.provenance?.output_sha256
    ?? `${result.construct_length}:${result.construct_forward}:${result.construct_reverse}`;
}

export default function Design() {
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [params, setParams, clearParams] = useWorkspaceState<DesignParams>("design.params", DEFAULTS);
  const [result, setResult, clearResult] = useWorkspaceState<AssemblyResult | null>("design.result", null);
  const [assembledToken, setAssembledToken, clearAssembledToken] = useWorkspaceState(
    "design.assembledToken", "",
  );
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const [saved, setSaved, clearSaved] = useWorkspaceState("design.saved", "");
  const [exportOpen, setExportOpen] = useState(false);
  const [copied, setCopied] = useState(false);
  const [experience, setExperience] = useWorkspaceState<"guided" | "expert">(
    "design.experience", "guided",
  );
  const location = useLocation();
  const navigate = useNavigate();
  const [hybridDetail, setHybridDetail, clearHybridDetail] = useWorkspaceState<"simple" | "detailed">(
    "design.hybridDetail", "simple",
  );

  useEffect(() => {
    const handed = location.state as { sequence?: string; isCoding?: boolean } | null;
    const sequence = handed?.sequence;
    if (sequence) {
      setParams((current) => ({
        ...current,
        sequence,
        ...(typeof handed.isCoding === "boolean" ? { is_coding: handed.isCoding } : {}),
      }));
      setResult(null);
      setAssembledToken("");
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
      setAssembledToken("");
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
      setAssembledToken((current) => current === assemblyTokenFor(data) ? current : "");
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

  async function exportConstruct(
    filetype: "genbank" | "fasta" | "oligos" | "all-sequences" | "sbol3",
  ) {
    const safe = (params.name || "construct").replace(/\s+/g, "_");
    const names = { genbank: `${safe}.gb`, fasta: `${safe}.fasta`,
                    oligos: `${safe}_oligos.fasta`,
                    "all-sequences": `${safe}_all_sequences.fasta`,
                    sbol3: `${safe}.sbol.json` };
    try {
      await api.download(
        `/api/design/assembly/export/?filetype=${filetype}`, params, names[filetype],
      );
    } catch {
      setError("The download failed. Try designing again first.");
    }
  }

  const verified = result !== null && result.verification.length === 0;
  const requiresFragmentAssembly = (result?.fragment_count ?? 0) >= 2;
  const assemblyComplete = verified && (
    !requiresFragmentAssembly
    || (assembledToken !== "" && assembledToken === assemblyTokenFor(result))
  );
  const preflightPassed = verified && (result?.preflight?.can_export ?? true);
  const canExport = preflightPassed && assemblyComplete;

  function clearWorkspace() {
    clearParams();
    clearResult();
    clearSaved();
    clearAssembledToken();
    setError("");
    setExportOpen(false);
    setCopied(false);
    clearHybridDetail();
  }

  function inspectHybridization() {
    if (!result || !assemblyComplete) return;
    navigate("/hybridize", {
      state: {
        tool: "hybridization",
        first: result.construct_forward,
        second: result.construct_reverse,
        name: params.name || "designed construct",
        leftEnzyme: params.left_enzyme,
        rightEnzyme: params.right_enzyme,
        orfStart: result.ssd.orf_start,
        insertAnnotations: result.ssd.segments.map((segment) => ({
          name: segment.name.toLowerCase() === "insert" ? `${params.name || "construct"} target` : segment.name,
          type: "misc_feature",
          start: segment.start, end: segment.end, direction: 1,
          color: segmentColour(segment.name),
        })),
        autoRun: true,
      },
    });
  }

  async function copyConstruct() {
    if (!result) return;
    await navigator.clipboard.writeText(result.construct_forward);
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1800);
  }

  function simulateAssembly() {
    if (!result || !verified || !requiresFragmentAssembly) return;
    setAssembledToken(assemblyTokenFor(result));
  }

  const status = busy
    ? "Designing…"
    : result === null
      ? ""
      : verified && !requiresFragmentAssembly
        ? `Single-fragment design ready: one ${result.construct_length} bp duplex; no fragment assembly is required.`
        : assemblyComplete
        ? `ESD assembly complete: ${result.fragment_count} fragments reconstruct one ${result.construct_length} bp duplex.`
        : verified
          ? `ESD plan ready: ${result.fragment_count} fragments, ${result.oligo_count} oligos to order.`
        : "Design failed verification. Do not order these oligos.";

  return (
    <>
      <LiveStatus message={status} />

      <div className="topbar design-topbar">
        <div className="grow design-heading">
          <h1>Design a construct</h1>
          <CoreWorkflowTrail active="design" />
        </div>
        <button
          className="btn btn-primary design-save"
          onClick={() => design(true)}
          disabled={busy || !canExport}
        >
          <Icon name="check" size={18} />
          Save project
        </button>
        <div className="design-export">
          <button
            className="btn btn-outline"
            onClick={() => setExportOpen((open) => !open)}
            disabled={!canExport}
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
              <button role="menuitem" onClick={() => void exportConstruct("all-sequences")}>
                All sequences FASTA
              </button>
              <button role="menuitem" onClick={() => void download("order-sheet")}>Oligo CSV</button>
              <button role="menuitem" onClick={() => void exportConstruct("oligos")}>Oligo FASTA</button>
              <button role="menuitem" onClick={() => void download("protocol")}>Protocol</button>
              <button role="menuitem" onClick={() => void exportConstruct("genbank")}>GenBank</button>
              <button role="menuitem" onClick={() => void exportConstruct("sbol3")}>SBOL 3</button>
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
            <div className="card-head">
              <h2 style={{ flex: 1 }}>Insert</h2>
              <div className="mode-switch" role="group" aria-label="Design detail level">
                <button className={experience === "guided" ? "active" : ""} onClick={() => {
                  setExperience("guided");
                  setParams((current) => ({ ...current, cleavage_site: "Thrombin", include_his_tag: true, include_linkers: true, remove_stop: false, target_oligo_length: 90, overhang_length: 4 }));
                  setResult(null);
                  setAssembledToken("");
                }}>Guided</button>
                <button className={experience === "expert" ? "active" : ""} onClick={() => setExperience("expert")}>Expert</button>
              </div>
            </div>
            <div className="card-body">
              {experience === "guided" && (
                <div className="notice notice-info compact">
                  Validated defaults add a 6×His tag, flexible linkers and a Thrombin site, using 90 nt oligos with 4 nt assembly junctions.
                </div>
              )}
              <InsertForm params={params} catalogue={catalogue} onChange={set} expert={experience === "expert"} />
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
                      <strong>
                        {requiresFragmentAssembly ? "ESD plan verified" : "Single-fragment design verified"}
                      </strong>
                      <span>
                        {requiresFragmentAssembly
                          ? "The fragment set is internally consistent and ready for assembly simulation."
                          : "No inter-fragment assembly is required; proceed to duplex hybridization."}
                      </span>
                    </div>
                  ) : (
                    <div className="design-verdict-copy">
                      <strong>Do not order</strong>
                      <span>{result.verification.join(" ")}</span>
                    </div>
                  )}
                  {verified && (
                    <div className="design-verdict-date">
                      <span>Design checked on</span>
                      <strong>{new Intl.DateTimeFormat("en-GB", {
                        day: "2-digit", month: "short", year: "numeric",
                      }).format(new Date())}</strong>
                    </div>
                  )}
                </div>

                <PreflightPanel report={result.preflight} />

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

                {requiresFragmentAssembly && (
                  <div className={`card esd-assembly-card ${assemblyComplete ? "complete" : ""}`}>
                    <div className="card-head">
                      <div className="grow">
                        <h2>ESD fragment assembly</h2>
                        <span className="label">1. Assemble the designed fragments</span>
                      </div>
                      <span className={`pill ${assemblyComplete ? "pill-ok" : "pill-warn"}`}>
                        {assemblyComplete ? "Assembled" : "Pending"}
                      </span>
                    </div>
                    <div className="card-body">
                      <div className="esd-assembly-scroll" aria-label="ESD fragment assembly order">
                        <div className="esd-fragment-chain">
                          {result.fragments.map((fragment, index) => (
                            <div className="esd-fragment-step" key={fragment.index}>
                              <div className="esd-fragment-node">
                                <strong>{fragment.name}</strong>
                                <small>{fragment.forward_length} / {fragment.reverse_length} nt</small>
                              </div>
                              {index < result.fragments.length - 1 && (
                                <div className="esd-junction-node">
                                  <span>{result.junction_overhangs[index]}</span>
                                  <Icon name="arrowRight" size={16} />
                                </div>
                              )}
                            </div>
                          ))}
                          <div className={`esd-product-node ${assemblyComplete ? "complete" : ""}`}>
                            <Icon name={assemblyComplete ? "check" : "helix"} size={20} />
                            <div>
                              <strong>Assembled duplex</strong>
                              <small>{result.construct_length} bp</small>
                            </div>
                          </div>
                        </div>
                      </div>
                      <div className="esd-assembly-action">
                        <div>
                          <strong>{assemblyComplete ? "Exact reconstruction confirmed" : "Assembly simulation required"}</strong>
                          <span>
                            {assemblyComplete
                              ? "The assembled forward and reverse strands are now fixed for hybridization."
                              : "G-Synth will ligate the fragments in order and verify both reconstructed strands."}
                          </span>
                        </div>
                        <button
                          className="btn btn-primary"
                          type="button"
                          onClick={simulateAssembly}
                          disabled={!preflightPassed || assemblyComplete}
                        >
                          <Icon name={assemblyComplete ? "check" : "plate"} size={17} />
                          {assemblyComplete ? "Assembly complete" : "Assemble fragments"}
                        </button>
                      </div>
                    </div>
                  </div>
                )}

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
                    <div style={{ flex: 1 }}>
                      <h2>Hybridization</h2>
                      <span className="label">
                        {requiresFragmentAssembly
                          ? "2. Verify the assembled duplex before cloning"
                          : "Verify the designed duplex before cloning"}
                      </span>
                    </div>
                    <div className="seg-toggle" role="group" aria-label="Hybridization detail level">
                      {(["simple", "detailed"] as const).map((level) => (
                        <button
                          key={level}
                          type="button"
                          className={hybridDetail === level ? "on" : ""}
                          aria-pressed={hybridDetail === level}
                          onClick={() => setHybridDetail(level)}
                        >
                          {level[0].toUpperCase() + level.slice(1)}
                        </button>
                      ))}
                    </div>
                  </div>
                  <div className="card-body">
                    {requiresFragmentAssembly && !assemblyComplete && (
                      <div className="notice notice-info compact">
                        Complete the ESD fragment assembly above to fix the duplex used here and in cloning.
                      </div>
                    )}
                    <div className="design-hybrid-enzyme-row">
                      <EnzymePicker
                        id="hybrid-left-enzyme"
                        label="Left cloning enzyme"
                        enzymes={catalogue?.enzymes ?? []}
                        value={params.left_enzyme}
                        onChange={(value) => set("left_enzyme", value)}
                        disabled={!catalogue}
                      />
                      <EnzymePicker
                        id="hybrid-right-enzyme"
                        label="Right cloning enzyme"
                        enzymes={catalogue?.enzymes ?? []}
                        value={params.right_enzyme}
                        onChange={(value) => set("right_enzyme", value)}
                        disabled={!catalogue}
                      />
                      <p className="field-hint">
                        Inherited from this design. Changing either enzyme invalidates the
                        current molecules and returns you to Update design; G-Synth never
                        relabels an existing sticky end as a different enzyme.
                      </p>
                    </div>

                    {hybridDetail === "simple" ? (
                      <div className="design-hybrid-simple">
                        <DesignEnd end={result.terminal_ends[0]} />
                        <div className="design-hybrid-core">
                          <span className="label">Antiparallel duplex</span>
                          <strong>{result.construct_length} paired columns</strong>
                          <div className="design-hybrid-core-bar" aria-hidden="true" />
                          <small>
                            {result.duplex.mismatches.length === 0
                              ? "Both reconstructed strands are complementary across the intended duplex."
                              : `${result.duplex.mismatches.length} mismatch positions require review.`}
                          </small>
                        </div>
                        <DesignEnd end={result.terminal_ends[1]} />
                      </div>
                    ) : (
                      <DuplexView duplex={result.duplex} />
                    )}

                    <div className="design-hybrid-actions">
                      <span className="design-hybrid-next-copy">
                        Hybridization is the required verification gate before cloning.
                      </span>
                      <button className="btn btn-primary" onClick={inspectHybridization} disabled={!canExport}>
                        Verify hybridization <Icon name="arrowRight" size={17} />
                      </button>
                    </div>
                  </div>
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Oligos to order</h2>
                    <button
                      className="btn btn-outline btn-small"
                      type="button"
                      onClick={() => void exportConstruct("all-sequences")}
                      disabled={!canExport}
                    >
                      <Icon name="arrowRight" size={16} />
                      Export all sequences
                    </button>
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

function DesignEnd({ end }: { end: TerminalEnd | undefined }) {
  if (!end || !end.overhang) {
    return (
      <div className="hybrid-end-card blunt">
        <span className="hybrid-end-label">{end?.side ?? "terminal"} end</span>
        <strong>Blunt</strong>
        <small>{end?.enzyme ?? "No protruding strand"}</small>
      </div>
    );
  }
  return (
    <div className="hybrid-end-card sticky">
      <span className="hybrid-end-label">{end.side} cohesive end</span>
      <div className="hybrid-end-title"><strong>{end.kind} overhang</strong><span>{end.overhang.length} nt</span></div>
      <code>5′-{end.overhang}-3′</code>
      <small>Generated by {end.enzyme}</small>
    </div>
  );
}
