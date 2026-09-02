import { useCallback, useEffect, useState } from "react";

import {
  ApiError,
  api,
  type LigationReaction,
  type PrimerSet,
  type Project,
  type ProjectSummary,
  type VerifyReport,
} from "../api/client";
import Icon from "../components/Icon";
import CoverageMap from "../components/CoverageMap";
import LiveStatus from "../components/LiveStatus";
import PreflightPanel from "../components/PreflightPanel";
import ReferenceAlignment from "../components/ReferenceAlignment";
import TraceView from "../components/TraceView";
import { useWorkspaceState } from "../state/WorkspaceStateContext";

/**
 * The end of the workflow: you built it, now check it is what you designed.
 *
 * All three tools on this page operate on one saved construct, because that
 * is how they are used — the ligation ratios, the primers to order, and the
 * reads that come back are all about the same molecule. Making the user
 * paste it three times would be the wrong shape.
 */

type Tab = "reads" | "primers" | "ligation";

export type ProjectCapabilities = {
  insertStart?: number;
  insertEnd?: number;
  backboneLength?: number;
  hasRegion: boolean;
  hasLigationContext: boolean;
  circular: boolean;
};

/** Derive workflow capabilities from explicit metadata, never sequence length. */
export function projectCapabilities(project: Project | null): ProjectCapabilities {
  const insertStart = project?.data.insert_start;
  const insertEnd = project?.data.insert_end;
  const backboneLength = project?.data.backbone_length;
  const hasRegion = Number.isInteger(insertStart)
    && Number.isInteger(insertEnd)
    && insertStart! >= 0
    && insertEnd! > insertStart!
    && insertEnd! <= (project?.sequence.length ?? 0);
  const hasLigationContext = hasRegion
    && typeof backboneLength === "number"
    && Number.isFinite(backboneLength)
    && backboneLength > 0;

  return {
    insertStart,
    insertEnd,
    backboneLength,
    hasRegion,
    hasLigationContext,
    circular: project?.data.topology === "circular",
  };
}

export default function Verify() {
  const [projects, setProjects] = useState<ProjectSummary[]>([]);
  const [project, setProject, clearProject] = useWorkspaceState<Project | null>("verify.project", null);
  const [tab, setTab, clearTab] = useWorkspaceState<Tab>("verify.tab", "reads");
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);

  const [reads, setReads, clearReads] = useWorkspaceState("verify.reads", "");
  const [traceFiles, setTraceFiles, clearTraceFiles] = useWorkspaceState<File[]>("verify.traceFiles", []);
  const [report, setReport, clearReport] = useWorkspaceState<VerifyReport | null>("verify.report", null);
  const [primers, setPrimers, clearPrimers] = useWorkspaceState<PrimerSet | null>("verify.primers", null);
  const [ligation, setLigation, clearLigation] = useWorkspaceState<LigationReaction[] | null>("verify.ligation", null);
  const [vectorNg, setVectorNg, clearVectorNg] = useWorkspaceState("verify.vectorNg", 50);
  const [trim, setTrim, clearTrim] = useWorkspaceState("verify.trim", 30);

  useEffect(() => {
    api.listProjects()
      .then((page) => setProjects(page.results))
      .catch(() => setError("Could not load your projects."));
  }, []);

  const open = useCallback(async (id: number) => {
    setError("");
    setReport(null);
    setPrimers(null);
    setLigation(null);
    try {
      setProject(await api.getProject(id));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Could not open that project.");
    }
  }, []);

  /** The insert span enables primers; ligation additionally needs a backbone. */
  const {
    insertStart,
    insertEnd,
    backboneLength,
    hasRegion,
    hasLigationContext,
    circular,
  } = projectCapabilities(project);

  /** One FASTA-ish blob in, named reads out. Bare sequence is one read. */
  function parseReads(text: string): Record<string, string> {
    const out: Record<string, string> = {};
    const blocks = text.split(/^>/m).filter((b) => b.trim());
    if (blocks.length === 0) return out;
    if (!text.trimStart().startsWith(">")) {
      return { read: text.replace(/[^A-Za-z]/g, "") };
    }
    for (const block of blocks) {
      const [header, ...rest] = block.split(/\r?\n/);
      const name = header.trim().split(/\s+/)[0] || `read ${Object.keys(out).length + 1}`;
      out[name] = rest.join("").replace(/[^A-Za-z]/g, "");
    }
    return out;
  }

  async function runVerify() {
    if (!project) return;
    const parsed = parseReads(reads);
    if (!traceFiles.length && !Object.keys(parsed).length) {
      setError("Add an ABIF or SCF trace, or paste the bases.");
      return;
    }
    setBusy(true);
    setError("");
    try {
      // A trace carries the confidence of every base; letters do not. When
      // both are given the traces win — there is no reason to discard the
      // one piece of evidence that separates a mutation from a bad call.
      const common = {
        design: project.sequence,
        circular,
        region_start: hasRegion ? insertStart : null,
        region_end: hasRegion ? insertEnd : null,
        coding_start: hasRegion ? insertStart : null,
        coding_end: hasRegion ? insertEnd : null,
      };
      if (traceFiles.length) {
        setReport(await api.verifyTraces({ ...common, files: traceFiles }));
        return;
      }
      setReport(await api.verify({
        design: project.sequence,
        reads: parsed,
        circular,
        trim,
        region_start: hasRegion ? insertStart : null,
        region_end: hasRegion ? insertEnd : null,
        coding_start: hasRegion ? insertStart : null,
        coding_end: hasRegion ? insertEnd : null,
      }));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "The comparison failed.");
    } finally {
      setBusy(false);
    }
  }

  async function runPrimers() {
    if (!project || !hasRegion) return;
    setBusy(true);
    setError("");
    try {
      setPrimers(await api.primers({
        template: project.sequence,
        target_start: insertStart!,
        target_end: insertEnd!,
        circular,
        name: project.name.replace(/\s+/g, "_").slice(0, 20) || "seq",
      }));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Primer design failed.");
    } finally {
      setBusy(false);
    }
  }

  /** A primer set is ordered, not read on screen. */
  async function exportPrimers(filetype: "csv" | "fasta") {
    if (!project || !hasRegion) return;
    const safe = project.name.replace(/\s+/g, "_").slice(0, 20) || "seq";
    try {
      await api.download(
        `/api/design/primers/export/?filetype=${filetype}`,
        {
          template: project.sequence,
          target_start: insertStart!,
          target_end: insertEnd!,
          circular,
          name: safe,
        } as never,
        `${safe}_primers.${filetype}`,
      );
    } catch {
      setError("The download failed. Design the primers again first.");
    }
  }

  async function runLigation() {
    if (!project || !hasLigationContext) return;
    const insertLength = insertEnd! - insertStart!;
    setBusy(true);
    setError("");
    try {
      const result = await api.ligation({
        vector_length: backboneLength!,
        insert_length: insertLength,
        vector_ng: vectorNg,
        ratios: [1, 3, 5],
      });
      setLigation(result.reactions);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "The calculation failed.");
    } finally {
      setBusy(false);
    }
  }

  /** Each tab answers its own question, so each has its own verdict. */
  function verdict(): string {
    if (tab === "reads") {
      if (!report) return "";
      switch (report.verification_state) {
        case "fully_verified": return `Fully verified: the complete requested region agrees with the design.`;
        case "differences_detected": return `${report.differences.length} sequence difference${report.differences.length === 1 ? "" : "s"} detected.`;
        case "reads_unplaced": return "Reads supplied, but none could be placed on this design.";
        case "partial_match": return `Partial match: covered bases agree, but only ${report.coverage}% of the region was read.`;
        default: return "Verification has not produced a conclusive result.";
      }
    }
    if (tab === "primers") {
      if (!primers) return "";
      return primers.covers_target
        ? `${primers.primers.length} primers, together reading the whole insert on both strands.`
        : "The primers do not cover the whole insert.";
    }
    if (!ligation) return "";
    return `Amounts worked out for ${ligation.length} ligation reactions.`;
  }

  const status = busy
    ? tab === "reads"
      ? "Comparing the reads to the design…"
      : tab === "primers"
        ? "Designing primers…"
        : "Working out the amounts…"
    : verdict();

  function clearWorkspace() {
    clearProject();
    clearTab();
    clearReads();
    clearTraceFiles();
    clearReport();
    clearPrimers();
    clearLigation();
    clearVectorNg();
    clearTrim();
    setError("");
  }

  return (
    <>
      <LiveStatus message={status} />

      <div className="topbar">
        <div className="grow">
          <h1>Validate a construct</h1>
          <p className="sub">
            Keep ligation calculations, sequencing primers and returned reads
            anchored to the same saved reference.
          </p>
        </div>
        <button className="btn btn-outline" onClick={clearWorkspace} disabled={busy}>
          Clear
        </button>
      </div>

      <div
        className="content"
        style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}
        aria-busy={busy}
      >
        {error && <div className="notice notice-error" role="alert" id="verify-error">{error}</div>}

        <div className="design-layout">
          {/* ── Pick the construct ─────────────────────────────────────── */}
          <div className="card">
            <div className="card-head"><h2>Construct</h2></div>
            <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.85rem" }}>
              <div className="field">
                <label htmlFor="proj">Saved project</label>
                <select
                  id="proj"
                  value={project?.id ?? ""}
                  onChange={(e) => e.target.value && void open(Number(e.target.value))}
                >
                  <option value="">Choose a construct…</option>
                  {projects.map((p) => (
                    <option key={p.id} value={p.id}>
                      {p.name} · {p.module.replace(/_/g, " ")}
                    </option>
                  ))}
                </select>
                {projects.length === 0 && (
                  <span className="label">
                    Save a design or a plasmid first
                  </span>
                )}
              </div>

              {project && (
                <div className="vector-brief">
                  <div className="vector-facts">
                    <span><b>{project.sequence.length.toLocaleString()}</b> bp</span>
                    <span>{circular ? "circular" : "linear"}</span>
                    {hasRegion && (
                      <span>insert <b>{insertEnd! - insertStart!}</b> bp</span>
                    )}
                  </div>
                  {!hasRegion && (
                    <p className="note vector-note">
                      This project does not record where the insert sits, so
                      primers and ligation amounts cannot be worked out from
                      it. Reads can still be compared against the whole
                      sequence.
                    </p>
                  )}
                  {hasRegion && !hasLigationContext && (
                    <p className="note vector-note">
                      Reads and sequencing primers are available for this
                      assembly. Ligation amounts require a saved cloned
                      plasmid, because the vector backbone length is part of
                      the molar calculation.
                    </p>
                  )}
                </div>
              )}

              {tab === "reads" && project && (
                <div className="field verify-evidence-field">
                  <span className="field-label">Sanger trace evidence</span>
                  <label className="trace-upload" htmlFor="traces">
                    <span className="trace-upload-icon" aria-hidden="true">
                      <Icon name="microscope" size={22} />
                    </span>
                    <span className="trace-upload-copy">
                      <strong>
                        {traceFiles.length
                          ? `${traceFiles.length} trace${traceFiles.length === 1 ? "" : "s"} selected`
                          : "Choose AB1 or SCF files"}
                      </strong>
                      <small>
                        {traceFiles.length
                          ? traceFiles.map((file) => file.name).join(" · ")
                          : "Forward and reverse traces can be uploaded together"}
                      </small>
                    </span>
                    <span className="btn btn-outline trace-upload-action" aria-hidden="true">
                      Browse
                    </span>
                  </label>
                  <input
                    id="traces"
                    className="sr-only"
                    type="file"
                    accept=".ab1,.scf,application/octet-stream"
                    multiple
                    // "Add an ABIF or SCF trace, or paste the bases" is a complaint
                    // about these two fields; it is attached to them so it is
                    // read when either is reached, not only when it appears.
                    aria-describedby={error ? "traces-hint verify-error" : "traces-hint"}
                    onChange={(e) => {
                      setTraceFiles(Array.from(e.target.files ?? []));
                      setReport(null);
                    }}
                  />
                  <span className="field-hint" id="traces-hint">
                    {traceFiles.length
                      ? "Trace quality, base calls and chromatogram peaks will be evaluated together."
                      : "Preferred: trace files retain the quality evidence needed to distinguish a mutation from a weak peak."}
                  </span>
                </div>
              )}

              {tab === "reads" && project && (
                <div className="field verify-text-reads">
                  <div className="verify-or" aria-hidden="true"><span>or</span></div>
                  <label htmlFor="reads">
                    {traceFiles.length ? "Use text reads instead" : "Paste text reads"}
                  </label>
                  <textarea
                    id="reads"
                    value={reads}
                    onChange={(e) => setReads(e.target.value)}
                    rows={7}
                    className="mono"
                    style={{ fontSize: "0.74rem" }}
                    placeholder={">T7-F\nGATCC...\n>T7-R\nCTAGG..."}
                    aria-describedby={error ? "reads-hint verify-error" : "reads-hint"}
                  />
                  <span className="field-hint" id="reads-hint">
                    FASTA or plain bases. Text reads do not contain chromatogram quality evidence.
                  </span>
                  {!traceFiles.length && (
                    <details className="advanced-control">
                      <summary>Text-read trimming</summary>
                      <div className="field advanced-control-body">
                        <label htmlFor="trim">Ignore bases at each read end</label>
                        <input
                          id="trim"
                          type="number"
                          min={0}
                          max={200}
                          value={trim}
                          onChange={(e) => setTrim(Number(e.target.value))}
                          aria-describedby="trim-hint"
                        />
                        <span className="field-hint" id="trim-hint">
                          {trim} bases at the start and {trim} at the end. Use 0 for pre-cleaned reads.
                        </span>
                      </div>
                    </details>
                  )}
                </div>
              )}

              {tab === "ligation" && project && hasLigationContext && (
                <div className="field">
                  <label htmlFor="vng">Vector in the reaction (ng)</label>
                  <input id="vng" type="number" min={1} max={1000} value={vectorNg}
                         onChange={(e) => setVectorNg(Number(e.target.value))} />
                </div>
              )}

              <button
                className="btn btn-primary"
                disabled={busy || !project
                  || (tab === "primers" && !hasRegion)
                  || (tab === "ligation" && !hasLigationContext)}
                onClick={() => {
                  if (tab === "reads") void runVerify();
                  else if (tab === "primers") void runPrimers();
                  else void runLigation();
                }}
              >
                {busy && <span className="spinner" />}
                {tab === "reads" ? "Compare to the design"
                  : tab === "primers" ? "Design primers"
                  : "Work out the amounts"}
              </button>
            </div>
          </div>

          {/* ── Results ────────────────────────────────────────────────── */}
          <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
            <div
              className="seg-toggle"
              style={{ alignSelf: "flex-start" }}
              role="group"
              aria-label="What to check"
            >
              {(["reads", "primers", "ligation"] as Tab[]).map((option) => (
                <button
                  key={option}
                  type="button"
                  className={tab === option ? "on" : ""}
                  aria-pressed={tab === option}
                  onClick={() => setTab(option)}
                >
                  {option === "reads" ? "Sequencing reads"
                    : option === "primers" ? "Primers"
                    : "Ligation"}
                </button>
              ))}
            </div>

            {tab === "reads" && (
              report ? (
                <>
                  <div className={`notice ${report.verification_state === "fully_verified" ? "notice-ok" : report.verification_state === "differences_detected" ? "notice-error" : "notice-info"}`}>
                    {report.verification_state === "fully_verified" ? (
                      <>
                        <strong>Fully verified.</strong> The complete requested
                        region is covered and agrees with the design.
                      </>
                    ) : report.verification_state === "differences_detected" ? (
                      <>
                        <strong>{report.differences.length} difference
                        {report.differences.length === 1 ? "" : "s"}.</strong>{" "}
                        The covered clone sequence is not what was designed.
                      </>
                    ) : report.verification_state === "reads_unplaced" ? (
                      <><strong>Reads could not be placed.</strong> Check that these reads belong to this construct, then review orientation and trimming.</>
                    ) : report.verification_state === "not_checked" ? (
                      <><strong>Not checked.</strong> No usable sequencing evidence was available.</>
                    ) : (
                      <>
                        <strong>Partial match.</strong> The reads agree
                        wherever they align, but cover only {report.coverage}%
                        of the requested region. Sequence the remaining gap
                        {report.gaps.length === 1 ? "" : "s"} before accepting
                        this clone.
                      </>
                    )}
                  </div>

                  <PreflightPanel report={report.preflight} />

                  {report.raw_consensus && report.consensus && (
                    <div className="card">
                      <div className="card-head"><h2>Forward/reverse consensus</h2></div>
                      <div className="card-body">
                        <div className="vector-facts" aria-label="Consensus sequencing metrics">
                          <span><b>{report.raw_consensus.coverage}%</b> raw consensus coverage</span>
                          <span><b>{report.raw_consensus.identity}%</b> raw consensus identity</span>
                          <span><b>{report.raw_consensus.bidirectional_overlap}%</b> covered by both strands</span>
                          <span><b>{report.raw_consensus.bidirectional_agreement}%</b> F/R overlap agreement</span>
                          <span><b>{report.consensus.coverage}%</b> Q{report.quality_cutoff ?? 13} consensus coverage</span>
                        </div>
                        <p className="note" style={{ marginTop: "0.7rem" }}>
                          Forward and reverse calls are oriented and merged before coverage is calculated.
                          Raw 100% coverage means every reference position has a consensus call; only the
                          bidirectional percentage is supported by both strands. The quality-gated value
                          determines whether the clone can be accepted as fully verified.
                        </p>
                      </div>
                    </div>
                  )}

                  {report.differences.length > 0 && (
                    <div className="card">
                      <div className="card-head"><h2>What differs</h2></div>
                      <div className="card-body">
                        <ul className="difference-list">
                          {report.differences.map((d) => {
                            const peaks = report.trace_windows?.find(
                              (w) => w.position === d.position,
                            );
                            return (
                              <li key={`${d.kind}-${d.position}-${d.found}`}
                                  className={
                                    (d.silent ? "silent " : "") +
                                    (d.confident === false ? "unconfident" : "")
                                  }>
                                {d.description}
                                {peaks && <TraceView window={peaks} />}
                              </li>
                            );
                          })}
                        </ul>
                        <p className="note" style={{ marginTop: "0.6rem" }}>
                          Positions are in the construct, counting from 1.
                          Silent changes leave the protein alone.
                          {report.differences.some((d) => d.confident === false) && (
                            <>
                              {" "}Differences marked Q&lt;20 sit on a peak the
                              basecaller was not sure of — read those again
                              before acting on them.
                            </>
                          )}
                        </p>
                      </div>
                    </div>
                  )}

                  <div className="card">
                    <div className="card-head">
                      <h2 style={{ flex: 1 }}>Reads</h2>
                      <span className="label">
                        {report.coverage}% covered
                        {report.fully_covered ? "" : ` · ${report.gaps.length} gap(s)`}
                      </span>
                    </div>
                    <div className="table-scroll">
                      <div className="card-body" style={{ paddingBottom: 0 }}>
                        <CoverageMap report={report} />
                      </div>
                      <table className="data">
                        <thead>
                          <tr>
                            <th>Read</th><th>Length</th><th>Aligned to</th>
                            <th>Strand</th><th>Identity</th><th>Differences</th>
                          </tr>
                        </thead>
                        <tbody>
                          {report.reads.map((r) => (
                            <tr key={r.name}>
                              <td className="mono">{r.name}</td>
                              <td className="num">{r.length}</td>
                              <td className="num">{r.start + 1}–{r.end}</td>
                              <td>{r.reverse_complemented ? "reverse" : "forward"}</td>
                              <td className="num">{r.identity}%</td>
                              <td className="num">{r.difference_count}</td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </div>

                  {report.warnings.length > 0 && (
                    <div className="card">
                      <div className="card-head"><h2>Notes</h2></div>
                      <div className="card-body">
                        <ul style={{ margin: 0, paddingLeft: "1.1rem", color: "var(--ink-soft)" }}>
                          {report.warnings.map((w) => (
                            <li key={w} style={{ marginBottom: "0.3rem", lineHeight: 1.5 }}>{w}</li>
                          ))}
                        </ul>
                      </div>
                    </div>
                  )}

                  {report.reads.some((read) => read.warnings.length > 0) && (
                    <div className="card">
                      <div className="card-head"><h2>Read quality notes</h2></div>
                      <div className="card-body">
                        <ul style={{ margin: 0, paddingLeft: "1.1rem", color: "var(--ink-soft)" }}>
                          {report.reads.flatMap((read) =>
                            read.warnings.map((warning) => (
                              <li key={`${read.name}-${warning}`} style={{ marginBottom: "0.3rem", lineHeight: 1.5 }}>
                                <strong>{read.name}:</strong> {warning}
                              </li>
                            )),
                          )}
                        </ul>
                      </div>
                    </div>
                  )}
                </>
              ) : (
                <div className="card">
                  <div className="empty">
                    <Icon name="microscope" size={38} className="glyph" />
                    <strong>No reads compared yet</strong>
                    <span>
                      {project
                        ? "Add trace files or paste text reads, then compare them with this reference. Either orientation is fine."
                        : "Choose the construct that these sequencing reads are meant to validate."}
                    </span>
                  </div>
                </div>
              )
            )}

            {tab === "primers" && (
              primers ? (
                <>
                  <div className={`notice ${primers.covers_target ? "notice-ok" : "notice-error"}`}>
                    {primers.covers_target ? (
                      <>
                        <strong>{primers.primers.length} primers.</strong> Together
                        they read the whole insert, on both strands.
                      </>
                    ) : (
                      <><strong>Incomplete.</strong> {primers.warnings.join(" ")}</>
                    )}
                  </div>
                  <div className="card">
                    <div className="card-head">
                      <h2 style={{ flex: 1 }}>Primers to order</h2>
                      <button className="btn btn-outline"
                              onClick={() => void exportPrimers("csv")}>
                        CSV
                      </button>
                      <button className="btn btn-outline"
                              onClick={() => void exportPrimers("fasta")}
                              title="For suppliers that take a FASTA upload">
                        FASTA
                      </button>
                    </div>
                    <div className="table-scroll">
                      <table className="data">
                        <thead>
                          <tr>
                            <th>Name</th><th>Sequence (5'→3')</th><th>Length</th>
                            <th>Tm</th><th>GC</th><th>Reads</th>
                          </tr>
                        </thead>
                        <tbody>
                          {primers.primers.map((p) => (
                            <tr key={p.name}>
                              <td className="mono">{p.name}</td>
                              <td className="mono seq-cell">{p.sequence}</td>
                              <td className="num">{p.length}</td>
                              <td className="num">{p.tm}</td>
                              <td className="num">{p.gc}</td>
                              <td className="num">
                                <Icon name={p.direction === 1 ? "arrowRight" : "arrowLeft"} size={14}
                                      title={p.direction === 1 ? "forward" : "reverse"} />{" "}
                                {p.reads_from + 1}–{p.reads_to}
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                    <div className="card-body" style={{ paddingTop: 0 }}>
                      <p className="note" style={{ margin: 0 }}>
                        Each primer sits back from what it reads: the first
                        fifty bases after a sequencing primer are noise.
                      </p>
                    </div>
                  </div>
                </>
              ) : (
                <div className="card">
                  <div className="empty">
                    <Icon name="target" size={38} className="glyph" />
                    <strong>No primers yet</strong>
                    <span>Pick a construct, then design primers that read its insert.</span>
                  </div>
                </div>
              )
            )}

            {tab === "ligation" && (
              ligation ? (
                <>
                  <div className="notice notice-info">
                    Set up all three. Nobody runs one ligation — they run a
                    small series and pick whichever plate gives colonies.
                  </div>
                  <div className="card">
                    <div className="card-head"><h2>Insert : vector</h2></div>
                    <div className="table-scroll">
                      <table className="data">
                        <thead>
                          <tr>
                            <th>Ratio</th><th>Vector</th><th>Insert</th>
                            <th>Vector</th><th>Insert</th>
                          </tr>
                        </thead>
                        <tbody>
                          {ligation.map((r) => (
                            <tr key={r.ratio}>
                              <td className="mono">{r.ratio}:1</td>
                              <td className="num">{r.vector_ng} ng</td>
                              <td className="num">{r.insert_ng} ng</td>
                              <td className="num">{r.vector_fmol} fmol</td>
                              <td className="num">{r.insert_fmol} fmol</td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                    <div className="card-body" style={{ paddingTop: 0 }}>
                      <p className="note" style={{ margin: 0 }}>
                        The ratio is molar, the amounts are mass. At equal
                        mass a 5.4 kb vector outnumbers a 150 bp insert
                        thirty-six to one — which is why that plate is empty.
                      </p>
                      {ligation[0]?.warnings.map((w) => (
                        <p key={w} className="note vector-note" style={{ marginTop: "0.5rem" }}>{w}</p>
                      ))}
                    </div>
                  </div>
                </>
              ) : (
                <div className="card">
                  <div className="empty">
                    <Icon name="scales" size={38} className="glyph" />
                    <strong>No amounts yet</strong>
                    <span>Pick a construct to work out how much insert to add.</span>
                  </div>
                </div>
              )
            )}
          </div>
        </div>

        {tab === "reads" && report && project && (report.trace_tracks?.length ?? 0) > 0 && (
          <div className="card alignment-card">
            <div className="card-body">
              <ReferenceAlignment reference={project.sequence} report={report} />
            </div>
          </div>
        )}
      </div>
    </>
  );
}
