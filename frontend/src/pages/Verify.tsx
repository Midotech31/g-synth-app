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

/**
 * The end of the workflow: you built it, now check it is what you designed.
 *
 * All three tools on this page operate on one saved construct, because that
 * is how they are used — the ligation ratios, the primers to order, and the
 * reads that come back are all about the same molecule. Making the user
 * paste it three times would be the wrong shape.
 */

type Tab = "reads" | "primers" | "ligation";

export default function Verify() {
  const [projects, setProjects] = useState<ProjectSummary[]>([]);
  const [project, setProject] = useState<Project | null>(null);
  const [tab, setTab] = useState<Tab>("reads");
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);

  const [reads, setReads] = useState("");
  const [report, setReport] = useState<VerifyReport | null>(null);
  const [primers, setPrimers] = useState<PrimerSet | null>(null);
  const [ligation, setLigation] = useState<LigationReaction[] | null>(null);
  const [vectorNg, setVectorNg] = useState(50);
  const [trim, setTrim] = useState(30);

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

  /** The insert's span inside the construct, when the payload knows it. */
  const data = (project?.data ?? {}) as Record<string, number | undefined>;
  const insertStart = data.insert_start;
  const insertEnd = data.insert_end;
  const hasRegion = insertStart !== undefined && insertEnd !== undefined;
  const circular = (project?.data as { topology?: string })?.topology === "circular";

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
    if (!Object.keys(parsed).length) {
      setError("Paste at least one read — FASTA, or just the bases.");
      return;
    }
    setBusy(true);
    setError("");
    try {
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
    if (!project || !hasRegion) return;
    const insertLength = insertEnd! - insertStart!;
    setBusy(true);
    setError("");
    try {
      const result = await api.ligation({
        vector_length: (data.backbone_length as number) || project.sequence.length,
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

  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>Check the clone</h1>
          <p className="sub">
            Ligation amounts, sequencing primers, and what the reads say when
            they come back — all against one saved construct.
          </p>
        </div>
      </div>

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
        {error && <div className="notice notice-error">{error}</div>}

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
                </div>
              )}

              {tab === "reads" && project && (
                <div className="field">
                  <label htmlFor="reads">Sequencing reads</label>
                  <textarea
                    id="reads"
                    value={reads}
                    onChange={(e) => setReads(e.target.value)}
                    rows={10}
                    className="mono"
                    style={{ fontSize: "0.74rem" }}
                    placeholder={">T7-F\nGATCC...\n>T7-R\nCTAGG..."}
                  />
                  <span className="label">
                    FASTA, or just the bases for a single read
                  </span>
                  <label htmlFor="trim">Ignore low-quality bases at each end</label>
                  <input
                    id="trim"
                    type="number"
                    min={0}
                    max={200}
                    value={trim}
                    onChange={(e) => setTrim(Number(e.target.value))}
                  />
                  <span className="label">
                    {trim} bases from the start and {trim} from the end · use 0 for a cleaned sequence
                  </span>
                </div>
              )}

              {tab === "ligation" && project && (
                <div className="field">
                  <label htmlFor="vng">Vector in the reaction (ng)</label>
                  <input id="vng" type="number" min={1} max={1000} value={vectorNg}
                         onChange={(e) => setVectorNg(Number(e.target.value))} />
                </div>
              )}

              <button
                className="btn btn-primary"
                disabled={busy || !project || (tab !== "reads" && !hasRegion)}
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
            <div className="seg-toggle" style={{ alignSelf: "flex-start" }}>
              {(["reads", "primers", "ligation"] as Tab[]).map((option) => (
                <button
                  key={option}
                  type="button"
                  className={tab === option ? "on" : ""}
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
                  <div className={`notice ${report.is_verified ? "notice-ok" : "notice-error"}`}>
                    {report.is_verified ? (
                      <>
                        <strong>It is the design.</strong> Every read agrees
                        with the construct over {report.coverage}% of the
                        region checked.
                      </>
                    ) : (
                      <>
                        <strong>{report.differences.length} difference
                        {report.differences.length === 1 ? "" : "s"}.</strong>{" "}
                        The clone is not what was designed.
                      </>
                    )}
                  </div>

                  {report.differences.length > 0 && (
                    <div className="card">
                      <div className="card-head"><h2>What differs</h2></div>
                      <div className="card-body">
                        <ul className="difference-list">
                          {report.differences.map((d) => (
                            <li key={`${d.kind}-${d.position}-${d.found}`}
                                className={d.silent ? "silent" : ""}>
                              {d.description}
                            </li>
                          ))}
                        </ul>
                        <p className="note" style={{ marginTop: "0.6rem" }}>
                          Positions are in the construct, counting from 1.
                          Silent changes leave the protein alone.
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
                    <div className="glyph">🔬</div>
                    <strong>No reads compared yet</strong>
                    <span>
                      Pick a construct and paste what the sequencing facility
                      sent back. Either orientation is fine.
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
                                {p.direction === 1 ? "→" : "←"} {p.reads_from + 1}–{p.reads_to}
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
                    <div className="glyph">🎯</div>
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
                    <div className="glyph">⚖️</div>
                    <strong>No amounts yet</strong>
                    <span>Pick a construct to work out how much insert to add.</span>
                  </div>
                </div>
              )
            )}
          </div>
        </div>
      </div>
    </>
  );
}
