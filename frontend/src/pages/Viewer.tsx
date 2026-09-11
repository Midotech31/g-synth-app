import { useEffect, useMemo, useRef, useState } from "react";
import { Link, useParams } from "react-router-dom";
import { SeqViz } from "seqviz";

import {
  ApiError,
  api,
  type Annotation,
  type DetectedFeature,
  type PreflightReport,
  type Project,
  type Provenance,
} from "../api/client";
import Icon from "../components/Icon";
import AnnotatedSequenceView from "../components/AnnotatedSequenceView";
import ExpandablePanel from "../components/ExpandablePanel";
import { featureLabel } from "../components/featureLabel";
import FeatureEvidence from "../components/FeatureEvidence";
import AnnotationEditor from "../components/AnnotationEditor";
import ConfirmDialog from "../components/ConfirmDialog";
import LiveStatus from "../components/LiveStatus";
import PreflightPanel from "../components/PreflightPanel";

/**
 * Find the annotation a map click landed on.
 *
 * SeqViz's onSelection reports the range clicked, not which annotation it
 * belongs to — a plasmid drawn at a few hundred pixels per thousand bases
 * routinely has several features under one click. The smallest annotation
 * containing the click point is the one a person meant: a 20 bp site inside
 * a 700 bp CDS is what the cursor was actually over.
 */
const COMPLEMENT: Record<string, string> = { A: "T", T: "A", G: "C", C: "G", R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K", B: "V", V: "B", D: "H", H: "D", N: "N" };

/** A reverse-strand feature is read 5'→3' opposite to how it is stored. */
function reverseComplement(seq: string): string {
  return seq
    .toUpperCase()
    .split("")
    .reverse()
    .map((base) => COMPLEMENT[base] ?? base)
    .join("");
}

export function annotationAt(
  annotations: Annotation[], start: number, end: number, sequenceLength = Infinity,
): Annotation | null {
  const covering = annotations.filter((a) => {
    if (a.start <= start && a.end >= end) return true;
    return Number.isFinite(sequenceLength)
      && a.end > sequenceLength
      && start < a.end - sequenceLength && end <= a.end - sequenceLength;
  });
  if (!covering.length) return null;
  return covering.reduce((smallest, a) =>
    a.end - a.start < smallest.end - smallest.start ? a : smallest,
  );
}

export function sequenceForAnnotation(sequence: string, annotation: Annotation): string {
  const span = annotation.end <= sequence.length
    ? sequence.slice(annotation.start, annotation.end)
    : sequence.slice(annotation.start) + sequence.slice(0, annotation.end - sequence.length);
  return annotation.direction === -1 ? reverseComplement(span) : span;
}

type ViewMode = "circular" | "linear" | "both" | "annotated";

export default function Viewer() {
  const { id } = useParams();
  const [project, setProject] = useState<Project | null>(null);
  const [error, setError] = useState("");
  const [actionError, setActionError] = useState("");
  const [mode, setMode] = useState<ViewMode>("circular");
  const [selected, setSelected] = useState<Annotation | null>(null);
  const [status, setStatus] = useState("");
  const [editorOpen, setEditorOpen] = useState(false);
  const [editingIndex, setEditingIndex] = useState<number | null>(null);
  const [initialAnnotation, setInitialAnnotation] = useState<Annotation | null>(null);
  const [savingAnnotation, setSavingAnnotation] = useState(false);
  const [deleteOpen, setDeleteOpen] = useState(false);
  const [detecting, setDetecting] = useState(false);
  const [detected, setDetected] = useState<DetectedFeature[]>([]);
  const [chosenMatches, setChosenMatches] = useState<Set<number>>(new Set());
  const [featureQuery, setFeatureQuery] = useState("");
  const announcedProjectId = useRef<number | null>(null);
  const lifecycle = useRef(0);
  const scanVersion = useRef(0);
  const savePending = useRef(false);

  useEffect(() => {
    let cancelled = false;
    lifecycle.current += 1;
    scanVersion.current += 1;
    savePending.current = false;
    setProject(null); setError(""); setSelected(null); setDetected([]); setActionError("");
    setEditorOpen(false); setDeleteOpen(false); setSavingAnnotation(false);
    setDetecting(false); setChosenMatches(new Set()); setStatus("");
    setFeatureQuery("");
    (async () => {
      try {
        const data = await api.getProject(Number(id));
        if (!cancelled) setProject(data);
      } catch (err) {
        if (!cancelled) {
          setError(err instanceof ApiError ? err.message : "Could not open that project.");
        }
      }
    })();
    return () => {
      cancelled = true;
      lifecycle.current += 1;
      scanVersion.current += 1;
    };
  }, [id]);

  const annotations: Annotation[] = useMemo(
    () => project?.data?.annotations ?? [],
    [project],
  );

  const saveAnnotationList = async (
    next: Annotation[],
    message: string,
    selectedIndex: number | null = null,
  ) => {
    if (!project || savePending.current) return false;
    const version = lifecycle.current;
    savePending.current = true;
    scanVersion.current += 1;
    setDetecting(false);
    setSavingAnnotation(true);
    setActionError("");
    try {
      const updated = await api.updateProjectAnnotations(project.id, next, project.updated_at);
      if (version !== lifecycle.current) return false;
      setProject(updated);
      const saved = updated.data?.annotations ?? [];
      setSelected(selectedIndex === null ? null : saved[selectedIndex] ?? null);
      setStatus(message);
      return true;
    } catch (err) {
      if (version === lifecycle.current) setActionError(err instanceof ApiError ? err.message : "Could not save those annotations.");
      return false;
    } finally {
      if (version === lifecycle.current) {
        savePending.current = false;
        setSavingAnnotation(false);
      }
    }
  };

  const openNewAnnotation = () => {
    setInitialAnnotation(null);
    setEditingIndex(null);
    setEditorOpen(true);
  };

  const detectFeatures = async () => {
    if (!project) return;
    const version = ++scanVersion.current;
    setDetecting(true);
    setActionError("");
    try {
      const result = await api.detectCommonFeatures(project.id);
      if (version !== scanVersion.current) return;
      setDetected(result.matches);
      setChosenMatches(new Set());
      setStatus(
        result.matches.length
          ? `${result.matches.length} exact motif ${result.matches.length === 1 ? "match" : "matches"} found for review.`
          : "No additional curated motif matches were found.",
      );
    } catch (err) {
      if (version === scanVersion.current) setActionError(err instanceof ApiError ? err.message : "Could not scan for common features.");
    } finally {
      if (version === scanVersion.current) setDetecting(false);
    }
  };

  useEffect(() => {
    if (!project) return;
    let cancelled = false;
    const version = ++scanVersion.current;
    setDetecting(true);
    setDetected([]); setChosenMatches(new Set());
    api.detectCommonFeatures(project.id).then((result) => {
      if (!cancelled && version === scanVersion.current) { setDetected(result.matches); setChosenMatches(new Set()); }
    }).catch(() => {
      if (!cancelled && version === scanVersion.current) setActionError("Automatic feature detection failed. Use Find common to retry.");
    }).finally(() => { if (!cancelled && version === scanVersion.current) setDetecting(false); });
    return () => { cancelled = true; };
  }, [project?.id, project?.updated_at]);

  // SeqViz wants its own shape; keep the mapping in one place.
  const seqvizAnnotations = useMemo(
    () =>
      annotations.map((a) => ({
        name: featureLabel(a),
        start: a.start,
        end: a.end,
        direction: a.direction as 1 | -1,
        color: a.color,
      })),
    [annotations],
  );

  useEffect(() => {
    // A circular view of a linear fragment is misleading — follow the record.
    if (project?.data?.topology === "linear") setMode("linear");
  }, [project]);

  useEffect(() => {
    // Set a render after the record lands, not with it: a live region that
    // arrives already holding its sentence is never announced, only one
    // already on the page whose contents then change.
    if (project && announcedProjectId.current !== project.id) {
      announcedProjectId.current = project.id;
      setStatus(
        `${project.name} opened. ${project.sequence.length.toLocaleString()} bases, ` +
        `${project.data?.annotations?.length ?? 0} features.`,
      );
    }
  }, [project]);

  if (error) {
    return (
      <div className="content">
        <div className="notice notice-error" role="alert">{error}</div>
        <p style={{ marginTop: "1rem" }}>
          <Link to="/projects" className="back-link">
            <Icon name="arrowLeft" size={15} /> Back to projects
          </Link>
        </p>
      </div>
    );
  }

  if (!project) {
    return (
      <div className="center-note" role="status" aria-live="polite" aria-busy="true">
        <span className="spinner" />
        <span>Opening sequence…</span>
      </div>
    );
  }

  const topology = project.data?.topology ?? "linear";
  // Normalize GC values across supported project payloads.
  const gc = project.data?.gc_content
    ?? project.data?.construct_gc
    ?? project.data?.gc
    ?? (project.sequence.length
      ? 100 * (project.sequence.match(/[GC]/gi)?.length ?? 0) / project.sequence.length
      : undefined);

  // Restore the complete saved design payload.
  const payload = (project.data ?? {}) as Record<string, unknown>;
  const preflight = payload.preflight as PreflightReport | undefined;
  const provenance = (project.provenance ?? payload.provenance) as Partial<Provenance> | undefined;
  const oligos = (payload.oligos as Record<string, string | number>[]) ?? [];
  const junctions = (payload.junctions as {
    name: string; enzyme: string; kind: string; overhang: string;
    context: string; site_regenerated: boolean;
  }[]) ?? [];
  const protein = typeof payload.protein === "string" ? payload.protein : "";
  const assembly = payload.assembly as { oligos?: Record<string, string | number>[] } | null;
  const allOligos = oligos.length ? oligos : assembly?.oligos ?? [];

  return (
    <>
      <LiveStatus message={status} />

      <div className="topbar">
        <div className="grow">
          <h1>{project.name}</h1>
          <p className="sub">
            {project.notes || `${topology} sequence`}
          </p>
        </div>
        <div className="map-view-switch" role="group" aria-label="Map view">
          {(["circular", "linear", "both", "annotated"] as ViewMode[]).map((m) => (
            <button
              key={m}
              className={`btn ${mode === m ? "btn-primary" : "btn-outline"}`}
              onClick={() => setMode(m)}
              aria-pressed={mode === m}
              disabled={m === "circular" && topology === "linear"}
              title={
                m === "circular" && topology === "linear"
                  ? "This record is linear"
                  : `Show the ${m} view`
              }
            >
              {m === "annotated" ? "Annotated" : m[0].toUpperCase() + m.slice(1)}
            </button>
          ))}
        </div>
      </div>

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
        {actionError && <div className="notice notice-error" role="alert">{actionError}</div>}
        <div className="card">
          <div className="card-body stat-row">
            <div className="stat">
              <div className="k">Length</div>
              <div className="v">
                {project.sequence.length.toLocaleString()}
                <small>bp</small>
              </div>
            </div>
            <div className="stat">
              <div className="k">GC content</div>
              <div className="v">
                {gc !== undefined ? gc.toFixed(1) : "—"}
                <small>%</small>
              </div>
            </div>
            <div className="stat">
              <div className="k">Topology</div>
              <div className="v" style={{ fontSize: "1.05rem", textTransform: "capitalize" }}>
                {topology}
              </div>
            </div>
            <div className="stat">
              <div className="k">Features</div>
              <div className="v">{annotations.length}</div>
            </div>
          </div>
        </div>

        <ExpandablePanel label="Plasmid explorer">
        <div className="viewer-layout">
          <div className={`card seq-stage ${mode === "annotated" ? "seq-stage-annotated" : ""}`}>
            {mode === "annotated" ? (
              <AnnotatedSequenceView
                sequence={project.sequence}
                annotations={annotations}
                selected={selected}
                preferredName={project.name}
                circular={topology === "circular"}
                onAnnotateRange={(range) => {
                  setInitialAnnotation({ ...range, name: "", type: "misc_feature", direction: 1, color: "#3F7A52" });
                  setEditingIndex(null); setEditorOpen(true);
                }}
                onSelect={(annotation) => setSelected(annotation)}
              />
            ) : (
              <SeqViz
              /* The full construct name is already the page heading. SeqViz
                 places its name inside the circular map where it can collide
                 with dense feature labels, so the centre is reserved for the
                 base-pair count. */
              name=""
              seq={project.sequence}
              annotations={seqvizAnnotations}
              viewer={mode}
              showComplement
              showIndex
              disableExternalFonts
              // Clicking a feature in the map selects it here too, so one
              // click either shows the same detail — an annotation is one
              // fact, not two independent views of it.
              onSelection={(sel) => {
                if (sel.type !== "ANNOTATION" || sel.start === undefined || sel.end === undefined) {
                  return;
                }
                const hit = annotationAt(
                  annotations, sel.start, sel.end, project.sequence.length,
                );
                if (hit) setSelected(hit);
              }}
              // The reverse direction: picking a feature from the list
              // highlights its span on the map, because "which one is that"
              // is the question a list of coordinates cannot answer alone.
              highlights={
                selected ? [{ start: selected.start, end: selected.end, color: selected.color }] : []
              }
                style={{ height: "100%", width: "100%" }}
              />
            )}
          </div>

          <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
            <div className="card">
              <div className="card-head">
                <h2 style={{ flex: 1 }}>Features</h2>
                <span className="label">{annotations.length}</span>
                <button
                  className="btn btn-outline"
                  onClick={() => void detectFeatures()}
                  disabled={detecting || savingAnnotation}
                  title="Find exact matches to curated common DNA motifs"
                >
                  {detecting ? "Scanning…" : "Find common"}
                </button>
                <button className="btn btn-primary" onClick={openNewAnnotation}>
                  Add feature
                </button>
              </div>
              {annotations.length > 0 && <div className="card-body field">
                <label htmlFor="feature-search">Find an annotation</label>
                <input id="feature-search" type="search" value={featureQuery}
                  onChange={(event) => setFeatureQuery(event.target.value)}
                  placeholder="Feature name or type" />
              </div>}
              {detected.length > 0 && (
                <div className="motif-review" aria-label="Common motif matches">
                  <div className="motif-review-intro">
                    <strong>Review exact motif matches</strong>
                    <p>
                      Sequence identity supports an annotation, but does not prove biological
                      activity. Choose only the features appropriate for this construct.
                    </p>
                  </div>
                  <div className="motif-match-list">
                    {detected.map((match, index) => {
                      const feature = match.annotation;
                      return (
                        <label className="motif-match" key={`${feature.name}-${feature.start}-${index}`}>
                          <input
                            type="checkbox"
                            checked={chosenMatches.has(index)}
                            onChange={(event) => setChosenMatches((current) => {
                              const next = new Set(current);
                              if (event.target.checked) next.add(index); else next.delete(index);
                              return next;
                            })}
                          />
                          <span className="dot" style={{ background: feature.color }} />
                          <span>
                            <strong>{feature.name}</strong>
                            <small>
                              {(feature.start + 1).toLocaleString()}–
                              {(feature.end > project.sequence.length
                                ? feature.end - project.sequence.length
                                : feature.end).toLocaleString()}
                              {feature.direction === -1 ? " · reverse" : " · forward"}
                            </small>
                            <code>{match.matched_sequence}</code>
                            <small>{match.basis}</small>
                          </span>
                        </label>
                      );
                    })}
                  </div>
                  <div className="motif-review-actions">
                    <button
                      className="btn btn-ghost"
                      onClick={() => {
                        setDetected([]);
                        setChosenMatches(new Set());
                      }}
                    >
                      Dismiss
                    </button>
                    <button
                      className="btn btn-primary"
                      disabled={!chosenMatches.size || savingAnnotation}
                      onClick={() => {
                        const additions = detected
                          .filter((_, index) => chosenMatches.has(index))
                          .map((match) => ({ ...match.annotation, inferred: true, basis: match.basis }));
                        const firstAdded = annotations.length;
                        void saveAnnotationList(
                          [...annotations, ...additions],
                          `${additions.length} reviewed ${additions.length === 1 ? "feature" : "features"} added.`,
                          firstAdded,
                        ).then((saved) => {
                          if (saved) {
                            setDetected([]);
                            setChosenMatches(new Set());
                          }
                        });
                      }}
                    >
                      Add selected ({chosenMatches.size})
                    </button>
                  </div>
                </div>
              )}
              {annotations.length === 0 ? (
                <div className="card-body" style={{ color: "var(--muted)", fontSize: "0.88rem" }}>
                  This record has no annotated features. Name a new insert with “Add feature”,
                  or review exact matches to curated motifs with “Find common”.
                </div>
              ) : (
                <div className="feature-list">
                  {annotations.filter((a) => `${a.name} ${a.type}`.toLowerCase().includes(featureQuery.trim().toLowerCase())).map((a, index) => (
                    <button
                      key={`${a.name}-${a.start}-${index}`}
                      className="feature-row"
                      onClick={() => setSelected(selected === a ? null : a)}
                      aria-pressed={selected === a}
                      style={{
                        background: selected === a ? "var(--accent-wash)" : "transparent",
                        border: "none",
                        borderBottom: "1px solid var(--line)",
                        textAlign: "left",
                        cursor: "pointer",
                        font: "inherit",
                        width: "100%",
                      }}
                    >
                      <span className="dot" style={{ background: a.color }} />
                      <span style={{ minWidth: 0 }}>
                        <span className="nm" style={{ display: "block" }}>{a.name}</span>
                        <span className="ty">
                          {a.type} · {a.direction === -1 ? "reverse" : a.direction === 1 ? "forward" : "unstranded"}{a.inferred ? " · candidate" : ""}
                        </span>
                      </span>
                      <span className="rg">
                        {(a.start + 1).toLocaleString()}–{((a.end - 1) % project.sequence.length + 1).toLocaleString()}
                        {a.end > project.sequence.length ? " · across origin" : ""}
                      </span>
                    </button>
                  ))}
                  {!annotations.some((a) => `${a.name} ${a.type}`.toLowerCase().includes(featureQuery.trim().toLowerCase()))
                    && <p className="card-body">No annotations match this search.</p>}
                </div>
              )}
            </div>

            {/* Clicking a feature — here or on the map itself — has to show
                something, or "interactive" is just a highlight with no
                content behind it. This is that content: what the feature
                is, where it sits, and the bases it actually spans. */}
            {selected && (
              <div className="card feature-detail">
                <div className="card-head">
                  <span className="dot" style={{ background: selected.color }} />
                  <h2 style={{ flex: 1 }}>{selected.name}</h2>
                  <button
                    className="btn btn-outline"
                    onClick={() => {
                      const index = annotations.indexOf(selected);
                      if (index >= 0) {
                        setEditingIndex(index);
                        setEditorOpen(true);
                      }
                    }}
                  >
                    Edit
                  </button>
                  <button className="btn btn-danger" onClick={() => setDeleteOpen(true)}>
                    Delete
                  </button>
                  <button
                    className="btn btn-ghost"
                    onClick={() => setSelected(null)}
                    title="Clear selection"
                    aria-label="Clear selection"
                  >
                    <Icon name="cross" size={14} />
                  </button>
                </div>
                <div className="card-body stat-row">
                  <div className="stat">
                    <div className="k">Type</div>
                    <div className="v" style={{ fontSize: "1.05rem" }}>{selected.type}</div>
                  </div>
                  <div className="stat">
                    <div className="k">Strand</div>
                    <div className="v" style={{ fontSize: "1.05rem" }}>
                      {selected.direction === -1 ? "reverse" : selected.direction === 1 ? "forward" : "—"}
                    </div>
                  </div>
                  <div className="stat">
                    <div className="k">Position</div>
                    <div className="v" style={{ fontSize: "1.05rem" }}>
                      {(selected.start + 1).toLocaleString()}–{((selected.end - 1) % project.sequence.length + 1).toLocaleString()}
                      {selected.end > project.sequence.length ? " · across origin" : ""}
                    </div>
                  </div>
                  <div className="stat">
                    <div className="k">Length</div>
                    <div className="v">
                      {(selected.end - selected.start).toLocaleString()}<small>bp</small>
                    </div>
                  </div>
                </div>
                <FeatureEvidence annotation={selected} />
                {selected.truncated && (
                  <p className="note" style={{ padding: "0 1.1rem 0.9rem", color: "var(--amber)" }}>
                    Truncated at the insert junction — this feature ran past the cut and only
                    part of it is on the plasmid.
                  </p>
                )}
                <div className="card-body" style={{ paddingTop: 0 }}>
                  <div className="seq-block">
                    {sequenceForAnnotation(project.sequence, selected)}
                  </div>
                  {selected.direction === -1 && (
                    <p className="note" style={{ marginTop: "0.5rem" }}>
                      Shown 5'→3' on the strand this feature reads from — the reverse complement
                      of that span in the sequence above.
                    </p>
                  )}
                </div>
              </div>
            )}

            {protein && (
              <div className="card">
                <div className="card-head">
                  <h2 style={{ flex: 1 }}>Protein</h2>
                  <span className="label">{protein.length} residues</span>
                </div>
                <div className="card-body">
                  <div className="seq-block">{protein}</div>
                </div>
              </div>
            )}

            {junctions.length > 0 && (
              <div className="card">
                <div className="card-head"><h2>Junctions</h2></div>
                <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.7rem" }}>
                  {junctions.map((j) => (
                    <div key={j.name} className="junction">
                      <div className="junction-head">
                        <strong>{j.name}</strong>
                        <span className="label">
                          {j.enzyme} · {j.kind} {j.overhang || "blunt"}
                        </span>
                        <span className="grow" />
                        <span className={j.site_regenerated ? "pill pill-ok" : "pill"}>
                          {j.site_regenerated ? "site regenerated" : "site lost"}
                        </span>
                      </div>
                      <div className="junction-seq">
                        <span>{j.context.slice(0, 12)}</span>
                        <span className="seam" />
                        <span>{j.context.slice(12)}</span>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {allOligos.length > 0 && (
              <div className="card">
                <div className="card-head">
                  <h2 style={{ flex: 1 }}>Oligos</h2>
                  <span className="label">{allOligos.length}</span>
                </div>
                <div className="table-scroll">
                  <table className="data">
                    <thead>
                      <tr><th>Name</th><th>Sequence (5'→3')</th><th>Length</th><th>Tm</th></tr>
                    </thead>
                    <tbody>
                      {allOligos.map((oligo) => (
                        <tr key={String(oligo.Name)}>
                          <td className="mono">{oligo.Name}</td>
                          <td className="mono seq-cell">{oligo["Sequence (5\'->3\')"]}</td>
                          <td className="num">{oligo["Length (nt)"]}</td>
                          <td className="num">{oligo["Tm (°C)"]}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            )}

            <div className="card">
              <div className="card-head">
                <h2 style={{ flex: 1 }}>Sequence</h2>
                <button
                  className="btn btn-outline"
                  onClick={() => void api.downloadUrl(
                    `/api/projects/${project.id}/export/`,
                    `${project.name.replace(/\s+/g, "_")}.gb`,
                  ).catch(() => setActionError("Could not download the GenBank file. Try again."))}
                >
                  GenBank
                </button>
              </div>
              <div className="card-body">
                <div className="seq-block">{project.sequence}</div>
              </div>
            </div>
          </div>
        </div>

        </ExpandablePanel>

        <div className="viewer-audit-grid">
          <PreflightPanel report={preflight} />

          {provenance?.output_sha256 && (
            <div className="card">
              <div className="card-head">
                <h2 style={{ flex: 1 }}>Reproducibility record</h2>
                <button className="btn btn-outline" onClick={() => {
                  void navigator.clipboard.writeText(JSON.stringify(provenance, null, 2))
                    .then(() => setStatus("Provenance manifest copied."))
                    .catch(() => setActionError("Clipboard access was denied. Select and copy the manifest manually."));
                }}>Copy manifest</button>
              </div>
              <div className="card-body provenance-grid">
                <div><span>Engine</span><strong>{provenance.engine_version ?? "unknown"}</strong></div>
                <div><span>Workflow</span><strong>{provenance.workflow ?? project.module}</strong></div>
                <div><span>Generated</span><strong>{provenance.generated_at ? new Date(provenance.generated_at).toLocaleString() : "not recorded"}</strong></div>
                <div><span>Output SHA-256</span><code title={provenance.output_sha256}>{provenance.output_sha256}</code></div>
                {provenance.vector_sha256 && <div><span>Vector SHA-256</span><code title={provenance.vector_sha256}>{provenance.vector_sha256}</code></div>}
                <div><span>Enzyme table SHA-256</span><code title={provenance.enzyme_table?.sha256}>{provenance.enzyme_table?.sha256 ?? "not recorded"}</code></div>
              </div>
            </div>
          )}
        </div>

        <p>
          <Link to="/projects" className="back-link">
            <Icon name="arrowLeft" size={15} /> Back to projects
          </Link>
        </p>
      </div>

      <AnnotationEditor
        open={editorOpen}
        annotation={editingIndex === null ? null : annotations[editingIndex] ?? null}
        initialAnnotation={initialAnnotation}
        sequenceLength={project.sequence.length}
        circular={topology === "circular"}
        saving={savingAnnotation}
        saveError={actionError}
        onCancel={() => setEditorOpen(false)}
        onSave={(annotation) => {
          const next = [...annotations];
          const index = editingIndex === null ? next.length : editingIndex;
          if (editingIndex === null) next.push(annotation); else next[editingIndex] = annotation;
          void saveAnnotationList(
            next,
            editingIndex === null
              ? `Feature “${annotation.name}” added.`
              : `Feature “${annotation.name}” updated.`,
            index,
          ).then((saved) => {
            if (saved) setEditorOpen(false);
          });
        }}
      />

      <ConfirmDialog
        open={deleteOpen && selected !== null}
        title={`Delete “${selected?.name ?? "feature"}”?`}
        body="This removes the annotation from this project and its future GenBank exports. It does not change any DNA bases."
        confirmLabel="Delete feature"
        onCancel={() => setDeleteOpen(false)}
        onConfirm={() => {
          if (!selected) return;
          const index = annotations.indexOf(selected);
          if (index < 0) return;
          const name = selected.name;
          setDeleteOpen(false);
          void saveAnnotationList(
            annotations.filter((_, annotationIndex) => annotationIndex !== index),
            `Feature “${name}” deleted.`,
          );
        }}
      />
    </>
  );
}
