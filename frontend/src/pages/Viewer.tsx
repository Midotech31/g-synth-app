import { useEffect, useMemo, useState } from "react";
import { Link, useParams } from "react-router-dom";
import { SeqViz } from "seqviz";

import { ApiError, api, type Annotation, type Project } from "../api/client";
import Icon from "../components/Icon";

/**
 * Find the annotation a map click landed on.
 *
 * SeqViz's onSelection reports the range clicked, not which annotation it
 * belongs to — a plasmid drawn at a few hundred pixels per thousand bases
 * routinely has several features under one click. The smallest annotation
 * containing the click point is the one a person meant: a 20 bp site inside
 * a 700 bp CDS is what the cursor was actually over.
 */
const COMPLEMENT: Record<string, string> = { A: "T", T: "A", G: "C", C: "G", N: "N" };

/** A reverse-strand feature is read 5'→3' opposite to how it is stored. */
function reverseComplement(seq: string): string {
  return seq
    .toUpperCase()
    .split("")
    .reverse()
    .map((base) => COMPLEMENT[base] ?? base)
    .join("");
}

function annotationAt(annotations: Annotation[], start: number, end: number): Annotation | null {
  const covering = annotations.filter((a) => a.start <= start && a.end >= end);
  if (!covering.length) return null;
  return covering.reduce((smallest, a) =>
    a.end - a.start < smallest.end - smallest.start ? a : smallest,
  );
}

type ViewMode = "circular" | "linear" | "both";

export default function Viewer() {
  const { id } = useParams();
  const [project, setProject] = useState<Project | null>(null);
  const [error, setError] = useState("");
  const [mode, setMode] = useState<ViewMode>("circular");
  const [selected, setSelected] = useState<Annotation | null>(null);

  useEffect(() => {
    let cancelled = false;
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
    };
  }, [id]);

  const annotations: Annotation[] = useMemo(
    () => project?.data?.annotations ?? [],
    [project],
  );

  // SeqViz wants its own shape; keep the mapping in one place.
  const seqvizAnnotations = useMemo(
    () =>
      annotations.map((a) => ({
        name: a.name,
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
  // Imports, assemblies and cloned plasmids use different payload names.
  // Keep old saved projects readable, and derive the value as a last resort.
  const gc = project.data?.gc_content
    ?? project.data?.construct_gc
    ?? project.data?.gc
    ?? (project.sequence.length
      ? 100 * (project.sequence.match(/[GC]/gi)?.length ?? 0) / project.sequence.length
      : undefined);

  // A saved design carries its whole payload. Showing only the map made the
  // oligos, the junctions and the protein unrecoverable — the parts someone
  // reopens a project *for*.
  const payload = (project.data ?? {}) as Record<string, unknown>;
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
      <div className="topbar">
        <div className="grow">
          <h1>{project.name}</h1>
          <p className="sub">
            {project.notes || `${topology} sequence`}
          </p>
        </div>
        <div style={{ display: "flex", gap: "0.4rem" }}>
          {(["circular", "linear", "both"] as ViewMode[]).map((m) => (
            <button
              key={m}
              className={`btn ${mode === m ? "btn-primary" : "btn-outline"}`}
              onClick={() => setMode(m)}
              disabled={m === "circular" && topology === "linear"}
              title={
                m === "circular" && topology === "linear"
                  ? "This record is linear"
                  : `Show the ${m} view`
              }
            >
              {m[0].toUpperCase() + m.slice(1)}
            </button>
          ))}
        </div>
      </div>

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
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

        <div className="viewer-layout">
          <div className="card seq-stage">
            <SeqViz
              name={project.name}
              seq={project.sequence}
              annotations={seqvizAnnotations}
              viewer={mode}
              showComplement
              showIndex
              // Clicking a feature in the map selects it here too, so one
              // click either shows the same detail — an annotation is one
              // fact, not two independent views of it.
              onSelection={(sel) => {
                if (sel.type !== "ANNOTATION" || sel.start === undefined || sel.end === undefined) {
                  return;
                }
                const hit = annotationAt(annotations, sel.start, sel.end);
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
          </div>

          <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
            <div className="card">
              <div className="card-head">
                <h2 style={{ flex: 1 }}>Features</h2>
                <span className="label">{annotations.length}</span>
              </div>
              {annotations.length === 0 ? (
                <div className="card-body" style={{ color: "var(--muted)", fontSize: "0.88rem" }}>
                  This record has no annotated features. FASTA files carry sequence only —
                  import a GenBank file to see a feature map.
                </div>
              ) : (
                <div className="feature-list">
                  {annotations.map((a, index) => (
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
                          {a.type} · {a.direction === -1 ? "reverse" : "forward"}
                        </span>
                      </span>
                      <span className="rg">
                        {(a.start + 1).toLocaleString()}–{a.end.toLocaleString()}
                      </span>
                    </button>
                  ))}
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
                    className="btn btn-ghost"
                    onClick={() => setSelected(null)}
                    title="Clear selection"
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
                      {(selected.start + 1).toLocaleString()}–{selected.end.toLocaleString()}
                    </div>
                  </div>
                  <div className="stat">
                    <div className="k">Length</div>
                    <div className="v">
                      {(selected.end - selected.start).toLocaleString()}<small>bp</small>
                    </div>
                  </div>
                </div>
                {selected.truncated && (
                  <p className="note" style={{ padding: "0 1.1rem 0.9rem", color: "var(--amber)" }}>
                    Truncated at the insert junction — this feature ran past the cut and only
                    part of it is on the plasmid.
                  </p>
                )}
                <div className="card-body" style={{ paddingTop: 0 }}>
                  <div className="seq-block">
                    {selected.direction === -1
                      ? reverseComplement(project.sequence.slice(selected.start, selected.end))
                      : project.sequence.slice(selected.start, selected.end)}
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
                  )}
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

        <p>
          <Link to="/projects" className="back-link">
            <Icon name="arrowLeft" size={15} /> Back to projects
          </Link>
        </p>
      </div>
    </>
  );
}
