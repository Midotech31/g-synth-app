import { useEffect, useMemo, useState } from "react";
import { Link, useParams } from "react-router-dom";
import { SeqViz } from "seqviz";

import { ApiError, api, type Annotation, type Project } from "../api/client";

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
        <div className="notice notice-error">{error}</div>
        <p style={{ marginTop: "1rem" }}>
          <Link to="/">← Back to projects</Link>
        </p>
      </div>
    );
  }

  if (!project) {
    return (
      <div className="center-note">
        <span className="spinner" />
        <span>Opening sequence…</span>
      </div>
    );
  }

  const topology = project.data?.topology ?? "linear";
  const gc = project.data?.gc_content;

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
                      onClick={() => setSelected(a)}
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

            <div className="card">
              <div className="card-head">
                <h2>Sequence</h2>
              </div>
              <div className="card-body">
                <div className="seq-block">{project.sequence}</div>
              </div>
            </div>
          </div>
        </div>

        <p>
          <Link to="/">← Back to projects</Link>
        </p>
      </div>
    </>
  );
}
