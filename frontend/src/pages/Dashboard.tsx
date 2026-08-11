import { useCallback, useEffect, useRef, useState } from "react";
import { Link, useNavigate } from "react-router-dom";

import { ApiError, api, type ProjectSummary } from "../api/client";
import ConfirmDialog from "../components/ConfirmDialog";
import Icon from "../components/Icon";
import LiveStatus from "../components/LiveStatus";

function formatDate(iso: string): string {
  const date = new Date(iso);
  return date.toLocaleDateString(undefined, { year: "numeric", month: "short", day: "numeric" });
}

export default function Dashboard() {
  const navigate = useNavigate();
  const [projects, setProjects] = useState<ProjectSummary[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState("");
  const [uploading, setUploading] = useState(false);
  const [dragOver, setDragOver] = useState(false);
  const [pending, setPending] = useState<ProjectSummary | null>(null);
  const [announcement, setAnnouncement] = useState("");
  const fileInput = useRef<HTMLInputElement>(null);
  const heading = useRef<HTMLHeadingElement>(null);

  const load = useCallback(async () => {
    setLoading(true);
    try {
      const page = await api.listProjects();
      setProjects(page.results);
      setError("");
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Could not load your projects.");
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    void load();
  }, [load]);

  const importFile = useCallback(
    async (file: File) => {
      setUploading(true);
      setError("");
      try {
        // Parse on the server (biopython), then persist what came back.
        const record = await api.parseFile(file);
        const project = await api.createProject({
          name: record.name,
          module: "plasmid_visualizer",
          sequence: record.sequence,
          notes: record.description,
          data: {
            annotations: record.annotations,
            topology: record.topology,
            gc_content: record.gc_content,
          },
        });
        navigate(`/projects/${project.id}`);
      } catch (err) {
        setError(err instanceof ApiError ? err.message : "Could not read that file.");
      } finally {
        setUploading(false);
      }
    },
    [navigate],
  );

  async function exportProject(project: ProjectSummary) {
    try {
      await api.downloadUrl(
        `/api/projects/${project.id}/export/`,
        `${project.name.replace(/\s+/g, "_")}.gb`,
      );
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Could not export that.");
    }
  }

  async function remove(project: ProjectSummary) {
    setPending(null);
    try {
      await api.deleteProject(project.id);
      setProjects((current) => current.filter((p) => p.id !== project.id));
      setAnnouncement(`Deleted “${project.name}”.`);
      // The button that opened the question went with the card, so there is
      // nothing to hand focus back to. The heading of the list they are
      // still in is the nearest honest place to leave them.
      heading.current?.focus();
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Could not delete that project.");
    }
  }

  const status = uploading
    ? "Reading your file…"
    : loading
      ? "Loading your projects…"
      : announcement ||
        `${projects.length} ${projects.length === 1 ? "project" : "projects"}.`;

  return (
    <>
      <LiveStatus message={status} />

      <div className="topbar">
        <div className="grow">
          <h1 ref={heading} tabIndex={-1}>Projects</h1>
          <p className="sub">
            Designs, plasmids and imported sequences. Only you can see them.
          </p>
        </div>
        <button
          className="btn btn-primary"
          onClick={() => fileInput.current?.click()}
          disabled={uploading}
        >
          {uploading && <span className="spinner" />}
          {uploading ? "Importing…" : "Import sequence"}
        </button>
      </div>

      <div
        className="content"
        style={{ display: "flex", flexDirection: "column", gap: "1.25rem" }}
        aria-busy={loading || uploading}
      >
        <input
          ref={fileInput}
          type="file"
          accept=".dna,.gb,.gbk,.genbank,.ape,.fa,.fasta,.fna,.seq"
          style={{ display: "none" }}
          onChange={(e) => {
            const file = e.target.files?.[0];
            e.target.value = "";
            if (file) void importFile(file);
          }}
        />

        {error && <div className="notice notice-error" role="alert">{error}</div>}

        <div
          className={`dropzone${dragOver ? " over" : ""}`}
          onDragOver={(e) => {
            e.preventDefault();
            setDragOver(true);
          }}
          onDragLeave={() => setDragOver(false)}
          onDrop={(e) => {
            e.preventDefault();
            setDragOver(false);
            const file = e.dataTransfer.files?.[0];
            if (file) void importFile(file);
          }}
        >
          <strong>Drop a SnapGene, GenBank or FASTA file here</strong>
          <span className="hint">
            .dna · .gb · .gbk · .fasta — features and topology are read automatically
          </span>
        </div>

        {loading ? (
          <div className="center-note" role="status">
            <span className="spinner" />
            <span>Loading your projects…</span>
          </div>
        ) : projects.length === 0 ? (
          <div className="card">
            <div className="empty">
              <Icon name="plate" size={38} className="glyph" />
              <strong>Nothing saved yet</strong>
              <span>
                Save a design or a plasmid from the workspace, or import a
                sequence here.
              </span>
            </div>
          </div>
        ) : (
          <div className="grid-projects">
            {projects.map((project) => (
              <div className="card project-card" key={project.id}>
                <div className="strip" />
                <Link to={`/projects/${project.id}`} className="body" style={{ textDecoration: "none", color: "inherit" }}>
                  <span className="name">{project.name}</span>
                  <span className="meta">
                    <span className="pill">{project.module.replace(/_/g, " ")}</span>
                    <span>{formatDate(project.updated_at)}</span>
                  </span>
                </Link>
                {/* Read out of context — in a list of every button on the
                    page — "Delete" alone does not say what of. */}
                <div className="foot">
                  <Link
                    to={`/projects/${project.id}`}
                    className="btn btn-ghost"
                    aria-label={`Open ${project.name}`}
                  >
                    Open
                  </Link>
                  <button
                    className="btn btn-ghost"
                    onClick={() => void exportProject(project)}
                    title="GenBank, with its features"
                    aria-label={`Export ${project.name}`}
                  >
                    Export
                  </button>
                  <button
                    className="btn btn-danger"
                    onClick={() => setPending(project)}
                    aria-label={`Delete ${project.name}`}
                  >
                    Delete
                  </button>
                </div>
              </div>
            ))}
          </div>
        )}
      </div>

      <ConfirmDialog
        open={pending !== null}
        title={`Delete “${pending?.name ?? ""}”?`}
        body="The sequence, its features and everything saved with it go for good. There is no copy on the server afterwards."
        confirmLabel="Delete project"
        onConfirm={() => pending && void remove(pending)}
        onCancel={() => setPending(null)}
      />
    </>
  );
}
