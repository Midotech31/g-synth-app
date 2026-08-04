import { useCallback, useEffect, useRef, useState } from "react";
import { Link, useNavigate } from "react-router-dom";

import { ApiError, api, type ProjectSummary } from "../api/client";
import Icon from "../components/Icon";

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
  const fileInput = useRef<HTMLInputElement>(null);

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

  async function remove(id: number, name: string) {
    if (!confirm(`Delete “${name}”? This cannot be undone.`)) return;
    try {
      await api.deleteProject(id);
      setProjects((current) => current.filter((p) => p.id !== id));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Could not delete that project.");
    }
  }

  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>Projects</h1>
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

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.25rem" }}>
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

        {error && <div className="notice notice-error">{error}</div>}

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
          <div className="center-note">
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
                <div className="foot">
                  <Link to={`/projects/${project.id}`} className="btn btn-ghost">
                    Open
                  </Link>
                  <button
                    className="btn btn-ghost"
                    onClick={() => void exportProject(project)}
                    title="GenBank, with its features"
                  >
                    Export
                  </button>
                  <button className="btn btn-danger" onClick={() => remove(project.id, project.name)}>
                    Delete
                  </button>
                </div>
              </div>
            ))}
          </div>
        )}
      </div>
    </>
  );
}
