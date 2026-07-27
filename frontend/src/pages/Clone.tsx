import { useCallback, useEffect, useRef, useState } from "react";
import { SeqViz } from "seqviz";

import {
  ApiError,
  api,
  type Annotation,
  type Catalogue,
  type CloneResult,
  type DesignParams,
} from "../api/client";
import InsertForm from "../components/InsertForm";

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

type Vector = {
  name: string;
  sequence: string;
  annotations: Annotation[];
  circular: boolean;
};

const EMPTY_VECTOR: Vector = {
  name: "vector",
  sequence: "",
  annotations: [],
  circular: true,
};

export default function Clone() {
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [params, setParams] = useState<DesignParams>(DEFAULTS);
  const [vector, setVector] = useState<Vector>(EMPTY_VECTOR);
  const [result, setResult] = useState<CloneResult | null>(null);
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const [saved, setSaved] = useState("");
  const fileInput = useRef<HTMLInputElement>(null);

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

  const setVectorField = useCallback(<K extends keyof Vector>(key: K, value: Vector[K]) => {
    setVector((current) => ({ ...current, [key]: value }));
    setResult(null);
    setSaved("");
  }, []);

  async function importFile(file: File) {
    setError("");
    try {
      const record = await api.parseFile(file);
      setVector({
        name: record.name || file.name,
        sequence: record.sequence,
        annotations: record.annotations,
        circular: record.topology === "circular",
      });
      setResult(null);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Could not read that file.");
    }
  }

  async function runClone(saveAsProject = false) {
    setBusy(true);
    setError("");
    setSaved("");
    try {
      const data = await api.clone({
        ...params,
        vector: vector.sequence,
        vector_name: vector.name,
        vector_annotations: vector.annotations,
        vector_is_circular: vector.circular,
        save_as_project: saveAsProject,
      });
      setResult(data);
      if (data.project_id) setSaved(`Saved to your projects (#${data.project_id}).`);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "The cloning failed.");
      setResult(null);
    } finally {
      setBusy(false);
    }
  }

  const vectorLength = vector.sequence.replace(/[^ACGTacgt]/g, "").length;
  const ready = vectorLength > 0 && params.sequence.trim().length > 0;

  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>Clone into a vector</h1>
          <p className="sub">
            Cut the vector, drop the construct in, and see the plasmid you end
            up with — junctions, reading frame and all.
          </p>
        </div>
        <button
          className="btn btn-primary"
          onClick={() => runClone(false)}
          disabled={busy || !ready}
          title={ready ? "Clone" : "Add a vector and an insert first"}
        >
          {busy && <span className="spinner" />}
          {busy ? "Cloning…" : "Clone"}
        </button>
      </div>

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
        {error && <div className="notice notice-error">{error}</div>}
        {saved && <div className="notice notice-info">{saved}</div>}

        <div className="design-layout">
          {/* ── Inputs ─────────────────────────────────────────────────── */}
          <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
            <div className="card">
              <div className="card-head">
                <h2 style={{ flex: 1 }}>Vector</h2>
                {vectorLength > 0 && (
                  <span className="label">
                    {vectorLength.toLocaleString()} bp · {vector.annotations.length} features
                  </span>
                )}
              </div>
              <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.85rem" }}>
                <div className="field">
                  <label htmlFor="vector-name">Name</label>
                  <input
                    id="vector-name"
                    type="text"
                    value={vector.name}
                    onChange={(e) => setVectorField("name", e.target.value)}
                    placeholder="pET-28a(+)"
                  />
                </div>

                <div>
                  <button
                    className="btn btn-outline"
                    onClick={() => fileInput.current?.click()}
                    style={{ width: "100%" }}
                  >
                    Import GenBank or FASTA…
                  </button>
                  <input
                    ref={fileInput}
                    type="file"
                    accept=".gb,.gbk,.genbank,.fa,.fasta,.fna,.txt"
                    style={{ display: "none" }}
                    onChange={(e) => {
                      const file = e.target.files?.[0];
                      if (file) void importFile(file);
                      e.target.value = "";
                    }}
                  />
                  <p className="note" style={{ marginTop: "0.4rem" }}>
                    GenBank keeps the vector's features, so they carry over onto
                    the recombinant map. FASTA gives sequence only.
                  </p>
                </div>

                <div className="field">
                  <label htmlFor="vector-seq">Sequence</label>
                  <textarea
                    id="vector-seq"
                    value={vector.sequence}
                    onChange={(e) => setVectorField("sequence", e.target.value)}
                    rows={5}
                    className="mono"
                    style={{ fontSize: "0.76rem" }}
                    placeholder="Paste the vector sequence, or import a file above."
                  />
                </div>

                <div className="checks">
                  <label>
                    <input
                      type="checkbox"
                      checked={vector.circular}
                      onChange={(e) => setVectorField("circular", e.target.checked)}
                    />
                    Circular plasmid
                  </label>
                </div>
              </div>
            </div>

            <div className="card">
              <div className="card-head"><h2>Insert</h2></div>
              <div className="card-body">
                <InsertForm
                  params={params}
                  catalogue={catalogue}
                  onChange={set}
                  idPrefix="clone-"
                />
              </div>
            </div>
          </div>

          {/* ── Results ────────────────────────────────────────────────── */}
          <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
            {!result ? (
              <div className="card">
                <div className="empty">
                  <div className="glyph">🧫</div>
                  <strong>No plasmid yet</strong>
                  <span>
                    {vectorLength === 0
                      ? "Import or paste a vector to begin."
                      : "Set the insert and its ends, then press Clone."}
                  </span>
                </div>
              </div>
            ) : (
              <>
                <div className={`notice ${result.is_clonable ? "notice-ok" : "notice-error"}`}>
                  {result.is_clonable ? (
                    <>
                      <strong>Clonable.</strong> {result.left_enzyme} and{" "}
                      {result.right_enzyme} each cut {result.vector_name} once, the
                      ends match, and the insert goes in one orientation.
                    </>
                  ) : (
                    <>
                      <strong>Will not clone.</strong> {result.problems.join(" ")}
                    </>
                  )}
                </div>

                <div className="card">
                  <div className="card-body stat-row">
                    <div className="stat">
                      <div className="k">Plasmid</div>
                      <div className="v">{result.length.toLocaleString()}<small>bp</small></div>
                    </div>
                    <div className="stat">
                      <div className="k">Backbone</div>
                      <div className="v">{result.backbone_length.toLocaleString()}<small>bp</small></div>
                    </div>
                    <div className="stat">
                      <div className="k">Insert</div>
                      <div className="v">{result.insert_length}<small>bp</small></div>
                    </div>
                    <div className="stat">
                      <div className="k">Removed</div>
                      <div className="v">{result.removed_length}<small>bp</small></div>
                    </div>
                    <div className="stat">
                      <div className="k">Protein</div>
                      <div className="v">{result.protein_length || "—"}<small>aa</small></div>
                    </div>
                  </div>
                </div>

                <div className="card seq-stage" style={{ minHeight: 460 }}>
                  <SeqViz
                    name={result.name}
                    seq={result.plasmid}
                    annotations={result.annotations.map((a) => ({
                      name: a.name,
                      start: a.start,
                      end: a.end,
                      direction: (a.direction === -1 ? -1 : 1) as 1 | -1,
                      color: a.color,
                    }))}
                    viewer="circular"
                    showIndex
                    style={{ height: "100%", width: "100%" }}
                  />
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Junctions</h2>
                    <button
                      className="btn btn-primary"
                      onClick={() => runClone(true)}
                      disabled={busy || !result.is_clonable}
                    >
                      Save plasmid
                    </button>
                  </div>
                  <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.9rem" }}>
                    {result.junctions.map((junction) => (
                      <div key={junction.name} className="junction">
                        <div className="junction-head">
                          <strong>{junction.name}</strong>
                          <span className="label">
                            {junction.enzyme} · {junction.kind} {junction.overhang || "blunt"}
                          </span>
                          <span className="grow" />
                          <span className={junction.site_regenerated ? "pill pill-ok" : "pill"}>
                            {junction.site_regenerated ? "site regenerated" : "site lost"}
                          </span>
                        </div>
                        <div className="junction-seq">
                          <span>{junction.context.slice(0, 12)}</span>
                          <span className="seam" />
                          <span>{junction.context.slice(12)}</span>
                        </div>
                      </div>
                    ))}
                    <p className="note" style={{ margin: 0 }}>
                      A regenerated site means the insert can be cut back out —
                      which is how the clone gets verified on a gel.
                    </p>
                  </div>
                </div>

                {result.protein && (
                  <div className="card">
                    <div className="card-head">
                      <h2 style={{ flex: 1 }}>Protein</h2>
                      <span className="label">{result.protein_length} residues</span>
                    </div>
                    <div className="card-body">
                      <div className="seq-block">{result.protein}</div>
                      <p className="note" style={{ marginTop: "0.6rem" }}>
                        Translated from the insert's ATG through the junction and
                        into the vector, to the first in-frame stop.
                      </p>
                    </div>
                  </div>
                )}

                {(result.warnings.length > 0 || result.insert.warnings.length > 0) && (
                  <div className="card">
                    <div className="card-head"><h2>Notes</h2></div>
                    <div className="card-body">
                      <ul style={{ margin: 0, paddingLeft: "1.1rem", color: "var(--ink-soft)" }}>
                        {[...new Set([...result.insert.warnings, ...result.warnings])].map((note) => (
                          <li key={note} style={{ marginBottom: "0.35rem", lineHeight: 1.5 }}>
                            {note}
                          </li>
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
