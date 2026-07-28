import { useCallback, useEffect, useRef, useState } from "react";
import { SeqViz } from "seqviz";

import {
  ApiError,
  api,
  type Annotation,
  type Catalogue,
  type CloneResult,
  type DesignParams,
  type VectorSpec,
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
  key: string;
  name: string;
  sequence: string;
  annotations: Annotation[];
  circular: boolean;
  /** Set when the sequence came from the catalogue rather than an import. */
  bundled: boolean;
};

const EMPTY_VECTOR: Vector = {
  key: "",
  name: "",
  sequence: "",
  annotations: [],
  circular: true,
  bundled: false,
};

export default function Clone() {
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [params, setParams] = useState<DesignParams>(DEFAULTS);
  const [vectors, setVectors] = useState<VectorSpec[]>([]);
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

  // Load the vector list, then the default vector's own sequence, so the
  // page is usable without importing anything.
  useEffect(() => {
    let cancelled = false;
    (async () => {
      try {
        const list = await api.vectors();
        if (cancelled) return;
        setVectors(list.vectors);
        await selectVector(list.default, list.vectors);
      } catch {
        if (!cancelled) setError("Could not load the vector catalogue.");
      }
    })();
    return () => {
      cancelled = true;
    };
    // selectVector is stable for the page's lifetime.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  /** Switch vector, pulling its bundled sequence when it has one. */
  async function selectVector(key: string, known: VectorSpec[] = vectors) {
    const spec = known.find((v) => v.key === key);
    setResult(null);
    setSaved("");

    if (!spec) {
      setVector({ ...EMPTY_VECTOR, key: "" });
      return;
    }

    // Follow the vector's own cloning pair — pET-21(+) has no NdeI site, so
    // leaving the G-Synth default selected would just fail.
    const pair = spec.recommended_pairs[0]?.split("/").map((p) => p.trim());
    if (pair?.length === 2) {
      setParams((current) => ({
        ...current,
        left_enzyme: pair[0],
        right_enzyme: pair[1],
      }));
    }

    if (!spec.has_sequence) {
      setVector({
        ...EMPTY_VECTOR, key: spec.key, name: spec.name,
      });
      return;
    }
    try {
      const record = await api.vectorSequence(spec.key);
      setVector({
        key: spec.key,
        name: record.name,
        sequence: record.sequence,
        annotations: record.annotations,
        circular: record.topology === "circular",
        bundled: true,
      });
    } catch {
      setError(`Could not load the sequence for ${spec.name}.`);
    }
  }

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
      setVector((current) => ({
        // Keep the catalogue entry selected: the imported sequence is then
        // checked against it, which is how a substitution gets caught.
        key: current.key,
        name: record.name || file.name,
        sequence: record.sequence,
        annotations: record.annotations,
        circular: record.topology === "circular",
        bundled: false,
      }));
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
        vector_key: vector.key,
        // A bundled sequence is already on the server; sending it back would
        // just be a megabyte of round trip.
        vector: vector.bundled ? "" : vector.sequence,
        vector_name: vector.name,
        vector_annotations: vector.bundled ? undefined : vector.annotations,
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
  const spec = vectors.find((v) => v.key === vector.key) ?? null;

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
                  <label htmlFor="vector-key">Backbone</label>
                  <select
                    id="vector-key"
                    value={vector.key}
                    onChange={(e) => void selectVector(e.target.value)}
                  >
                    {vectors.map((v) => (
                      <option key={v.key} value={v.key}>
                        {v.name} · {v.length.toLocaleString()} bp · {v.resistance}
                        {v.has_sequence ? "" : " (import needed)"}
                      </option>
                    ))}
                    <option value="">Something else — I'll supply it</option>
                  </select>
                  {spec && <span className="label">{spec.tag_summary}</span>}
                </div>

                {spec && (
                  <div className="vector-brief">
                    <p className="note" style={{ margin: 0 }}>{spec.summary}</p>
                    <div className="vector-facts">
                      <span><b>{spec.promoter}</b> promoter</span>
                      <span><b>{spec.resistance}</b></span>
                      <span>cloned with <b>{spec.recommended_pairs[0]}</b></span>
                    </div>
                    {spec.notes.map((note) => (
                      <p key={note} className="note vector-note">{note}</p>
                    ))}
                  </div>
                )}

                {!vector.key && (
                  <div className="field">
                    <label htmlFor="vector-name">Name</label>
                    <input
                      id="vector-name"
                      type="text"
                      value={vector.name}
                      onChange={(e) => setVectorField("name", e.target.value)}
                      placeholder="pLab-01"
                    />
                  </div>
                )}

                <div>
                  <button
                    className="btn btn-outline"
                    onClick={() => fileInput.current?.click()}
                    style={{ width: "100%" }}
                  >
                    {vector.bundled ? "Use my own copy instead…" : "Import GenBank or FASTA…"}
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
                    {vector.bundled
                      ? "Using the sequence that ships with G-Synth. Import your lab's own copy if it differs — it will be checked against this entry."
                      : "GenBank keeps the vector's features, so they carry over onto the recombinant map. FASTA gives sequence only."}
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
                {result.vector.check && !result.vector.check.matches && (
                  <div className="notice notice-error">
                    <strong>This is not {result.vector.spec?.name}.</strong>{" "}
                    {result.vector.check.problems.join(" ")}
                  </div>
                )}
                {result.vector.check?.notes.length ? (
                  <div className="notice notice-info">
                    {result.vector.check.notes.join(" ")}
                  </div>
                ) : null}

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
                      {result.tags.length > 0 && (
                        <div className="tag-outcomes">
                          {result.tags.map((tag) => (
                            <div
                              key={`${tag.name}-${tag.end}`}
                              className={tag.present ? "tag-row on" : "tag-row"}
                            >
                              <span className="pill">
                                {tag.end}-term {tag.name}
                              </span>
                              <span>
                                {tag.present
                                  ? `on the protein at residue ${tag.position}`
                                  : "not on this protein"}
                              </span>
                            </div>
                          ))}
                        </div>
                      )}
                      <p className="note" style={{ marginTop: "0.6rem" }}>
                        Translated from the insert's ATG through the junction and
                        into the vector, to the first in-frame stop.
                        {result.reversed_insert &&
                          " The insert reads on the minus strand of the vector's own numbering, as every pET cassette does."}
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
