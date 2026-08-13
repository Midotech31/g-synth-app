import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useLocation } from "react-router-dom";
import { SeqViz } from "seqviz";

import {
  ApiError,
  api,
  type Annotation,
  type Catalogue,
  type CloneParams,
  type CloneResult,
  type DesignParams,
  type VectorSpec,
} from "../api/client";
import InsertForm from "../components/InsertForm";
import JunctionDuplex from "../components/JunctionDuplex";
import Icon from "../components/Icon";
import LiveStatus from "../components/LiveStatus";

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

/** An insert that arrived from the PCR page already cut. Both strands are
 *  carried because the stagger between them is the overhang: from one strand
 *  alone half the geometry is invisible, and an insert cut for a different
 *  enzyme pair would look correct. */
type PreDigested = {
  top: string;
  bottom: string;
  leftEnzyme: string | null;
  rightEnzyme: string | null;
};

export default function Clone() {
  const location = useLocation() as { state?: { preDigested?: PreDigested } | null };
  const [preDigested, setPreDigested] = useState<PreDigested | null>(
    location.state?.preDigested ?? null,
  );
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [params, setParams] = useState<DesignParams>(DEFAULTS);
  const [vectors, setVectors] = useState<VectorSpec[]>([]);
  const [vector, setVector] = useState<Vector>(EMPTY_VECTOR);
  const [result, setResult] = useState<CloneResult | null>(null);
  const [fragment, setFragment] = useState(true);
  const [mapView, setMapView] = useState<"circular" | "linear" | "both">("circular");
  const [showSites, setShowSites] = useState(true);
  const [onlyUsedSites, setOnlyUsedSites] = useState(false);
  const [selected, setSelected] = useState<{
    name: string; start: number; end: number; direction: number; color: string;
    kind: "feature" | "site"; recognition?: string; cuts?: number; used?: boolean;
  } | null>(null);
  const [showEnds, setShowEnds] = useState(true);
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
      const data = await api.clone(clonePayload(saveAsProject));
      setResult(data);
      setSelected(null);        // a stale selection would name a feature from the last plasmid
      if (data.project_id) setSaved(`Saved to your projects (#${data.project_id}).`);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "The cloning failed.");
      setResult(null);
    } finally {
      setBusy(false);
    }
  }

  /** The exact molecule represented by the current inputs. Clone, save and
   * export all use this builder so a PCR-derived insert cannot silently turn
   * back into the ordinary insert form on one of those paths. */
  function clonePayload(saveAsProject = false): CloneParams {
    return {
      ...params,
      // A cut PCR product is an insert, not a gene: designing one around it
      // would add a second set of sites and tags outside ends that are already
      // sticky.
      ...(preDigested
        ? {
            sequence: preDigested.top,
            insert_reverse: preDigested.bottom,
            pre_digested: true,
            left_enzyme: preDigested.leftEnzyme ?? params.left_enzyme,
            right_enzyme: preDigested.rightEnzyme ?? params.right_enzyme,
          }
        : {}),
      vector_key: vector.key,
      // A bundled sequence is already on the server; sending it back would
      // just be a megabyte of round trip.
      vector: vector.bundled ? "" : vector.sequence,
      vector_name: vector.name,
      vector_annotations: vector.bundled ? undefined : vector.annotations,
      vector_is_circular: vector.circular,
      fragment,
      save_as_project: saveAsProject,
    };
  }

  /** Take the plasmid out of G-Synth: GenBank keeps the features. */
  async function exportPlasmid(filetype: "genbank" | "fasta") {
    const safe = (params.name || "construct").replace(/\s+/g, "_");
    try {
      await api.download(
        `/api/design/clone/export/?filetype=${filetype}`,
        clonePayload(),
        `${safe}.${filetype === "fasta" ? "fasta" : "gb"}`,
      );
    } catch {
      setError("The download failed. Try cloning again first.");
    }
  }

  /** What the map draws: the features, plus restriction sites when asked.
   *
   * A site that straddles the origin has an end past the sequence length,
   * which the viewer cannot place — it stays in the count but is not drawn.
   */
  /** Sites the filters keep, and of those the ones the viewer can place. */
  const chosenSites = useMemo(
    () => (result?.restriction_sites ?? []).filter(
      (site) => (onlyUsedSites ? site.used : true),
    ),
    [result, onlyUsedSites],
  );
  const visibleSites = useMemo(
    () => chosenSites.filter((site) => !site.wraps),
    [chosenSites],
  );

  const mapAnnotations = useMemo(() => {
    const base = (result?.annotations ?? []).map((a) => ({
      name: a.name,
      start: a.start,
      end: a.end,
      direction: (a.direction === -1 ? -1 : 1) as 1 | -1,
      color: a.color,
    }));
    if (!result || !showSites) return base;

    const sites = visibleSites.map((site) => ({
        name: site.name,
        start: site.start,
        end: site.end,
        direction: 1 as 1,
        color: site.color,
      }));
    return [...base, ...sites];
  }, [result, showSites, visibleSites]);

  /** Everything a click could land on, features and sites together, so one
   *  handler can find whichever the cursor was actually over. */
  const clickable = useMemo(() => {
    const features = (result?.annotations ?? []).map((a) => ({
      kind: "feature" as const, name: a.name, start: a.start, end: a.end,
      direction: a.direction, color: a.color,
    }));
    const sites = showSites ? visibleSites.map((s) => ({
      kind: "site" as const, name: s.name, start: s.start, end: s.end,
      direction: 1, color: s.color, recognition: s.recognition,
      cuts: s.cuts, used: s.used,
    })) : [];
    return [...features, ...sites];
  }, [result, showSites, visibleSites]);

  const vectorLength = vector.sequence.replace(/[^ACGTacgt]/g, "").length;
  const insertReady = preDigested
    ? preDigested.top.trim().length > 0 && preDigested.bottom.trim().length > 0
    : params.sequence.trim().length > 0;
  const ready = vectorLength > 0 && insertReady;
  const spec = vectors.find((v) => v.key === vector.key) ?? null;

  const status = busy
    ? "Cloning…"
    : result === null
      ? ""
      : result.is_clonable
        ? `Clonable: ${result.length.toLocaleString()} bp plasmid, ${result.validation.filter((c) => c.passed).length} of ${result.validation.length} checks passed.`
        : "This will not clone. Read the reasons above the result.";

  return (
    <>
      <LiveStatus message={status} />

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

      <div
        className="content"
        style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}
        aria-busy={busy}
      >
        {error && <div className="notice notice-error" role="alert">{error}</div>}
        {saved && <div className="notice notice-info" role="status">{saved}</div>}

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
                    {vector.bundled ? "Use my own copy instead…" : "Import SnapGene, GenBank or FASTA…"}
                  </button>
                  <input
                    ref={fileInput}
                    type="file"
                    accept=".dna,.gb,.gbk,.genbank,.fa,.fasta,.fna,.txt"
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
                      : "SnapGene .dna and GenBank both keep the vector's features, so they carry over onto the recombinant map. FASTA gives sequence only."}
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
              {preDigested && (
                /* Without this the page silently ignores the insert form, and
                   the reader has no way to tell why their edits do nothing. */
                <div className="notice notice-info" style={{ margin: "0 1.1rem", marginTop: "0.9rem" }}>
                  <strong>Using a cut PCR product.</strong>{" "}
                  {preDigested.top.length} bp with {preDigested.leftEnzyme} and{" "}
                  {preDigested.rightEnzyme} ends, ligated as supplied &mdash; the
                  options below are not applied to it.{" "}
                  <button
                    className="btn btn-ghost"
                    style={{ padding: "0.1rem 0.4rem", fontSize: "0.8rem" }}
                    onClick={() => {
                      setPreDigested(null);
                      setResult(null);
                      setSelected(null);
                      setSaved("");
                      setError("");
                    }}
                  >
                    Design an insert instead
                  </button>
                </div>
              )}
              <div className="card-body">
                <InsertForm
                  params={params}
                  catalogue={catalogue}
                  onChange={set}
                  showFragmentation={fragment}
                  idPrefix="clone-"
                />
                <div className="checks">
                  <label>
                    <input type="checkbox" checked={!fragment}
                           onChange={(e) => { setFragment(!e.target.checked); setResult(null); }} />
                    Clone the SSD duplex as it is, without fragmenting it
                  </label>
                </div>
              </div>
            </div>
          </div>

          {/* ── Results ────────────────────────────────────────────────── */}
          <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
            {!result ? (
              <div className="card">
                <div className="empty">
                  <Icon name="plate" size={38} className="glyph" />
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

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Map</h2>
                    <div className="seg-toggle" role="group" aria-label="Map view">
                      {(["circular", "linear", "both"] as const).map((option) => (
                        <button key={option} type="button"
                                className={mapView === option ? "on" : ""}
                                aria-pressed={mapView === option}
                                onClick={() => setMapView(option)}>
                          {option[0].toUpperCase() + option.slice(1)}
                        </button>
                      ))}
                    </div>
                  </div>
                  <div className="card-body" style={{ paddingBottom: "0.5rem" }}>
                    <div className="map-filters">
                      <label>
                        <input type="checkbox" checked={showSites}
                               onChange={(e) => setShowSites(e.target.checked)} />
                        Restriction sites
                      </label>
                      <label>
                        <input type="checkbox" checked={onlyUsedSites}
                               disabled={!showSites}
                               onChange={(e) => setOnlyUsedSites(e.target.checked)} />
                        Only the pair used for cloning
                      </label>
                      {showSites && (
                        <span className="label">
                          {visibleSites.length} of {chosenSites.length} sites drawn
                          {chosenSites.length > visibleSites.length && (
                            <>
                              {" · "}
                              {chosenSites.length - visibleSites.length} spans the
                              origin
                            </>
                          )}
                        </span>
                      )}
                    </div>
                  </div>
                  <div className="seq-stage" style={{ minHeight: 460 }}>
                    <SeqViz
                      name={result.name}
                      seq={result.plasmid}
                      annotations={mapAnnotations}
                      viewer={mapView}
                      showIndex
                      onSelection={(sel) => {
                        if (sel.type !== "ANNOTATION" || sel.start === undefined || sel.end === undefined) {
                          return;
                        }
                        const covering = clickable.filter(
                          (a) => a.start <= sel.start! && a.end >= sel.end!,
                        );
                        if (!covering.length) return;
                        const hit = covering.reduce((smallest, a) =>
                          a.end - a.start < smallest.end - smallest.start ? a : smallest,
                        );
                        setSelected(hit);
                      }}
                      highlights={
                        selected ? [{ start: selected.start, end: selected.end, color: selected.color }] : []
                      }
                      style={{ height: "100%", width: "100%" }}
                    />
                  </div>

                  {/* What was clicked. A restriction site and a vector
                      feature answer different questions — a site's recognition
                      sequence and cut count matter more than its strand — so
                      the panel shows what fits the kind rather than one shape
                      forced onto both. */}
                  {selected && (
                    <div className="card-body feature-detail" style={{ borderTop: "1px solid var(--line)" }}>
                      <div className="card-head" style={{ padding: 0, border: "none", marginBottom: "0.6rem" }}>
                        <span className="dot" style={{ background: selected.color }} />
                        {/* Inside the map card, under its own heading. */}
                        <h3 style={{ flex: 1, fontSize: "1rem" }}>
                          {selected.name}
                          {selected.kind === "site" && (
                            <span className="label" style={{ marginLeft: "0.5rem" }}>restriction site</span>
                          )}
                        </h3>
                        <button className="btn btn-ghost" onClick={() => setSelected(null)}
                                title="Clear selection" aria-label="Clear selection">
                          <Icon name="cross" size={14} />
                        </button>
                      </div>
                      <div className="stat-row">
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
                        {selected.kind === "feature" ? (
                          <div className="stat">
                            <div className="k">Strand</div>
                            <div className="v" style={{ fontSize: "1.05rem" }}>
                              {selected.direction === -1 ? "reverse" : "forward"}
                            </div>
                          </div>
                        ) : (
                          <>
                            <div className="stat">
                              <div className="k">Recognition</div>
                              <div className="v mono" style={{ fontSize: "1.05rem" }}>
                                {selected.recognition}
                              </div>
                            </div>
                            <div className="stat">
                              <div className="k">Cuts here</div>
                              <div className="v" style={{ fontSize: "1.05rem" }}>
                                {selected.cuts}{selected.used ? " · used for cloning" : ""}
                              </div>
                            </div>
                          </>
                        )}
                      </div>
                      <div className="seq-block" style={{ marginTop: "0.6rem" }}>
                        {selected.direction === -1
                          ? result.plasmid.slice(selected.start, selected.end)
                              .split("").reverse()
                              .map((b) => ({ A: "T", T: "A", G: "C", C: "G" } as Record<string, string>)[b] ?? b)
                              .join("")
                          : result.plasmid.slice(selected.start, selected.end)}
                      </div>
                    </div>
                  )}
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Checks</h2>
                    <span className="label">
                      {result.validation.filter((c) => c.passed).length}
                      {" of "}{result.validation.length} passed
                    </span>
                  </div>
                  <div className="card-body">
                    <ul className="check-list">
                      {result.validation.map((row) => (
                        <li key={row.check} className={row.passed ? "ok" : "bad"}>
                          <Icon name={row.passed ? "check" : "cross"} size={16} className="mark" />
                          <span>
                            <strong>{row.check}</strong>
                            <span className="detail">{row.detail}</span>
                          </span>
                        </li>
                      ))}
                    </ul>
                  </div>
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Junctions</h2>
                    <label className="inline-check">
                      <input type="checkbox" checked={showEnds}
                             onChange={(e) => setShowEnds(e.target.checked)} />
                      Show the ends before ligation
                    </label>
                    <button className="btn btn-outline"
                            onClick={() => void exportPlasmid("genbank")}
                            disabled={!result.is_clonable}
                            title="Opens in SnapGene, Benchling or ApE with its features">
                      GenBank
                    </button>
                    <button className="btn btn-outline"
                            onClick={() => void exportPlasmid("fasta")}
                            disabled={!result.is_clonable}>
                      FASTA
                    </button>
                    <button
                      className="btn btn-primary"
                      onClick={() => runClone(true)}
                      disabled={busy || !result.is_clonable}
                    >
                      Save plasmid
                    </button>
                  </div>
                  <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
                    {result.junction_views.map((view) => (
                      <JunctionDuplex key={view.name} view={view} showEnds={showEnds} />
                    ))}
                    <p className="note" style={{ margin: 0 }}>
                      Coloured bases are the overhang. Both pieces carry it, on
                      opposite strands — that is what lets them anneal, and a
                      site that survives the join is how the clone gets
                      verified on a gel.
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

                {(result.warnings.length > 0 || (result.insert?.warnings.length ?? 0) > 0) && (
                  <div className="card">
                    <div className="card-head"><h2>Notes</h2></div>
                    <div className="card-body">
                      <ul style={{ margin: 0, paddingLeft: "1.1rem", color: "var(--ink-soft)" }}>
                        {[...new Set([...(result.insert?.warnings ?? []), ...result.warnings])].map((note) => (
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
