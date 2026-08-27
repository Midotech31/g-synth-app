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
import PreflightPanel from "../components/PreflightPanel";
import { useWorkspaceState } from "../state/WorkspaceStateContext";

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
  const [preDigested, setPreDigested, clearPreDigested] = useWorkspaceState<PreDigested | null>("clone.preDigested", null);
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [params, setParams, clearParams] = useWorkspaceState<DesignParams>("clone.params", DEFAULTS);
  const [vectors, setVectors] = useState<VectorSpec[]>([]);
  const [vector, setVector] = useWorkspaceState<Vector>("clone.vector", EMPTY_VECTOR);
  const [vectorLoaded, setVectorLoaded] = useWorkspaceState("clone.vectorLoaded", false);
  const [result, setResult, clearResult] = useWorkspaceState<CloneResult | null>("clone.result", null);
  const [fragment, setFragment, clearFragment] = useWorkspaceState("clone.fragment", true);
  const [mapView, setMapView, clearMapView] = useWorkspaceState<"circular" | "linear" | "both">("clone.mapView", "circular");
  const [showSites, setShowSites, clearShowSites] = useWorkspaceState("clone.showSites", true);
  const [onlyUsedSites, setOnlyUsedSites, clearOnlyUsedSites] = useWorkspaceState("clone.onlyUsedSites", false);
  const [selected, setSelected, clearSelected] = useWorkspaceState<{
    name: string; start: number; end: number; direction: number; color: string;
    kind: "feature" | "site"; recognition?: string; cuts?: number; used?: boolean;
  } | null>("clone.selected", null);
  const [showEnds, setShowEnds, clearShowEnds] = useWorkspaceState("clone.showEnds", true);
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const [saved, setSaved, clearSaved] = useWorkspaceState("clone.saved", "");
  const [experience, setExperience] = useWorkspaceState<"guided" | "expert">("clone.experience", "guided");
  const fileInput = useRef<HTMLInputElement>(null);

  useEffect(() => {
    api.catalogue().then(setCatalogue).catch(() => {
      setError("Could not load the enzyme catalogue.");
    });
  }, []);

  useEffect(() => {
    if (!location.state?.preDigested) return;
    setPreDigested(location.state.preDigested);
    setResult(null);
  }, [location.state?.preDigested, setPreDigested, setResult]);

  // Load the vector list, then the default vector's own sequence, so the
  // page is usable without importing anything.
  useEffect(() => {
    let cancelled = false;
    (async () => {
      try {
        const list = await api.vectors();
        if (cancelled) return;
        setVectors(list.vectors);
        if (!vectorLoaded) {
          await selectVector(list.default, list.vectors);
          setVectorLoaded(true);
        }
      } catch {
        if (!cancelled) setError("Could not load the vector catalogue.");
      }
    })();
    return () => {
      cancelled = true;
    };
    // Capture the restored flag once per mount. Updating it after the first
    // catalogue load must not fetch the same catalogue a second time.
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

  function clearWorkspace() {
    clearPreDigested();
    clearParams();
    clearResult();
    clearFragment();
    clearMapView();
    clearShowSites();
    clearOnlyUsedSites();
    clearSelected();
    clearShowEnds();
    clearSaved();
    setError("");
  }

  async function downloadWorksheet() {
    const safe = (params.name || "construct").replace(/\s+/g, "_");
    try {
      await api.download(
        "/api/design/clone/worksheet/",
        clonePayload(),
        `${safe}_bench_worksheet.txt`,
      );
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "The bench worksheet could not be generated.");
    }
  }

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
        <button className="btn btn-outline" onClick={clearWorkspace} disabled={busy}>
          Clear
        </button>
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
      