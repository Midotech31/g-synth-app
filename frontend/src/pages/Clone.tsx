import { useCallback, useEffect, useRef, useState } from "react";
import { useLocation } from "react-router-dom";

import {
  ApiError,
  api,
  type Annotation,
  type Catalogue,
  type CloneParams,
  type CloneResult,
  type DesignParams,
  type JunctionView,
  type ValidationCheck,
  type VectorSpec,
} from "../api/client";
import InsertSettingsSummary from "../components/InsertSettingsSummary";
import VectorConfiguration from "../components/VectorConfiguration";
import InsertForm from "../components/InsertForm";
import ConstructWorkbench from "../components/ConstructWorkbench";
import CoreWorkflowTrail from "../components/CoreWorkflowTrail";
import EnzymePicker from "../components/EnzymePicker";
import JunctionDuplex from "../components/JunctionDuplex";
import LigationOutcome from "../components/LigationOutcome";
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
  orfStart?: number | null;
  insertAnnotations?: Annotation[];
  name?: string;
  origin?: "design" | "hybridization" | "pcr";
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
  const [productAnnotations, setProductAnnotations, clearProductAnnotations] = useWorkspaceState<Annotation[] | null>("clone.productAnnotations", null);
  const [fragment, setFragment, clearFragment] = useWorkspaceState("clone.fragment", true);
  const [assemblyDetail, setAssemblyDetail, clearAssemblyDetail] = useWorkspaceState<"simple" | "detailed">("clone.assemblyDetail", "simple");
  const [ligationCommitted, setLigationCommitted, clearLigationCommitted] = useWorkspaceState("clone.ligationCommitted", false);
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const [saved, setSaved, clearSaved] = useWorkspaceState("clone.saved", "");
  const [experience, setExperience] = useWorkspaceState<"guided" | "expert">("clone.experience", "guided");
  const fileInput = useRef<HTMLInputElement>(null);

  const requestVersion = useRef(0);
  const vectorRequestVersion = useRef(0);
  const [vectorLoading, setVectorLoading] = useState(false);

  const invalidate = useCallback(() => {
    requestVersion.current += 1;
    setResult(null);
    setProductAnnotations(null);
    setLigationCommitted(false);
    setSaved("");
    setError("");
    setBusy(false);
  }, [setResult, setProductAnnotations, setLigationCommitted, setSaved]);

  useEffect(() => {
    api.catalogue().then(setCatalogue).catch(() => {
      setError("Could not load the enzyme catalogue.");
    });
  }, []);

  useEffect(() => {
    if (!location.state?.preDigested) return;
    setPreDigested(location.state.preDigested);
    setParams((current) => ({
      ...current,
      name: location.state?.preDigested?.name ?? current.name,
      left_enzyme: location.state?.preDigested?.leftEnzyme ?? current.left_enzyme,
      right_enzyme: location.state?.preDigested?.rightEnzyme ?? current.right_enzyme,
    }));
    invalidate();
  }, [location.state?.preDigested, invalidate, setParams, setPreDigested]);

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
    const version = ++vectorRequestVersion.current;
    invalidate();
    setVectorLoading(false);
    setVector({ ...EMPTY_VECTOR, key: spec?.key ?? "", name: spec?.name ?? "" });

    if (!spec) {
      setVector({ ...EMPTY_VECTOR, key: "" });
      return;
    }

    // Follow the vector's own cloning pair — pET-21(+) has no NdeI site, so
    // leaving the G-Synth default selected would just fail.
    const pair = spec.recommended_pairs[0]?.split("/").map((p) => p.trim());
    if (pair?.length === 2 && !preDigested && !location.state?.preDigested) {
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
    setVectorLoading(true);
    try {
      const record = await api.vectorSequence(spec.key);
      if (vectorRequestVersion.current !== version) return;
      setVector({
        key: spec.key,
        name: record.name,
        sequence: record.sequence,
        annotations: record.annotations,
        circular: record.topology === "circular",
        bundled: true,
      });
    } catch {
      if (vectorRequestVersion.current === version) setError(`Could not load the sequence for ${spec.name}.`);
    } finally {
      if (vectorRequestVersion.current === version) setVectorLoading(false);
    }
  }

  const set = useCallback(
    <K extends keyof DesignParams>(key: K, value: DesignParams[K]) => {
      setParams((current) => ({ ...current, [key]: value }));
      invalidate();
    },
    [invalidate, setParams],
  );

  const setVectorField = useCallback(<K extends keyof Vector>(key: K, value: Vector[K]) => {
    vectorRequestVersion.current += 1;
    setVectorLoading(false);
    setVector((current) => ({
      ...current, [key]: value,
      // Editing bases invalidates imported coordinates and the catalogue shortcut.
      ...(key === "sequence" ? { bundled: false, annotations: [] } : {}),
    }));
    invalidate();
  }, [invalidate, setVector]);

  function setTransferredEnzyme(side: "leftEnzyme" | "rightEnzyme", enzyme: string) {
    setPreDigested((current) => current ? { ...current, [side]: enzyme } : current);
    setParams((current) => ({
      ...current,
      [side === "leftEnzyme" ? "left_enzyme" : "right_enzyme"]: enzyme,
    }));
    invalidate();
  }

  async function importFile(file: File) {
    const version = ++vectorRequestVersion.current;
    invalidate();
    setVectorLoading(true);
    try {
      const record = await api.parseFile(file);
      if (vectorRequestVersion.current !== version) return;
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
    } catch (err) {
      if (vectorRequestVersion.current === version) setError(err instanceof ApiError ? err.message : "Could not read that file.");
    } finally {
      if (vectorRequestVersion.current === version) setVectorLoading(false);
    }
  }

  async function runClone(saveAsProject = false) {
    const version = ++requestVersion.current;
    setBusy(true);
    setError("");
    setSaved("");
    if (!saveAsProject) setLigationCommitted(false);
    try {
      const data = await api.clone(clonePayload(saveAsProject, saveAsProject));
      if (requestVersion.current !== version) return;
      setResult(data);
      setProductAnnotations(data.annotations);
      if (data.project_id) setSaved(`Saved to your projects (#${data.project_id}).`);
    } catch (err) {
      if (requestVersion.current !== version) return;
      setError(err instanceof ApiError ? err.message : "The cloning failed.");
      setResult(null);
    } finally {
      if (requestVersion.current === version) setBusy(false);
    }
  }

  /** The exact molecule represented by the current inputs. Clone, save and
   * export all use this builder so a PCR-derived insert cannot silently turn
   * back into the ordinary insert form on one of those paths. */
  function clonePayload(saveAsProject = false, includeReviewedAnnotations = false): CloneParams {
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
            orf_start: preDigested.orfStart ?? null,
            insert_annotations: preDigested.insertAnnotations,
          }
        : {}),
      vector_key: vector.key,
      // A bundled sequence is already on the server; sending it back would
      // just be a megabyte of round trip.
      vector: vector.bundled ? "" : vector.sequence,
      vector_name: vector.name,
      vector_annotations: vector.bundled ? undefined : vector.annotations,
      product_annotations: includeReviewedAnnotations
        ? productAnnotations ?? result?.annotations
        : undefined,
      vector_is_circular: vector.circular,
      fragment,
      save_as_project: saveAsProject,
    };
  }

  /** Take the plasmid out of G-Synth: GenBank keeps the features. */
  async function exportPlasmid(filetype: "genbank" | "fasta" | "sbol3") {
    const safe = (params.name || "construct").replace(/\s+/g, "_");
    try {
      await api.download(
        `/api/design/clone/export/?filetype=${filetype}`,
        clonePayload(false, true),
        `${safe}.${filetype === "fasta" ? "fasta" : filetype === "sbol3" ? "sbol.json" : "gb"}`,
      );
    } catch {
      setError("The download failed. Try cloning again first.");
    }
  }

  const vectorLength = vector.sequence.replace(/[^ACGTacgt]/g, "").length;
  const insertReady = preDigested
    ? preDigested.top.trim().length > 0 && preDigested.bottom.trim().length > 0
    : params.sequence.trim().length > 0;
  const ready = !vectorLoading && vectorLength > 0 && insertReady;
  const spec = vectors.find((v) => v.key === vector.key) ?? null;
  const validationStatus = (check: ValidationCheck) =>
    check.status ?? (check.passed ? "pass" : "block");
  const validationCounts = result?.validation.reduce(
    (counts, check) => {
      counts[validationStatus(check)] += 1;
      return counts;
    },
    { pass: 0, review: 0, block: 0 },
  ) ?? { pass: 0, review: 0, block: 0 };

  const status = busy
    ? "Checking insert and vector ends…"
    : result === null
      ? ""
      : result.is_clonable
        ? ligationCommitted
          ? `Ligation simulated: ${result.length.toLocaleString()} bp plasmid assembled.`
          : result.preflight?.can_export === false
            ? "The ends are compatible, but expression validation is blocked. Correct the failed check before ligation."
            : `Ends compatible: ${validationCounts.pass} passed${validationCounts.review ? `, ${validationCounts.review} to review` : ""}; ready for simulated ligation.`
        : "This will not clone. Read the reasons above the result.";

  function clearWorkspace() {
    invalidate();
    clearPreDigested();
    clearParams();
    clearResult();
    clearProductAnnotations();
    clearFragment();
    clearAssemblyDetail();
    clearLigationCommitted();
    clearSaved();
    setError("");
  }

  async function downloadWorksheet() {
    const safe = (params.name || "construct").replace(/\s+/g, "_");
    try {
      await api.download(
        "/api/design/clone/worksheet/",
        clonePayload(false, true),
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
          <h1>Restriction enzyme cloning simulation</h1>
          <p className="sub">
            Digest vector and insert in silico, verify both exposed ends and
            orientation, then ligate the compatible product.
          </p>
          <CoreWorkflowTrail active="cloning" />
        </div>
        <button className="btn btn-outline" onClick={clearWorkspace} disabled={busy}>
          Clear
        </button>
        <button
          className="btn btn-primary"
          onClick={() => runClone(false)}
          disabled={busy || !ready}
          title={ready ? "Simulate restriction digestion" : "Add a vector and an insert first"}
        >
          {busy && <span className="spinner" />}
          {busy ? "Digesting…" : result ? "Recheck digestion" : "Simulate digestion"}
        </button>
      </div>

      <div
        className="content"
        style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}
        aria-busy={busy}
      >
        {error && <div className="notice notice-error" role="alert">{error}</div>}
        {saved && <div className="notice notice-info" role="status">{saved}</div>}

        <div className="design-layout clone-layout">
          <div className="clone-inputs">
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
                </div>

                <VectorConfiguration spec={spec}
                  leftEnzyme={preDigested?.leftEnzyme ?? params.left_enzyme}
                  rightEnzyme={preDigested?.rightEnzyme ?? params.right_enzyme}
                  result={result} bundled={vector.bundled} loading={vectorLoading} />

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
                    accept=".dna,.gb,.gbk,.genbank,.fa,.fasta,.fna,.jsonld,.sbol,.ttl,.rdf,.xml,.txt"
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
                      : "SnapGene, GenBank and SBOL 3 preserve features; FASTA supplies sequence only."}
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
              <div className="card-head">
                <h2 style={{ flex: 1 }}>Insert</h2>
                <div className="mode-switch" role="group" aria-label="Cloning detail level">
                  <button className={experience === "guided" ? "active" : ""} onClick={() => {
                    setExperience("guided");
                    setParams((current) => ({ ...current, cleavage_site: "Thrombin", include_his_tag: true, include_linkers: true, remove_stop: false, target_oligo_length: 90, overhang_length: 4 }));
                    invalidate();
                  }}>Guided</button>
                  <button className={experience === "expert" ? "active" : ""} onClick={() => setExperience("expert")}>Expert</button>
                </div>
              </div>
              {preDigested && (
                <div className="notice notice-info" style={{ margin: "0 1.1rem", marginTop: "0.9rem" }}>
                  <strong>
                    {preDigested.origin === "hybridization"
                      ? "Using the duplex validated in Hybridization."
                      : preDigested.origin === "design"
                        ? "Using the duplex produced by Design."
                        : preDigested.origin === "pcr"
                          ? "Using the digested PCR product."
                          : "Using a pre-digested duplex."}
                  </strong>{" "}
                  {preDigested.top.length} bp with{" "}
                  {preDigested.leftEnzyme ?? params.left_enzyme} and{" "}
                  {preDigested.rightEnzyme ?? params.right_enzyme} end assignments.
                  G-Synth will compare the observed
                  strand geometry with the vector before enabling ligation.{" "}
                  <button
                    className="btn btn-ghost"
                    style={{ padding: "0.1rem 0.4rem", fontSize: "0.8rem" }}
                    onClick={() => {
                      setPreDigested(null);
                      invalidate();
                    }}
                  >
                    Change insert or enzymes
                  </button>
                  <div className="transferred-enzyme-grid">
                    <EnzymePicker
                      id="transferred-left-enzyme"
                      label="Left restriction enzyme"
                      enzymes={catalogue?.enzymes ?? []}
                      value={preDigested.leftEnzyme ?? params.left_enzyme}
                      onChange={(value) => setTransferredEnzyme("leftEnzyme", value)}
                      disabled={!catalogue}
                    />
                    <EnzymePicker
                      id="transferred-right-enzyme"
                      label="Right restriction enzyme"
                      enzymes={catalogue?.enzymes ?? []}
                      value={preDigested.rightEnzyme ?? params.right_enzyme}
                      onChange={(value) => setTransferredEnzyme("rightEnzyme", value)}
                      disabled={!catalogue}
                    />
                    <p className="field-hint">
                      Changing an assignment never changes the transferred bases;
                      the simulation must prove that the new cut geometry is compatible.
                    </p>
                  </div>
                </div>
              )}
              <div className="card-body">
                {experience === "guided" && !preDigested && (
                  <InsertSettingsSummary params={params} fragment={fragment} />
                )}
                {!preDigested && (
                  <InsertForm
                    params={params}
                    catalogue={catalogue}
                    onChange={set}
                    showFragmentation={fragment}
                    idPrefix="clone-"
                    expert={experience === "expert"}
                  />
                )}
                {experience === "expert" && <div className="checks">
                  <label>
                    <input type="checkbox" checked={!fragment}
                           onChange={(e) => { setFragment(!e.target.checked); invalidate(); }} />
                    Clone the SSD duplex as it is, without fragmenting it
                  </label>
                </div>}
              </div>
            </div>
          </div>

          <div className="clone-results">
            {!result ? (
              <div className="card">
                <div className="empty">
                  <Icon name="plate" size={38} className="glyph" />
                  <strong>No plasmid yet</strong>
                  <span>
                    {vectorLength === 0
                      ? "Import or paste a vector to begin."
                      : "Set the insert and its ends, then press Simulate digestion."}
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

                <div className={`notice ${result.is_clonable && result.preflight?.can_export !== false ? "notice-ok" : "notice-error"}`}>
                  {result.is_clonable && result.preflight?.can_export !== false ? (
                    <>
                      <strong>{ligationCommitted ? "Ligation simulated." : "Ready to ligate."}</strong>{" "}
                      {result.left_enzyme} and {result.right_enzyme} each cut{" "}
                      {result.vector_name} once, the ends match, and the insert goes
                      in one orientation{ligationCommitted ? "." : "; inspect the ends below before joining them."}
                    </>
                  ) : result.is_clonable ? (
                    <>
                      <strong>Expression validation blocked.</strong>{" "}
                      The DNA ends can ligate, but the construct is not confirmed for expression.
                    </>
                  ) : (
                    <>
                      <strong>Will not clone.</strong> {result.problems.join(" ")}
                    </>
                  )}
                </div>

                <PreflightPanel report={result.preflight} />

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

                <section className={`card frame-assessment frame-${result.reading_frame.status}`} aria-labelledby="frame-assessment-title">
                  <div className="card-head">
                    <div style={{ flex: 1 }}>
                      <h2 id="frame-assessment-title">Expression reading frame</h2>
                      <span className="label">RBS → start codon → junctions → terminus</span>
                    </div>
                    <span className="frame-verdict">
                      {result.reading_frame.status === "pass"
                        ? "Confirmed"
                        : result.reading_frame.status === "block"
                          ? "Invalid"
                          : result.reading_frame.status === "not_applicable"
                            ? "Not applicable"
                            : "Review required"}
                    </span>
                  </div>
                  <div className="card-body frame-assessment-body">
                    <p className="frame-summary">{result.reading_frame.summary}</p>
                    {result.reading_frame.status !== "not_applicable" && (
                      <div className="frame-path" aria-label="Expression cassette frame evidence">
                        <div className={result.reading_frame.promoter_name ? "verified" : "unknown"}>
                          <span>Promoter</span>
                          <strong>{result.reading_frame.promoter_name ?? "Unconfirmed"}</strong>
                          {result.reading_frame.promoter_source === "sequence_motif" && <small>detected motif · review</small>}
                        </div>
                        <Icon name="arrowRight" size={16} />
                        <div className={result.reading_frame.rbs_name ? "verified" : "unknown"}>
                          <span>RBS</span>
                          <strong>{result.reading_frame.rbs_name ?? "Unconfirmed"}</strong>
                          {result.reading_frame.rbs_spacing_nt !== null && (
                            <small>
                              {result.reading_frame.rbs_spacing_nt} nt to ATG
                              {result.reading_frame.rbs_source === "sequence_motif" ? " · detected" : ""}
                            </small>
                          )}
                        </div>
                        <Icon name="arrowRight" size={16} />
                        <div className={result.reading_frame.start_codon === "ATG" ? "verified" : "invalid"}>
                          <span>Start</span>
                          <strong>{result.reading_frame.start_codon ?? "—"}</strong>
                          {result.reading_frame.translation_start !== null && (
                            <small>
                              base {result.reading_frame.translation_start + 1}
                              {result.reading_frame.start_source === "sequence_candidate" ? " · candidate" : ""}
                            </small>
                          )}
                        </div>
                        <Icon name="arrowRight" size={16} />
                        <div className={result.reading_frame.status === "block" ? "invalid" : "verified"}>
                          <span>ORF</span>
                          <strong>{result.reading_frame.protein_length || "—"} aa</strong>
                          {result.reading_frame.right_junction_phase !== null && <small>junction offset {result.reading_frame.right_junction_phase}/3</small>}
                        </div>
                        <Icon name="arrowRight" size={16} />
                        <div className={result.reading_frame.stop_codon ? "verified" : "invalid"}>
                          <span>Stop</span>
                          <strong>{result.reading_frame.stop_codon ?? "Missing"}</strong>
                          {result.reading_frame.stop_context && <small>{result.reading_frame.stop_context.replace("_", " ")}</small>}
                        </div>
                      </div>
                    )}
                    <div className="frame-checks">
                      {result.reading_frame.checks.map((check) => (
                        <div key={check.code} className={check.status}>
                          <Icon name={check.status === "pass" ? "check" : check.status === "review" ? "target" : "cross"} size={15} />
                          <span><strong>{check.label}</strong><small>{check.detail}</small></span>
                        </div>
                      ))}
                    </div>
                  </div>
                </section>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Compatibility checks</h2>
                    <span className="label">
                      {validationCounts.pass} passed
                      {validationCounts.review > 0 && ` · ${validationCounts.review} review`}
                    </span>
                  </div>
                  <div className="card-body">
                    <ul className="check-list">
                      {result.validation.map((check) => {
                        const state = validationStatus(check);
                        return (
                          <li key={check.check} className={state === "pass" ? "ok" : state}>
                            <Icon
                              name={state === "pass" ? "check" : state === "review" ? "target" : "cross"}
                              size={16}
                              className="mark"
                            />
                            <span>
                              <strong>{check.check}</strong>
                              <span className="detail">{check.detail}</span>
                            </span>
                          </li>
                        );
                      })}
                    </ul>
                  </div>
                </div>

                <div className="card clone-assembly-card">
                  <div className="card-head">
                    <div style={{ flex: 1 }}>
                      <h2>Insert–vector end compatibility</h2>
                      <span className="label">Vector + duplex insert → recombinant product</span>
                    </div>
                    <div className="seg-toggle" role="group" aria-label="Junction detail level">
                      {(["simple", "detailed"] as const).map((level) => (
                        <button
                          key={level}
                          type="button"
                          className={assemblyDetail === level ? "on" : ""}
                          aria-pressed={assemblyDetail === level}
                          onClick={() => setAssemblyDetail(level)}
                        >
                          {level[0].toUpperCase() + level.slice(1)}
                        </button>
                      ))}
                    </div>
                  </div>
                  <div className="card-body clone-assembly-body">
                    <div className="clone-molecule-flow" aria-label="Cloning assembly flow">
                      <div><span className="label">Vector</span><strong>{result.vector_name}</strong><small>{result.backbone_length.toLocaleString()} bp after digestion</small></div>
                      <span className="clone-flow-operator" aria-hidden="true">+</span>
                      <div><span className="label">Insert</span><strong>{params.name || "construct"}</strong><small>{result.insert_length.toLocaleString()} bp duplex</small></div>
                      <Icon name="arrowRight" size={20} />
                      <div className={ligationCommitted ? "product ready" : "product"}><span className="label">Product</span><strong>{ligationCommitted ? `${result.length.toLocaleString()} bp plasmid` : "Waiting for ligation"}</strong><small>{ligationCommitted ? "joined in silico" : "review both junctions first"}</small></div>
                    </div>

                    {assemblyDetail === "simple" ? (
                      <div className="clone-junction-grid">
                        {result.junction_views.map((view) => <CloneJunctionSummary key={view.name} view={view} />)}
                      </div>
                    ) : (
                      <div className="clone-junction-detail">
                        {result.junction_views.map((view) => (
                          <JunctionDuplex key={view.name} view={view} showEnds ligated={ligationCommitted} />
                        ))}
                      </div>
                    )}

                    <div className="clone-ligation-action">
                      <div>
                        <strong>{result.is_clonable ? "Every exposed insert end has the complementary vector end." : "At least one junction is incompatible."}</strong>
                        <span>{result.is_clonable ? "The button simulates joining the phosphodiester backbone; it does not represent an experimental ligation." : "Ligation remains disabled until enzyme sites, polarity and overhang sequences all agree."}</span>
                      </div>
                      <button
                        className="btn btn-primary"
                        onClick={() => setLigationCommitted(true)}
                        disabled={ligationCommitted || !result.is_clonable || result.preflight?.can_export === false}
                      >
                        <Icon name="plate" size={17} />
                        {ligationCommitted ? "Ends ligated" : "Ligate compatible ends"}
                      </button>
                    </div>
                  </div>
                </div>

                {ligationCommitted && result.is_clonable && result.preflight?.can_export !== false && (
                  <LigationOutcome
                    vectorName={result.vector_name}
                    backboneLength={result.backbone_length}
                    insertName={params.name || "construct"}
                    insertLength={result.insert_length}
                    productName={result.name}
                    productLength={result.length}
                    junctions={result.junction_views}
                  />
                )}

                <ConstructWorkbench
                  result={result}
                  vector={vector}
                  constructName={params.name || result.name}
                  catalogue={catalogue}
                  ligated={ligationCommitted}
                  annotations={productAnnotations ?? result.annotations}
                  onAnnotationsChange={setProductAnnotations}
                  onExport={(filetype) => void exportPlasmid(filetype)}
                  onWorksheet={() => void downloadWorksheet()}
                  onSave={() => void runClone(true)}
                  busy={busy}
                />

                {ligationCommitted && result.protein && (
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
                              <span className="pill">{tag.end}-term {tag.name}</span>
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

function CloneJunctionSummary({ view }: { view: JunctionView }) {
  const [start, end] = view.overhang_span;
  const vectorEnd = view.overhang
    ? view.joined_top.slice(start, end).trim()
    : "blunt";
  const insertEnd = view.overhang
    ? view.joined_bottom.slice(start, end).trim()
    : "blunt";
  const pairing = view.overhang
    ? view.joined_pairs.slice(start, end).trim()
    : "flush";

  return (
    <div className={view.compatible ? "clone-junction-summary compatible" : "clone-junction-summary incompatible"}>
      <div className="clone-junction-summary-head">
        <div>
          <span className="label">{view.name}</span>
          <strong>{view.enzyme}</strong>
        </div>
        <span className={view.compatible ? "pill pill-ok" : "pill pill-bad"}>
          {view.compatible ? "complementary" : "blocked"}
        </span>
      </div>
      <div className="clone-junction-seqs" aria-label={`${view.name} cohesive ends`}>
        <div>
          <span>Vector end</span>
          <code>{view.overhang ? `5′-${vectorEnd}-3′` : "blunt"}</code>
        </div>
        <span className="clone-pair-mark" aria-hidden="true">
          {view.compatible ? pairing : "×".repeat(Math.max(view.overhang.length, 1))}
        </span>
        <div>
          <span>Insert end</span>
          <code>{view.overhang ? `3′-${insertEnd}-5′` : "blunt"}</code>
        </div>
      </div>
      <p>{view.compatible ? `${view.kind} cohesive ends align with opposite polarity.` : view.reason}</p>
    </div>
  );
}
