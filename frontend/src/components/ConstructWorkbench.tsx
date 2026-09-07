import { useEffect, useMemo, useState } from "react";
import { SeqViz } from "seqviz";

import type {
  Annotation,
  Catalogue,
  CloneResult,
  RestrictionSite,
} from "../api/client";
import AnnotatedSequenceView from "./AnnotatedSequenceView";
import ExpandablePanel from "./ExpandablePanel";
import { useDetectedFeatures } from "../hooks/useDetectedFeatures";
import FeatureEvidence from "./FeatureEvidence";
import AnnotationEditor from "./AnnotationEditor";
import GelSimulation from "./GelSimulation";
import Icon from "./Icon";

type MoleculeStage = "vector" | "insert" | "product";
type InspectionView = "map" | "sequence" | "annotations" | "enzymes" | "oligos" | "gel";
type MapView = "circular" | "linear" | "both";

type VectorSnapshot = {
  name: string;
  sequence: string;
  annotations: Annotation[];
  circular: boolean;
};

type SelectedItem = {
  kind: "feature" | "site";
  name: string;
  type?: string;
  start: number;
  end: number;
  direction: number;
  color: string;
  annotationIndex?: number;
  inferred?: boolean;
  basis?: string;
  regulatory_class?: string;
  recognition?: string;
  cuts?: number;
  used?: boolean;
  wraps?: boolean;
};

type Props = {
  result: CloneResult;
  vector: VectorSnapshot;
  constructName: string;
  catalogue: Catalogue | null;
  ligated: boolean;
  annotations: Annotation[];
  onAnnotationsChange: (annotations: Annotation[]) => void;
  onExport: (filetype: "genbank" | "fasta" | "sbol3") => void;
  onWorksheet: () => void;
  onSave: () => void;
  busy: boolean;
};

const COMPLEMENT: Record<string, string> = {
  A: "T", T: "A", G: "C", C: "G", R: "Y", Y: "R", S: "S",
  W: "W", K: "M", M: "K", B: "V", D: "H", H: "D", V: "B", N: "N",
};

const IUPAC: Record<string, string> = {
  A: "A", C: "C", G: "G", T: "T", R: "[AG]", Y: "[CT]", S: "[GC]",
  W: "[AT]", K: "[GT]", M: "[AC]", B: "[CGT]", D: "[AGT]",
  H: "[ACT]", V: "[ACG]", N: "[ACGT]",
};

function reverseComplement(sequence: string): string {
  return sequence
    .toUpperCase()
    .split("")
    .reverse()
    .map((base) => COMPLEMENT[base] ?? "N")
    .join("");
}

function circularSpan(sequence: string, start: number, end: number): string {
  if (!sequence.length || end <= start) return "";
  return Array.from({ length: end - start }, (_, index) =>
    sequence[(start + index) % sequence.length]).join("");
}

function selectedSequence(sequence: string, item: SelectedItem): string {
  const span = circularSpan(sequence, item.start, item.end);
  return item.direction === -1 ? reverseComplement(span) : span;
}

export function restrictionSitesForSequence(
  sequence: string,
  circular: boolean,
  catalogue: Catalogue | null,
  usedNames: Set<string>,
): RestrictionSite[] {
  if (!catalogue || !sequence.length) return [];
  const upper = sequence.toUpperCase();
  const sites: RestrictionSite[] = [];

  for (const enzyme of catalogue.enzymes) {
    const motif = enzyme.recognition.toUpperCase();
    if (!motif || motif.length > upper.length) continue;
    const searchable = circular ? upper + upper.slice(0, motif.length - 1) : upper;
    const positions = new Set<number>();
    for (const strandMotif of new Set([motif, reverseComplement(motif)])) {
      const pattern = strandMotif.split("").map((base) => IUPAC[base] ?? base).join("");
      for (const match of searchable.matchAll(new RegExp(`(?=${pattern})`, "g"))) {
        const position = match.index ?? -1;
        if (position >= 0 && position < upper.length) positions.add(position);
      }
    }
    const unique = [...positions].sort((left, right) => left - right);
    for (const position of unique) {
      sites.push({
        name: enzyme.name,
        type: "restriction_site",
        start: position,
        end: position + motif.length,
        direction: 1,
        color: usedNames.has(enzyme.name) ? "#9E3D3D" : "#C97634",
        recognition: enzyme.recognition,
        cuts: unique.length,
        used: usedNames.has(enzyme.name),
        wraps: position + motif.length > upper.length,
      });
    }
  }
  return sites;
}

function featureForStage(
  annotation: Annotation,
  index: number,
  insertStart: number,
  insertEnd: number,
): (Annotation & { sourceIndex: number }) | null {
  if (annotation.end <= insertStart || annotation.start >= insertEnd) return null;
  return {
    ...annotation,
    start: Math.max(annotation.start, insertStart) - insertStart,
    end: Math.min(annotation.end, insertEnd) - insertStart,
    translation_start: annotation.translation_start === undefined
      ? undefined
      : Math.max(annotation.translation_start, insertStart) - insertStart,
    translation_end: annotation.translation_end === undefined
      ? undefined
      : Math.min(annotation.translation_end, insertEnd) - insertStart,
    sourceIndex: index,
  };
}

function toSelectedFeature(annotation: Annotation, annotationIndex?: number): SelectedItem {
  return {
    kind: "feature",
    name: annotation.name,
    type: annotation.type,
    start: annotation.start,
    end: annotation.end,
    direction: annotation.direction,
    color: annotation.color,
    annotationIndex,
    inferred: annotation.inferred,
    basis: annotation.basis,
    regulatory_class: annotation.regulatory_class,
    wraps: false,
  };
}

function toSelectedSite(site: RestrictionSite): SelectedItem {
  return {
    kind: "site",
    name: site.name,
    type: site.type,
    start: site.start,
    end: site.end,
    direction: 1,
    color: site.color,
    recognition: site.recognition,
    cuts: site.cuts,
    used: site.used,
    wraps: site.wraps,
  };
}

export default function ConstructWorkbench({
  result,
  vector,
  constructName,
  catalogue,
  ligated,
  annotations,
  onAnnotationsChange,
  onExport,
  onWorksheet,
  onSave,
  busy,
}: Props) {
  const [stage, setStage] = useState<MoleculeStage>("insert");
  const [view, setView] = useState<InspectionView>("map");
  const [mapView, setMapView] = useState<MapView>("linear");
  const [selected, setSelected] = useState<SelectedItem | null>(null);
  const [showSites, setShowSites] = useState(true);
  const [onlyUsedSites, setOnlyUsedSites] = useState(false);
  const [showMultiCutters, setShowMultiCutters] = useState(false);
  const [enzymeQuery, setEnzymeQuery] = useState("");
  const [editorOpen, setEditorOpen] = useState(false);
  const [editingIndex, setEditingIndex] = useState<number | null>(null);
  const [initialAnnotation, setInitialAnnotation] = useState<Annotation | null>(null);
  const [annotationHistory, setAnnotationHistory] = useState<Annotation[][]>([]);

  useEffect(() => {
    if (!ligated) return;
    setStage("product");
    setView("map");
    setMapView("circular");
    setSelected(null);
  }, [ligated]);

  useEffect(() => {
    if (!ligated && stage === "product") setStage("insert");
    if (!ligated && view === "gel") setView("map");
  }, [ligated, stage, view]);

  const insertSequence = useMemo(
    () => circularSpan(result.plasmid, result.insert_start, result.insert_end),
    [result.insert_end, result.insert_start, result.plasmid],
  );

  const insertAnnotations = useMemo(
    () => annotations
      .map((annotation, index) => featureForStage(
        annotation, index, result.insert_start, result.insert_end,
      ))
      .filter((annotation): annotation is Annotation & { sourceIndex: number } => annotation !== null),
    [annotations, result.insert_end, result.insert_start],
  );

  const stageSequence = stage === "vector"
    ? vector.sequence
    : stage === "insert"
      ? insertSequence
      : result.plasmid;
  const baseAnnotations = stage === "vector"
    ? vector.annotations
    : stage === "insert"
      ? insertAnnotations
      : annotations;
  const stageCircular = stage === "vector" ? vector.circular : stage === "product";
  const scan = useDetectedFeatures(stageSequence, baseAnnotations, stageCircular);
  const [showDetected, setShowDetected] = useState(true);
  const [showInspector, setShowInspector] = useState(true);
  const stageAnnotations = useMemo(() => [...baseAnnotations,
    ...(showDetected ? scan.matches.map((match) => match.annotation) : [])],
  [baseAnnotations, showDetected, scan.matches]);
  const stageName = stage === "vector"
    ? vector.name || result.vector_name
    : stage === "insert"
      ? constructName || result.name
      : result.name;

  const usedNames = useMemo(
    () => new Set([result.left_enzyme, result.right_enzyme]),
    [result.left_enzyme, result.right_enzyme],
  );
  const stageSites = useMemo(() => {
    if (stage === "product") return result.restriction_sites;
    return restrictionSitesForSequence(stageSequence, stageCircular, catalogue, usedNames);
  }, [catalogue, result.restriction_sites, stage, stageCircular, stageSequence, usedNames]);
  const visibleSites = useMemo(
    () => stageSites.filter((site) => {
      if (onlyUsedSites) return site.used;
      return showMultiCutters || site.cuts === 1 || site.used;
    }),
    [onlyUsedSites, showMultiCutters, stageSites],
  );
  const stageEnzymeCount = useMemo(
    () => new Set(stageSites.map((site) => site.name)).size,
    [stageSites],
  );
  const visibleEnzymeSites = useMemo(() => {
    const query = enzymeQuery.trim().toLowerCase();
    return visibleSites.filter((site) => !query
      || site.name.toLowerCase().includes(query)
      || site.recognition.toLowerCase().includes(query));
  }, [enzymeQuery, visibleSites]);

  const mapSiteParts = useMemo(() => visibleSites.flatMap((site) => site.wraps
    ? [
        { ...site, end: stageSequence.length },
        { ...site, start: 0, end: site.end - stageSequence.length },
      ]
    : [site]), [stageSequence.length, visibleSites]);

  const mapAnnotations = useMemo(() => {
    const features = stageAnnotations.map((annotation) => ({
      name: annotation.name,
      start: annotation.start,
      end: annotation.end,
      direction: annotation.direction as 1 | -1 | 0,
      color: annotation.color,
    }));
    if (!showSites) return features;
    return [...features, ...mapSiteParts.map((site) => ({
      name: site.name,
      start: site.start,
      end: site.end,
      direction: 1 as const,
      color: site.color,
    }))];
  }, [mapSiteParts, showSites, stageAnnotations]);

  const clickable = useMemo(() => {
    const features = stageAnnotations.map((annotation, index) => toSelectedFeature(
      annotation,
      stage === "product"
        ? index
        : stage === "insert"
          ? insertAnnotations[index]?.sourceIndex
          : undefined,
    ));
    const sites = showSites ? visibleSites.map(toSelectedSite) : [];
    return [...features, ...sites];
  }, [insertAnnotations, visibleSites, showSites, stage, stageAnnotations]);

  const oligos = result.assembly?.oligos ?? [];
  const canExport = ligated && result.is_clonable && result.preflight?.can_export !== false;

  function chooseStage(next: MoleculeStage) {
    if (next === "product" && !ligated) return;
    setStage(next);
    setSelected(null);
    if (next === "insert") setMapView("linear");
    else if (mapView === "linear") setMapView("circular");
  }

  function chooseView(next: InspectionView) {
    if (next === "gel" && !ligated) return;
    setView(next);
  }

  function updateAnnotations(next: Annotation[]) {
    setAnnotationHistory((current) => [...current.slice(-19), annotations]);
    onAnnotationsChange(next);
  }

  function openNewAnnotation() {
    setInitialAnnotation(null);
    setEditingIndex(null);
    setEditorOpen(true);
  }

  function editSelectedAnnotation() {
    if (stage !== "product" || selected?.kind !== "feature" || selected.annotationIndex === undefined) return;
    setEditingIndex(selected.annotationIndex);
    setEditorOpen(true);
  }

  function removeSelectedAnnotation() {
    if (stage !== "product" || selected?.kind !== "feature" || selected.annotationIndex === undefined) return;
    updateAnnotations(annotations.filter((_, index) => index !== selected.annotationIndex));
    setSelected(null);
  }

  const productLockedMessage = "Ligate the compatible ends before inspecting or exporting a recombinant product.";

  return (
    <ExpandablePanel className="card construct-workbench" label="Construct workbench">
      <div className="workbench-head">
        <div>
          <span className="eyebrow">Construct workbench</span>
          <h2>Inspect one molecule without losing its context</h2>
          <p>Map, sequence, features, enzymes and oligos share one selection and coordinate system.</p>
        </div>
        <div className="workbench-actions" aria-label="Construct actions">
          <button type="button" className="btn btn-outline" onClick={() => onExport("genbank")} disabled={!canExport}>
            GenBank
          </button>
          <button type="button" className="btn btn-outline" onClick={() => onExport("fasta")} disabled={!canExport}>
            FASTA
          </button>
          <button type="button" className="btn btn-outline" onClick={() => onExport("sbol3")} disabled={!canExport}>
            SBOL 3
          </button>
          <button type="button" className="btn btn-outline" onClick={onWorksheet} disabled={!canExport}>
            Bench worksheet
          </button>
          <button type="button" className="btn btn-primary" onClick={onSave} disabled={!canExport || busy}>
            {busy ? "Saving…" : "Save plasmid"}
          </button>
        </div>
      </div>

      <div className="workbench-molecules" role="tablist" aria-label="Molecule">
        <button type="button" role="tab" aria-selected={stage === "vector"} className={stage === "vector" ? "active" : ""} onClick={() => chooseStage("vector")}>
          <span>Vector</span><strong>{result.vector_name}</strong><small>{vector.sequence.length.toLocaleString()} bp source</small>
        </button>
        <button type="button" role="tab" aria-selected={stage === "insert"} className={stage === "insert" ? "active" : ""} onClick={() => chooseStage("insert")}>
          <span>Insert</span><strong>{constructName || result.name}</strong><small>{result.insert_length.toLocaleString()} bp duplex</small>
        </button>
        <button
          role="tab"
          type="button"
          aria-selected={stage === "product"}
          aria-disabled={!ligated}
          className={`${stage === "product" ? "active" : ""}${!ligated ? " locked" : ""}`}
          onClick={() => chooseStage("product")}
          title={!ligated ? productLockedMessage : "Inspect the ligated recombinant product"}
        >
          <span>Product</span><strong>{ligated ? result.name : "Awaiting ligation"}</strong><small>{ligated ? `${result.length.toLocaleString()} bp circular` : "not created yet"}</small>
        </button>
      </div>

      <div className="workbench-view-tabs" role="tablist" aria-label="Inspection view">
        {([
          ["map", "Map"],
          ["sequence", "Sequence"],
          ["annotations", "Annotations"],
          ["enzymes", "Enzymes"],
          ["oligos", "Primers & oligos"],
          ["gel", "Gel"],
        ] as [InspectionView, string][]).map(([key, label]) => (
          <button
            key={key}
            role="tab"
            type="button"
            aria-selected={view === key}
            aria-disabled={key === "gel" && !ligated}
            className={view === key ? "active" : ""}
            onClick={() => chooseView(key)}
            title={key === "gel" && !ligated ? productLockedMessage : undefined}
          >
            {label}
            {key === "annotations" && <small>{stageAnnotations.length}</small>}
            {key === "enzymes" && <small>{stageEnzymeCount}</small>}
            {key === "oligos" && <small>{oligos.length}</small>}
          </button>
        ))}
      </div>

      <div className="feature-detection-toolbar">
        <label><input type="checkbox" checked={showDetected} onChange={(event) => setShowDetected(event.target.checked)} /> Show detected candidates</label>
        <span role="status">{scan.scanning ? "Scanning both strands…" : scan.error || `${scan.matches.length} additional candidates · review before saving`}</span>
        <button type="button" className="btn btn-ghost" onClick={scan.rescan} disabled={scan.scanning}>Rescan</button>
        <button type="button" className="btn btn-outline" aria-pressed={showInspector} onClick={() => setShowInspector(!showInspector)}>{showInspector ? "Hide details" : "Show details"}</button>
      </div>
      <div className={`workbench-layout${showInspector ? "" : " inspector-hidden"}`} data-view={view} role="tabpanel" aria-label={`${stageName} ${view} view`}>
        <div className="workbench-canvas">
          <div className="workbench-canvas-head">
            <div>
              <strong>{stageName}</strong>
              <span>{stageSequence.length.toLocaleString()} bp · {stageCircular ? "circular" : "linear"} · {stageAnnotations.length} features</span>
            </div>
            {view === "map" && (
              <div className="seg-toggle" role="group" aria-label="Map layout">
                {(["circular", "linear", "both"] as MapView[]).map((option) => (
                  <button
                    key={option}
                    type="button"
                    className={mapView === option ? "on" : ""}
                    aria-pressed={mapView === option}
                    disabled={!stageCircular && option !== "linear"}
                    onClick={() => setMapView(option)}
                  >
                    {option[0].toUpperCase() + option.slice(1)}
                  </button>
                ))}
              </div>
            )}
          </div>

          {view === "map" && (
            <>
              <div className="workbench-map-filters">
                <label><input type="checkbox" checked={showSites} onChange={(event) => setShowSites(event.target.checked)} /> Restriction sites</label>
                <label><input type="checkbox" checked={onlyUsedSites} disabled={!showSites} onChange={(event) => setOnlyUsedSites(event.target.checked)} /> Cloning pair only</label>
                <label><input type="checkbox" checked={showMultiCutters} disabled={!showSites || onlyUsedSites} onChange={(event) => setShowMultiCutters(event.target.checked)} /> Include multi-cutters</label>
                <span>{showSites ? `${visibleSites.length} sites shown` : "feature map only"}</span>
              </div>
              <div className="workbench-map-stage">
                <SeqViz
                  name=""
                  seq={stageSequence}
                  annotations={mapAnnotations}
                  viewer={stageCircular ? mapView : "linear"}
                  showIndex
                  showComplement={mapView !== "circular" || !stageCircular}
                  disableExternalFonts
                  onSelection={(selection) => {
                    if (selection.type !== "ANNOTATION" || selection.start === undefined || selection.end === undefined) return;
                    const covering = clickable.filter((item) => [0, stageSequence.length].some((offset) =>
                      item.start <= selection.start! + offset && item.end >= selection.end! + offset));
                    if (!covering.length) return;
                    setSelected(covering.reduce((smallest, item) =>
                      item.end - item.start < smallest.end - smallest.start ? item : smallest));
                  }}
                  highlights={selected ? [{ start: selected.start, end: selected.end, color: selected.color }] : []}
                  style={{ height: "100%", width: "100%" }}
                />
              </div>
            </>
          )}

          {view === "sequence" && (
            <AnnotatedSequenceView
              key={stage}
              sequence={stageSequence}
              annotations={stageAnnotations}
              selected={selected?.kind === "feature"
                ? stageAnnotations.find((annotation) => annotation.start === selected.start && annotation.end === selected.end && annotation.name === selected.name) ?? null
                : selected ? { ...selected, type: "restriction_site" } : null}
              preferredName={stage === "product" ? result.name : stageName}
              circular={stageCircular}
              onAnnotateRange={stage === "vector" ? undefined : (range) => {
                const offset = stage === "insert" ? result.insert_start : 0;
                setInitialAnnotation({ name: "", type: "misc_feature", direction: 1, color: "#3F7A52",
                  start: range.start + offset, end: range.end + offset });
                setEditingIndex(null); setEditorOpen(true);
              }}
              onSelect={(annotation) => {
                const localIndex = stageAnnotations.indexOf(annotation);
                const sourceIndex = stage === "product"
                  ? localIndex
                  : stage === "insert"
                    ? insertAnnotations[localIndex]?.sourceIndex
                    : undefined;
                setSelected(toSelectedFeature(annotation, sourceIndex));
              }}
            />
          )}

          {view === "annotations" && (
            <div className="workbench-list-view">
              <div className="workbench-list-toolbar">
                <div>
                  <strong>{stageAnnotations.length} annotated features</strong>
                  <span>{stage === "product" && ligated ? "Product annotations are included in saved projects and GenBank exports." : "Switch to the ligated product to edit final construct annotations."}</span>
                </div>
                <div>
                  <button
                    type="button"
                    className="btn btn-outline"
                    disabled={!annotationHistory.length || stage !== "product"}
                    onClick={() => {
                      const previous = annotationHistory[annotationHistory.length - 1];
                      setAnnotationHistory((current) => current.slice(0, -1));
                      onAnnotationsChange(previous);
                      setSelected(null);
                    }}
                  >
                    Undo
                  </button>
                  <button type="button" className="btn btn-primary" disabled={stage !== "product" || !ligated} onClick={openNewAnnotation}>
                    Add feature
                  </button>
                </div>
              </div>
              <div className="workbench-feature-list">
                {stageAnnotations.map((annotation, index) => (
                  <button
                    type="button"
                    key={`${annotation.name}-${annotation.start}-${index}`}
                    className={selected?.kind === "feature" && selected.start === annotation.start && selected.end === annotation.end ? "active" : ""}
                    onClick={() => {
                      const sourceIndex = stage === "product"
                        ? index
                        : stage === "insert"
                          ? insertAnnotations[index]?.sourceIndex
                          : undefined;
                      setSelected(toSelectedFeature(annotation, sourceIndex));
                    }}
                  >
                    <i style={{ background: annotation.color }} />
                    <span><strong>{annotation.name}</strong><small>{annotation.type} · {annotation.direction === -1 ? "reverse" : annotation.direction === 1 ? "forward" : "unstranded"}{annotation.inferred ? " · detected" : ""}</small></span>
                    <code>{(annotation.start + 1).toLocaleString()}–{annotation.end.toLocaleString()}</code>
                  </button>
                ))}
              </div>
            </div>
          )}

          {view === "enzymes" && (
            <div className="workbench-list-view">
              <div className="workbench-list-toolbar">
                <div>
                  <strong>{visibleEnzymeSites.length.toLocaleString()} of {stageSites.length.toLocaleString()} restriction sites shown</strong>
                  <span>{stageEnzymeCount.toLocaleString()} enzymes detected. Unique cutters and the cloning pair are shown by default.</span>
                </div>
                <div className="workbench-enzyme-controls">
                  <label className="workbench-enzyme-search">
                    <span>Filter enzymes</span>
                    <input
                      type="search"
                      value={enzymeQuery}
                      placeholder="HindIII or AAGCTT"
                      onChange={(event) => setEnzymeQuery(event.target.value)}
                    />
                  </label>
                  <label><input type="checkbox" checked={onlyUsedSites} onChange={(event) => setOnlyUsedSites(event.target.checked)} /> Cloning pair only</label>
                  <label><input type="checkbox" checked={showMultiCutters} disabled={onlyUsedSites} onChange={(event) => setShowMultiCutters(event.target.checked)} /> Include multi-cutters</label>
                </div>
              </div>
              <div className="workbench-enzyme-table" aria-label="Restriction sites">
                <div className="head" aria-hidden="true"><span>Enzyme</span><span>Recognition</span><span>Position</span><span>Frequency</span><span>Use</span></div>
                {visibleEnzymeSites.map((site, index) => (
                  <button
                    type="button"
                    key={`${site.name}-${site.start}-${index}`}
                    aria-label={`${site.name}, ${site.recognition}, position ${(site.start + 1).toLocaleString()}, ${site.cuts === 1 ? "unique" : `${site.cuts} cuts`}, ${site.used ? "cloning enzyme" : "diagnostic enzyme"}`}
                    onClick={() => setSelected(toSelectedSite(site))}
                  >
                    <strong>{site.name}</strong><code>{site.recognition}</code><span>{(site.start + 1).toLocaleString()}{site.wraps ? " ↻" : ""}</span><span>{site.cuts === 1 ? "unique" : `${site.cuts} cuts`}</span><span className={site.used ? "pill pill-ok" : "pill"}>{site.used ? "cloning" : "diagnostic"}</span>
                  </button>
                ))}
                {!stageSites.length && <p className="empty-inline">No catalogued restriction sites were detected in this molecule.</p>}
                {stageSites.length > 0 && !visibleEnzymeSites.length && <p className="empty-inline">No restriction sites match the current enzyme filters.</p>}
              </div>
            </div>
          )}

          {view === "oligos" && (
            <div className="workbench-list-view">
              <div className="workbench-list-toolbar">
                <div>
                  <strong>{oligos.length ? `${oligos.length} synthesis oligos` : "No synthesis oligo set attached"}</strong>
                  <span>ESD oligos remain tied to the insert; PCR primer design is performed from the PCR workspace.</span>
                </div>
              </div>
              {oligos.length ? (
                <div className="table-scroll">
                  <table className="data">
                    <thead><tr><th>Name</th><th>Sequence (5′→3′)</th><th>Length</th><th>Tm</th></tr></thead>
                    <tbody>{oligos.map((oligo, index) => (
                      <tr key={`${String(oligo.Name ?? "oligo")}-${index}`}>
                        <td className="mono">{oligo.Name}</td>
                        <td className="mono seq-cell">{oligo["Sequence (5'->3')"]}</td>
                        <td className="num">{oligo["Length (nt)"]}</td>
                        <td className="num">{oligo["Tm (°C)"]}</td>
                      </tr>
                    ))}</tbody>
                  </table>
                </div>
              ) : <p className="empty-inline">A pre-digested insert has no ESD oligo plan in this cloning record.</p>}
            </div>
          )}

          {view === "gel" && result.gel && (
            <div className="workbench-gel"><GelSimulation simulation={result.gel} /></div>
          )}
        </div>

        <aside hidden={!showInspector} className="workbench-inspector" aria-label="Selection inspector">
          <div className="workbench-inspector-section">
            <span className="eyebrow">Current molecule</span>
            <strong>{stageName}</strong>
            <dl>
              <div><dt>Length</dt><dd>{stageSequence.length.toLocaleString()} bp</dd></div>
              <div><dt>Topology</dt><dd>{stageCircular ? "circular" : "linear"}</dd></div>
              <div><dt>Orientation</dt><dd>{stage === "insert" ? (result.reversed_insert ? "reverse in vector coordinates" : "forward") : "reference"}</dd></div>
              {stage === "insert" && <div><dt>Directionality</dt><dd>{result.left_enzyme === result.right_enzyme ? "orientation requires review" : `fixed by ${result.left_enzyme}/${result.right_enzyme}`}</dd></div>}
            </dl>
          </div>

          {selected ? (
            <div className="workbench-inspector-section selection">
              <div className="selection-title"><i style={{ background: selected.color }} /><span><strong>{selected.name}</strong><small>{selected.kind === "site" ? "restriction site" : selected.type ?? "feature"}</small></span><button type="button" className="btn btn-ghost" onClick={() => setSelected(null)} aria-label="Clear selection"><Icon name="cross" size={14} /></button></div>
              <dl>
                <div><dt>Coordinates</dt><dd>{(selected.start + 1).toLocaleString()}–{(stageCircular ? ((selected.end - 1) % stageSequence.length) + 1 : selected.end).toLocaleString()}{selected.end > stageSequence.length && stageCircular ? " (across origin)" : ""}</dd></div>
                <div><dt>Length</dt><dd>{(selected.end - selected.start).toLocaleString()} bp</dd></div>
                {selected.kind === "feature" ? <div><dt>Strand</dt><dd>{selected.direction === -1 ? "reverse" : selected.direction === 1 ? "forward" : "unstranded"}</dd></div> : <>
                  <div><dt>Recognition</dt><dd className="mono">{selected.recognition}</dd></div>
                  <div><dt>Frequency</dt><dd>{selected.cuts === 1 ? "unique" : `${selected.cuts} cuts`}</dd></div>
                </>}
              </dl>
              <FeatureEvidence annotation={selected} />
              <button type="button" className="btn btn-outline" onClick={() => setView("sequence")}>View sequence</button>
              {stage === "product" && selected.kind === "feature" && !baseAnnotations.some((item) => item.name === selected.name && item.start === selected.start && item.end === selected.end) && (
                <button type="button" className="btn btn-primary" onClick={() => {
                  const match = scan.matches.find((item) => item.annotation.name === selected.name && item.annotation.start === selected.start && item.annotation.end === selected.end);
                  if (!match) return;
                  updateAnnotations([...annotations, match.annotation]);
                  setSelected(toSelectedFeature(match.annotation, annotations.length));
                }}>Keep candidate annotation</button>
              )}
              <div className="workbench-selected-sequence">{selectedSequence(stageSequence, selected)}</div>
              {stage === "product" && selected.kind === "feature" && baseAnnotations.some((item) => item.name === selected.name && item.start === selected.start && item.end === selected.end) && (
                <div className="workbench-selection-actions">
                  <button type="button" className="btn btn-outline" onClick={editSelectedAnnotation}>Edit</button>
                  <button type="button" className="btn btn-danger" onClick={removeSelectedAnnotation}>Remove</button>
                </div>
              )}
            </div>
          ) : (
            <div className="workbench-inspector-empty">
              <Icon name="target" size={26} />
              <strong>Select biological evidence</strong>
              <span>Choose a feature or restriction site on the map, sequence, annotation list or enzyme table. The same coordinates remain selected across views.</span>
            </div>
          )}

          {!ligated && (
            <div className="workbench-lock-note">
              <Icon name="scales" size={20} />
              <span><strong>Product intentionally unavailable</strong>Review the two junctions and perform the explicit ligation step first.</span>
            </div>
          )}
        </aside>
      </div>

      <AnnotationEditor
        open={editorOpen}
        annotation={editingIndex === null ? null : annotations[editingIndex] ?? null}
        initialAnnotation={initialAnnotation}
        sequenceLength={result.plasmid.length}
        circular
        saving={false}
        onCancel={() => setEditorOpen(false)}
        onSave={(annotation) => {
          const next = [...annotations];
          const index = editingIndex === null ? next.length : editingIndex;
          if (editingIndex === null) next.push(annotation); else next[editingIndex] = annotation;
          updateAnnotations(next);
          setSelected(toSelectedFeature(stage === "insert"
            ? { ...annotation, start: annotation.start - result.insert_start, end: annotation.end - result.insert_start }
            : annotation, index));
          setEditorOpen(false);
        }}
      />
    </ExpandablePanel>
  );
}
