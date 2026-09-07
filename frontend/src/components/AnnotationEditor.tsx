import { useEffect, useId, useRef, useState } from "react";

import type { Annotation } from "../api/client";

const TYPES = [
  ["misc_feature", "Feature / insert"],
  ["CDS", "Coding sequence (CDS)"],
  ["mat_peptide", "Mature peptide / insert"],
  ["promoter", "Promoter"],
  ["RBS", "Ribosome-binding site (mRNA)"],
  ["regulatory", "Regulatory region"],
  ["protein_bind", "Protein-binding site"],
  ["terminator", "Terminator"],
  ["rep_origin", "Replication origin"],
] as const;

const DEFAULT_COLOUR = "#3F7A52";

export type AnnotationDraft = {
  name: string;
  type: string;
  start: string;
  end: string;
  direction: string;
  color: string;
  wraps: boolean;
};

export function draftFromAnnotation(
  annotation: Annotation | null,
  sequenceLength: number,
): AnnotationDraft {
  if (!annotation) {
    return {
      name: "",
      type: "misc_feature",
      start: "1",
      end: String(sequenceLength),
      direction: "1",
      color: DEFAULT_COLOUR,
      wraps: false,
    };
  }
  return {
    name: annotation.name,
    type: annotation.type,
    start: String(annotation.start + 1),
    end: String(annotation.end > sequenceLength ? annotation.end - sequenceLength : annotation.end),
    direction: String(annotation.direction),
    color: annotation.color || DEFAULT_COLOUR,
    wraps: annotation.end > sequenceLength,
  };
}

export function annotationFromDraft(
  draft: AnnotationDraft,
  sequenceLength: number,
  circular: boolean,
  previous?: Annotation | null,
): { annotation?: Annotation; error?: string } {
  const name = draft.name.trim();
  const startOneBased = Number(draft.start);
  const endOneBased = Number(draft.end);
  if (!name) return { error: "Give this feature a name." };
  if (!Number.isInteger(startOneBased) || startOneBased < 1 || startOneBased > sequenceLength) {
    return { error: `Start must be a whole number from 1 to ${sequenceLength.toLocaleString()}.` };
  }
  if (!Number.isInteger(endOneBased) || endOneBased < 1 || endOneBased > sequenceLength) {
    return { error: `End must be a whole number from 1 to ${sequenceLength.toLocaleString()}.` };
  }
  if (draft.wraps && !circular) return { error: "Only a circular sequence can cross its origin." };
  if (draft.wraps && endOneBased >= startOneBased) {
    return { error: "For an origin-crossing feature, the end coordinate must be before the start." };
  }
  if (!draft.wraps && endOneBased < startOneBased) {
    return { error: "End must be at or after start, or select ‘Crosses origin’." };
  }

  const start = startOneBased - 1;
  const end = draft.wraps ? sequenceLength + endOneBased : endOneBased;
  const annotation: Annotation = {
    name,
    type: draft.type,
    start,
    end,
    direction: Number(draft.direction),
    color: draft.color.toUpperCase(),
  };
  if (previous?.truncated) annotation.truncated = true;
  const unchangedSpan = previous?.start === start && previous.end === end
    && previous.direction === annotation.direction && previous.type === draft.type;
  if (previous?.inferred) {
    annotation.inferred = true;
    annotation.basis = unchangedSpan ? previous.basis : "Annotation edited; sequence evidence and biological function require review.";
  } else if (unchangedSpan && previous?.basis) annotation.basis = previous.basis;
  if (previous?.type === draft.type && previous.regulatory_class) annotation.regulatory_class = previous.regulatory_class;
  if (draft.type === "RBS") annotation.regulatory_class = "ribosome_binding_site";
  if (draft.type === "CDS") {
    annotation.translation_start = unchangedSpan ? previous?.translation_start ?? start : start;
    annotation.translation_end = unchangedSpan ? previous?.translation_end ?? end : end;
  }
  return { annotation };
}

type Props = {
  open: boolean;
  annotation: Annotation | null;
  initialAnnotation?: Annotation | null;
  sequenceLength: number;
  circular: boolean;
  saving: boolean;
  onSave: (annotation: Annotation) => void;
  onCancel: () => void;
};

export default function AnnotationEditor({
  open,
  annotation,
  initialAnnotation,
  sequenceLength,
  circular,
  saving,
  onSave,
  onCancel,
}: Props) {
  const [draft, setDraft] = useState<AnnotationDraft>(() =>
    draftFromAnnotation(annotation, sequenceLength));
  const [error, setError] = useState("");
  const titleId = useId();
  const firstField = useRef<HTMLInputElement>(null);

  useEffect(() => {
    if (!open) return;
    setDraft(draftFromAnnotation(annotation ?? initialAnnotation ?? null, sequenceLength));
    setError("");
    const opener = document.activeElement as HTMLElement | null;
    const timer = window.setTimeout(() => firstField.current?.focus(), 0);
    return () => { window.clearTimeout(timer); if (opener?.isConnected) opener.focus(); };
  }, [annotation, initialAnnotation, open, sequenceLength]);

  if (!open) return null;

  const set = <K extends keyof AnnotationDraft>(key: K, value: AnnotationDraft[K]) =>
    setDraft((current) => ({ ...current, [key]: value }));

  return (
    <div className="modal-backdrop" onMouseDown={(event) => {
      if (event.target === event.currentTarget && !saving) onCancel();
    }}>
      <form
        className="modal annotation-editor"
        role="dialog"
        aria-modal="true"
        aria-labelledby={titleId}
        onKeyDown={(event) => {
          if (event.key === "Escape") { event.stopPropagation(); if (!saving) onCancel(); }
          if (event.key !== "Tab") return;
          event.stopPropagation();
          const controls = [...event.currentTarget.querySelectorAll<HTMLElement>('button:not([disabled]), input:not([disabled]), select:not([disabled])')];
          const first = controls[0]; const last = controls[controls.length - 1];
          if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last?.focus(); }
          if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first?.focus(); }
        }}
        onSubmit={(event) => {
          event.preventDefault();
          const result = annotationFromDraft(draft, sequenceLength, circular, annotation);
          if (!result.annotation) {
            setError(result.error ?? "Check the annotation fields.");
            return;
          }
          setError("");
          onSave(result.annotation);
        }}
      >
        <h2 id={titleId}>{annotation ? "Edit feature" : "Annotate a feature"}</h2>
        <p>
          Name a new insert or describe any sequence span. Coordinates are 1-based and inclusive,
          as they appear on the sequence ruler.
        </p>
        {error && <div className="notice notice-error" role="alert">{error}</div>}
        <div className="annotation-editor-grid">
          <div className="field annotation-name-field">
            <label htmlFor="annotation-name">Feature name</label>
            <input
              id="annotation-name"
              ref={firstField}
              type="text"
              maxLength={200}
              value={draft.name}
              placeholder="e.g. Insulin glargine B-chain"
              onChange={(event) => set("name", event.target.value)}
            />
          </div>
          <div className="field">
            <label htmlFor="annotation-type">Type</label>
            <select
              id="annotation-type"
              value={draft.type}
              onChange={(event) => set("type", event.target.value)}
            >
              {!TYPES.some(([type]) => type === draft.type) && <option value={draft.type}>{draft.type}</option>}
              {TYPES.map(([value, label]) => <option key={value} value={value}>{label}</option>)}
            </select>
          </div>
          <div className="field">
            <label htmlFor="annotation-start">Start base</label>
            <input
              id="annotation-start"
              type="number"
              min={1}
              max={sequenceLength}
              value={draft.start}
              onChange={(event) => set("start", event.target.value)}
            />
          </div>
          <div className="field">
            <label htmlFor="annotation-end">End base</label>
            <input
              id="annotation-end"
              type="number"
              min={1}
              max={sequenceLength}
              value={draft.end}
              onChange={(event) => set("end", event.target.value)}
            />
          </div>
          <div className="field">
            <label htmlFor="annotation-strand">Strand</label>
            <select
              id="annotation-strand"
              value={draft.direction}
              onChange={(event) => set("direction", event.target.value)}
            >
              <option value="1">Forward (+)</option>
              <option value="-1">Reverse (−)</option>
              <option value="0">Unstranded</option>
            </select>
          </div>
          <div className="field">
            <label htmlFor="annotation-colour">Colour</label>
            <input
              id="annotation-colour"
              className="annotation-colour-input"
              type="color"
              value={draft.color}
              onChange={(event) => set("color", event.target.value)}
            />
          </div>
        </div>
        {circular && (
          <label className="checkbox-row annotation-wrap-check">
            <input
              type="checkbox"
              checked={draft.wraps}
              onChange={(event) => set("wraps", event.target.checked)}
            />
            <span><strong>Crosses origin</strong><small>Use when the feature starts near the end and finishes near base 1.</small></span>
          </label>
        )}
        <div className="modal-actions">
          <button type="button" className="btn btn-outline" onClick={onCancel} disabled={saving}>
            Cancel
          </button>
          <button type="submit" className="btn btn-primary" disabled={saving}>
            {saving ? "Saving…" : annotation ? "Save changes" : "Add feature"}
          </button>
        </div>
      </form>
    </div>
  );
}
