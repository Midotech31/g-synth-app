import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import AnnotationEditor, { annotationFromDraft } from "../AnnotationEditor";

describe("editable sequence annotations", () => {
  it("prefills a selected range while retaining the add-annotation operation", () => {
    const onSave = vi.fn();
    render(<AnnotationEditor open annotation={null}
      initialAnnotation={{ name: "", type: "misc_feature", start: 58, end: 63, direction: 1, color: "#3F7A52" }}
      sequenceLength={160} circular={false} saving={false} onCancel={vi.fn()} onSave={onSave} />);
    expect(screen.getByLabelText("Start base")).toHaveValue(59);
    expect(screen.getByLabelText("End base")).toHaveValue(63);
    fireEvent.change(screen.getByLabelText("Feature name"), { target: { value: "Selected target" } });
    fireEvent.click(screen.getByRole("button", { name: "Add feature" }));
    expect(onSave).toHaveBeenCalledWith(expect.objectContaining({ name: "Selected target", start: 58, end: 63 }));
  });

  it("converts displayed 1-based inclusive coordinates to engine coordinates", () => {
    const result = annotationFromDraft({
      name: "New insulin insert", type: "CDS", start: "11", end: "40",
      direction: "1", color: "#3f7a52", wraps: false,
    }, 100, false);

    expect(result.annotation).toMatchObject({
      name: "New insulin insert", start: 10, end: 40, direction: 1,
      translation_start: 10, translation_end: 40, color: "#3F7A52",
    });
  });

  it("represents a circular origin-crossing feature with an extended end", () => {
    const result = annotationFromDraft({
      name: "wrapped insert", type: "misc_feature", start: "91", end: "12",
      direction: "-1", color: "#0e6e77", wraps: true,
    }, 100, true);

    expect(result.annotation).toMatchObject({ start: 90, end: 112, direction: -1 });
  });

  it("requires an explicit origin crossing when end precedes start", () => {
    const result = annotationFromDraft({
      name: "bad span", type: "misc_feature", start: "91", end: "12",
      direction: "1", color: "#0e6e77", wraps: false,
    }, 100, true);

    expect(result.error).toMatch(/Crosses origin/i);
  });

  it("lets a scientist name and add an unannotated insert", () => {
    const onSave = vi.fn();
    render(
      <AnnotationEditor
        open
        annotation={null}
        sequenceLength={120}
        circular={false}
        saving={false}
        onSave={onSave}
        onCancel={vi.fn()}
      />,
    );

    fireEvent.change(screen.getByLabelText("Feature name"), {
      target: { value: "Insulin glargine A-chain" },
    });
    fireEvent.change(screen.getByLabelText("End base"), { target: { value: "63" } });
    fireEvent.click(screen.getByRole("button", { name: "Add feature" }));

    expect(onSave).toHaveBeenCalledWith(expect.objectContaining({
      name: "Insulin glargine A-chain", start: 0, end: 63,
    }));
  });
});

it("renaming a CDS preserves its translation origin", () => {
  const previous = { name: "CDS", type: "CDS", start: 10, end: 31, translation_start: 11,
    translation_end: 29, direction: 1, color: "#0E6E77", inferred: true, basis: "Sequence evidence." };
  const draft = { name: "Renamed CDS", type: "CDS", start: "11", end: "31", direction: "1", color: "#0E6E77", wraps: false };
  const result = annotationFromDraft(draft, 100, false, previous);
  expect(result.annotation?.translation_start).toBe(11);
  expect(result.annotation?.translation_end).toBe(29);
  expect(result.annotation?.inferred).toBe(true);
  expect(result.annotation?.basis).toBe("Sequence evidence.");
});
