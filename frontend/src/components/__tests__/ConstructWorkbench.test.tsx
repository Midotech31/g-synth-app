import { fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import type { Annotation, Catalogue, CloneResult } from "../../api/client";
import ConstructWorkbench, { restrictionSitesForSequence } from "../ConstructWorkbench";

vi.mock("seqviz", () => ({
  SeqViz: ({ onSelection }: { onSelection?: (selection: { type: string; start: number; end: number }) => void }) => (
    <>
    <button type="button" onClick={() => onSelection?.({ type: "ANNOTATION", start: 4, end: 8 })}>
      Select map feature
    </button>
    <button type="button" onClick={() => onSelection?.({ type: "ANNOTATION", start: 0, end: 3 })}>
      Select site across origin
    </button>
    </>
  ),
}));

const annotations: Annotation[] = [{
  name: "Therapeutic insert",
  type: "CDS",
  start: 4,
  end: 16,
  direction: 1,
  color: "#0E6E77",
}];

const catalogue: Catalogue = {
  enzymes: [
    {
      name: "HindIII",
      recognition: "AAGCTT",
      overhang: "AGCT",
      overhang_type: "5'",
      supplies_start_codon: false,
      common: true,
    },
    {
      name: "AccBSI",
      recognition: "CCGCTC",
      overhang: "",
      overhang_type: "blunt",
      supplies_start_codon: false,
      common: false,
    },
  ],
  common_pairs: [],
  cleavage_sites: [],
};

const result = {
  plasmid: "AAAATGCCAAGCTTGGGGTTTT",
  name: "pGS-test",
  vector_name: "pTest",
  length: 22,
  gc: 45,
  topology: "circular",
  insert_start: 4,
  insert_end: 16,
  insert_length: 12,
  backbone_length: 10,
  removed_length: 0,
  left_enzyme: "NdeI",
  right_enzyme: "HindIII",
  protein: "MP",
  protein_length: 2,
  reversed_insert: false,
  tags: [],
  vector: { recognised: true, spec: null, check: null },
  annotations,
  junctions: [],
  orfs: [],
  junction_views: [],
  restriction_sites: [{
    name: "HindIII", type: "restriction_site", start: 8, end: 14,
    direction: 1, color: "#9E3D3D", recognition: "AAGCTT", cuts: 1,
    used: true, wraps: false,
  }],
  validation: [],
  warnings: [],
  problems: [],
  is_clonable: true,
  insert: null,
  assembly: null,
  preflight: { workflow: "cloning", verdict: "ready", can_export: true, checks: [], diagnostics: [] },
} as unknown as CloneResult;

const vector = {
  name: "pTest",
  sequence: "AAAATGCCAAGCTTGGGGTTTT",
  annotations,
  circular: true,
};

function renderWorkbench(ligated: boolean, onAnnotationsChange = vi.fn()) {
  return render(
    <ConstructWorkbench
      result={result}
      vector={vector}
      constructName="Therapeutic insert"
      catalogue={catalogue}
      ligated={ligated}
      annotations={annotations}
      onAnnotationsChange={onAnnotationsChange}
      onExport={vi.fn()}
      onWorksheet={vi.fn()}
      onSave={vi.fn()}
      busy={false}
    />,
  );
}

describe("construct workbench", () => {
  it("inspects the entire restriction site when its map segment crosses the origin", async () => {
    const plasmid = "CTTGGGGGGGGGGAAG";
    render(<ConstructWorkbench result={{ ...result, plasmid, length: plasmid.length,
      restriction_sites: restrictionSitesForSequence(plasmid, true, catalogue, new Set(["HindIII"])) }}
      vector={vector} constructName="Origin example" catalogue={catalogue} ligated
      annotations={[]} onAnnotationsChange={vi.fn()} onExport={vi.fn()}
      onWorksheet={vi.fn()} onSave={vi.fn()} busy={false} />);
    fireEvent.click(screen.getByRole("button", { name: "Select site across origin" }));
    const inspector = within(screen.getByRole("complementary", { name: "Selection inspector" }));
    expect(inspector.getByText("14–3 (across origin)")).toBeInTheDocument();
    expect(inspector.getAllByText("AAGCTT").length).toBeGreaterThan(0);
    fireEvent.click(inspector.getByRole("button", { name: "View sequence" }));
    expect(screen.getByRole("region", { name: "Coordinate-level annotated sequence" })).toBeInTheDocument();
  });

  it("finds HindIII and preserves its exact recognition coordinates", () => {
    const sites = restrictionSitesForSequence("CCAAGCTTGG", false, catalogue, new Set(["HindIII"]));
    expect(sites).toHaveLength(1);
    expect(sites[0]).toMatchObject({
      name: "HindIII", recognition: "AAGCTT", start: 2, end: 8, cuts: 1, used: true,
    });
  });

  it("detects non-palindromic enzyme sites on the reverse strand", () => {
    const sites = restrictionSitesForSequence("TTTGAGCGGAAA", false, catalogue, new Set());
    expect(sites).toEqual(expect.arrayContaining([
      expect.objectContaining({ name: "AccBSI", recognition: "CCGCTC", start: 3, end: 9 }),
    ]));
  });

  it("keeps the product and downstream evidence locked until ligation", () => {
    renderWorkbench(false);
    expect(screen.getByRole("tab", { name: /Product.*Awaiting ligation/ })).toHaveAttribute("aria-disabled", "true");
    expect(screen.getByRole("tab", { name: "Gel" })).toHaveAttribute("aria-disabled", "true");
    expect(screen.getByRole("button", { name: "Save plasmid" })).toBeDisabled();
    expect(screen.getByText("Product intentionally unavailable")).toBeInTheDocument();
  });

  it("shares a map selection with the inspector and permits reviewed product annotations", async () => {
    const onAnnotationsChange = vi.fn();
    renderWorkbench(true, onAnnotationsChange);

    await waitFor(() => expect(screen.getByRole("tab", { name: /Product.*pGS-test/ })).toHaveAttribute("aria-selected", "true"));
    fireEvent.click(screen.getByRole("button", { name: "Select map feature" }));
    expect(screen.getAllByText("Therapeutic insert").length).toBeGreaterThan(0);
    expect(screen.getByText("5–16")).toBeInTheDocument();

    fireEvent.click(screen.getByRole("tab", { name: /Annotations.*1/ }));
    fireEvent.click(screen.getByRole("button", { name: "Add feature" }));
    const dialog = within(screen.getByRole("dialog"));
    fireEvent.change(dialog.getByLabelText("Feature name"), { target: { value: "Sequencing primer span" } });
    fireEvent.change(dialog.getByLabelText("Start base"), { target: { value: "2" } });
    fireEvent.change(dialog.getByLabelText("End base"), { target: { value: "6" } });
    fireEvent.click(dialog.getByRole("button", { name: "Add feature" }));

    expect(onAnnotationsChange).toHaveBeenCalledWith(expect.arrayContaining([
      expect.objectContaining({ name: "Sequencing primer span", start: 1, end: 6 }),
    ]));
  });
});
