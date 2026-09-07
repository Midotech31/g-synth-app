import { act, fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { beforeEach, expect, it, vi } from "vitest";

import { api, type CloneResult, type VectorRecord, type VectorSpec } from "../api/client";
import { WorkspaceStateProvider, useWorkspaceState } from "../state/WorkspaceStateContext";
import Clone from "./Clone";

const spec: VectorSpec = {
  key: "pET-21a", name: "pET-21a(+)", length: 100, resistance: "Ampicillin", promoter: "T7lac",
  host: "E. coli", supplier: "Reference", summary: "Old enzyme-dependent catalogue guidance",
  unique_sites: [], recommended_pairs: ["NdeI / XhoI"], tags: [],
  notes: ["Cloning at NdeI uses the vector's ATG and removes the T7·Tag."],
  reference: "", has_sequence: true, supplies_translation_start: true,
  expression_capable: true, tag_summary: "C-terminal His-tag",
};
const vector = { key: spec.key, name: spec.name, sequence: "ACGT".repeat(25),
  annotations: [{ name: "old feature", start: 1, end: 5, type: "misc_feature", direction: 1 }], circular: true, bundled: true };

function StateEvidence() {
  const [annotations] = useWorkspaceState("clone.productAnnotations", null);
  const [committed] = useWorkspaceState("clone.ligationCommitted", false);
  return <output data-testid="derived-state">{JSON.stringify({ annotations, committed })}</output>;
}

function renderPage(extra: Record<string, unknown> = {}, initialState?: object) {
  sessionStorage.setItem("gsynth.workspace.anonymous", JSON.stringify({
    "clone.vectorLoaded": true, "clone.vector": vector, ...extra,
  }));
  vi.spyOn(api, "catalogue").mockResolvedValue({ enzymes: ["NdeI", "XhoI", "BamHI", "EcoRI"].map(name => ({
    name, recognition: "ACGT", overhang: "", overhang_type: "blunt", supplies_start_codon: name === "NdeI", common: true,
  })), common_pairs: [], cleavage_sites: [] });
  vi.spyOn(api, "vectors").mockResolvedValue({ default: spec.key, vectors: [spec, { ...spec, key: "other", name: "Other reference" }] });
  render(<MemoryRouter initialEntries={[{ pathname: "/clone", state: initialState }]}>
    <WorkspaceStateProvider><Clone /><StateEvidence /></WorkspaceStateProvider>
  </MemoryRouter>);
  return screen.findByText(/Catalogue reference/);
}

beforeEach(() => { sessionStorage.clear(); vi.restoreAllMocks(); });

it("shows selected enzymes, removes catalogue outcome claims, and clears all derived state", async () => {
  await renderPage({ "clone.productAnnotations": vector.annotations, "clone.ligationCommitted": true, "clone.saved": "Previously saved" });
  const configuration = screen.getByRole("region", { name: "Cloning configuration" });
  expect(within(configuration).getByText("NdeI / XhoI")).toBeInTheDocument();
  expect(screen.queryByText(spec.notes[0])).not.toBeInTheDocument();
  expect(screen.queryByText(spec.summary)).not.toBeInTheDocument();
  fireEvent.change(screen.getByLabelText("5' enzyme"), { target: { value: "BamHI" } });
  fireEvent.change(screen.getByLabelText("3' enzyme"), { target: { value: "EcoRI" } });
  expect(within(configuration).getByText("BamHI / EcoRI")).toBeInTheDocument();
  expect(within(configuration).queryByText("NdeI / XhoI")).not.toBeInTheDocument();
  expect(screen.queryByText("Previously saved")).not.toBeInTheDocument();
  expect(screen.getByTestId("derived-state")).toHaveTextContent('{"annotations":null,"committed":false}');
});

it("uses transferred end assignments in both the summary and request", async () => {
  const clone = vi.spyOn(api, "clone").mockReturnValue(new Promise(() => {}));
  await renderPage({}, { preDigested: { top: "ACGT", bottom: "ACGT", leftEnzyme: "BamHI", rightEnzyme: "EcoRI" } });
  const configuration = screen.getByRole("region", { name: "Cloning configuration" });
  expect(within(configuration).getByText("BamHI / EcoRI")).toBeInTheDocument();
  fireEvent.change(screen.getByLabelText("Left restriction enzyme"), { target: { value: "NdeI" } });
  expect(within(configuration).getByText("NdeI / EcoRI")).toBeInTheDocument();
  fireEvent.click(screen.getByRole("button", { name: "Simulate digestion" }));
  expect(clone).toHaveBeenCalledWith(expect.objectContaining({ left_enzyme: "NdeI", right_enzyme: "EcoRI", pre_digested: true }));
});

it("submits edited vector bases and discards annotations on the previous sequence", async () => {
  const clone = vi.spyOn(api, "clone").mockReturnValue(new Promise(() => {}));
  await renderPage();
  fireEvent.change(screen.getByLabelText("Sequence", { exact: true }), { target: { value: "ACGTACGT" } });
  fireEvent.click(screen.getByRole("button", { name: "Simulate digestion" }));
  expect(clone).toHaveBeenCalledWith(expect.objectContaining({ vector: "ACGTACGT", vector_annotations: [], vector_key: spec.key }));
});

it("ignores a simulation response that arrives after the enzymes change", async () => {
  let finish!: (value: CloneResult) => void;
  vi.spyOn(api, "clone").mockReturnValue(new Promise(resolve => { finish = resolve; }));
  await renderPage();
  fireEvent.click(screen.getByRole("button", { name: "Simulate digestion" }));
  fireEvent.change(screen.getByLabelText("5' enzyme"), { target: { value: "BamHI" } });
  await act(async () => finish({ project_id: 999 } as CloneResult));
  expect(screen.getByText("No plasmid yet")).toBeInTheDocument();
  expect(screen.queryByText(/#999/)).not.toBeInTheDocument();
  expect(screen.getByRole("button", { name: "Simulate digestion" })).toBeEnabled();
});

it("keeps the latest vector when catalogue requests finish out of order", async () => {
  const pending: ((record: VectorRecord) => void)[] = [];
  vi.spyOn(api, "vectorSequence").mockImplementation(() => new Promise(resolve => pending.push(resolve)));
  await renderPage();
  fireEvent.change(screen.getByLabelText("Backbone"), { target: { value: "other" } });
  expect(screen.getByRole("button", { name: "Simulate digestion" })).toBeDisabled();
  fireEvent.change(screen.getByLabelText("Backbone"), { target: { value: spec.key } });
  await act(async () => pending[1]({ ...vector, topology: "circular", sequence: "AAAA" } as unknown as VectorRecord));
  await act(async () => pending[0]({ ...vector, key: "other", topology: "circular", sequence: "CCCC" } as unknown as VectorRecord));
  await waitFor(() => expect(screen.getByLabelText("Sequence", { exact: true })).toHaveValue("AAAA"));
  expect(screen.getByLabelText("Backbone")).toHaveValue(spec.key);
});
