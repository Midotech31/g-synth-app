import { act, fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { Link, MemoryRouter, Route, Routes } from "react-router-dom";
import { beforeEach, expect, it, vi } from "vitest";
import { api, type Project, type DetectedFeature } from "../api/client";
import Viewer from "./Viewer";

vi.mock("seqviz", () => ({ SeqViz: () => <div>Sequence map</div> }));
const project = (id: number): Project => ({ id, name: `Record ${id}`, module: "general",
  sequence: "ACGT".repeat(10), notes: "", data: { topology: "circular", annotations: [] },
  created_at: "2026-09-11T00:00:00Z", updated_at: "2026-09-11T00:00:00Z" });
const match: DetectedFeature = { annotation: { name: "Old candidate", type: "misc_feature", start: 0,
  end: 4, direction: 1, color: "#000000", inferred: true }, matched_sequence: "ACGT", basis: "Review required" };

function open() {
  render(<MemoryRouter initialEntries={["/projects/1"]}><Link to="/projects/2">Next project</Link>
    <Routes><Route path="/projects/:id" element={<Viewer />} /></Routes></MemoryRouter>);
}
beforeEach(() => {
  vi.restoreAllMocks();
  vi.spyOn(api, "getProject").mockImplementation(async id => project(id));
  vi.spyOn(api, "detectCommonFeatures").mockResolvedValue({ matches: [], method: "test" });
});

it("does not replace the next project with a late annotation save", async () => {
  let finish!: (value: Project) => void;
  const save = vi.spyOn(api, "updateProjectAnnotations").mockReturnValue(new Promise(resolve => { finish = resolve; }));
  open();
  await screen.findByRole("heading", { name: "Record 1" });
  fireEvent.click(screen.getByRole("button", { name: "Add feature" }));
  fireEvent.change(screen.getByLabelText("Feature name"), { target: { value: "New feature" } });
  fireEvent.click(within(screen.getByRole("dialog")).getByRole("button", { name: "Add feature" }));
  await waitFor(() => expect(save).toHaveBeenCalledWith(1, expect.any(Array), project(1).updated_at));
  fireEvent.click(screen.getByRole("link", { name: "Next project" }));
  await screen.findByRole("heading", { name: "Record 2" });
  await act(async () => finish(project(1)));
  expect(screen.getByRole("heading", { name: "Record 2" })).toBeInTheDocument();
  expect(screen.queryByRole("dialog")).not.toBeInTheDocument();
});

it("ignores a manual scan that finishes after navigation", async () => {
  let finish!: (value: { matches: DetectedFeature[]; method: string }) => void;
  vi.mocked(api.detectCommonFeatures).mockResolvedValueOnce({ matches: [], method: "test" })
    .mockImplementationOnce(() => new Promise(resolve => { finish = resolve; }));
  open();
  await screen.findByRole("button", { name: "Find common" });
  fireEvent.click(screen.getByRole("button", { name: "Find common" }));
  fireEvent.click(screen.getByRole("link", { name: "Next project" }));
  await screen.findByRole("heading", { name: "Record 2" });
  await act(async () => finish({ matches: [match], method: "test" }));
  expect(screen.queryByText("Old candidate")).not.toBeInTheDocument();
});

it("requires choosing candidates explicitly after a manual scan", async () => {
  vi.mocked(api.detectCommonFeatures).mockResolvedValue({ matches: [match], method: "test" });
  open();
  await screen.findByText("Old candidate");
  fireEvent.click(screen.getByRole("button", { name: "Find common" }));
  await waitFor(() => expect(screen.getByRole("button", { name: "Add selected (0)" })).toBeDisabled());
});
