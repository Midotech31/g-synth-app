import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { api, type PcrPrimer, type PcrResult } from "../api/client";
import { WorkspaceStateProvider } from "../state/WorkspaceStateContext";
import Pcr from "./Pcr";

const TEMPLATE = "ATGACCACCAGCAAACTGGGCAAAGGCCTGGGCTATATTGGCAACAACGGCGCGCACATGGGCTTAAACTTAGCATTACTGGGCCTGGCGAGCCTGCTGGGCAAAGGCATTAGCAAACTGGGC";

const primer = (direction: 1 | -1): PcrPrimer => ({
  name: direction === 1 ? "product_F" : "product_R",
  sequence: direction === 1
    ? "GGAGGTCATATGACCACCAGCAAACTGGGC"
    : "GGAGGTCTCGAGGCCCAGTTTGCTAATGCC",
  tail: direction === 1 ? "GGAGGTCATATG" : "GGAGGTCTCGAG",
  anneals: direction === 1 ? "ACCACCAGCAAACTGGGC" : "GCCCAGTTTGCTAATGCC",
  direction,
  start: direction === 1 ? 3 : 105,
  end: direction === 1 ? 21 : 123,
  length: 30,
  anneal_length: 18,
  tm: 60,
  tm_full: 73,
  gc: 57,
  enzyme: direction === 1 ? "NdeI" : "XhoI",
  restriction_site: direction === 1 ? "CATATG" : "CTCGAG",
  has_gc_clamp: true,
  warnings: [],
});

const result = (source: "automatic" | "custom"): PcrResult => ({
  forward: primer(1),
  reverse: primer(-1),
  product: TEMPLATE,
  product_length: TEMPLATE.length,
  amplified_region: TEMPLATE.slice(3),
  template_start: 3,
  template_end: TEMPLATE.length,
  annealing_temperature: 55,
  left_enzyme: "NdeI",
  right_enzyme: "XhoI",
  insert_orf_start: 0,
  problems: [],
  warnings: [],
  is_clean: true,
  primer_source: source,
  digest: null,
});

function renderPage() {
  vi.spyOn(api, "catalogue").mockResolvedValue({
    enzymes: [
      { name: "NdeI", recognition: "CATATG", overhang: "TA", overhang_type: "5'", supplies_start_codon: true, common: true },
      { name: "XhoI", recognition: "CTCGAG", overhang: "TCGA", overhang_type: "5'", supplies_start_codon: false, common: true },
    ],
    common_pairs: [],
    cleavage_sites: [],
  });
  const pcr = vi.spyOn(api, "pcr")
    .mockResolvedValueOnce(result("automatic"))
    .mockResolvedValueOnce(result("custom"));
  render(
    <MemoryRouter>
      <WorkspaceStateProvider><Pcr /></WorkspaceStateProvider>
    </MemoryRouter>,
  );
  return pcr;
}

beforeEach(() => {
  window.sessionStorage.clear();
  vi.restoreAllMocks();
});

describe("editable PCR primers", () => {
  it("requires edited primers to be revalidated before presenting them as validated", async () => {
    const pcr = renderPage();
    fireEvent.change(screen.getByLabelText("Sequence to amplify (A/C/G/T)"), {
      target: { value: TEMPLATE },
    });
    fireEvent.click(screen.getByRole("button", { name: "Design PCR" }));

    await screen.findByText("G-Synth design");
    fireEvent.click(screen.getByRole("button", { name: "Edit primers" }));
    const forward = screen.getByLabelText("Forward primer (5′→3′)");
    fireEvent.change(forward, { target: { value: "AACCGGCATATGACCACCAGCAAACTGGGC" } });

    expect(screen.getByText("Edited primers need validation")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "Revalidate primers" }));

    await screen.findByText("Edited · validated");
    await waitFor(() => expect(pcr).toHaveBeenCalledTimes(2));
    expect(pcr.mock.calls[1][0]).toMatchObject({
      forward_primer: "AACCGGCATATGACCACCAGCAAACTGGGC",
      reverse_primer: "GGAGGTCTCGAGGCCCAGTTTGCTAATGCC",
    });
  });
});
