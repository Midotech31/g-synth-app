import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { api } from "../api/client";
import { WorkspaceStateProvider } from "../state/WorkspaceStateContext";
import Optimise from "./Optimise";

const HOSTS = [
  {
    key: "ecoli",
    name: "Escherichia coli",
    source: "FDA HIVE-CUTs/CoCoPUTs September 2021 snapshot",
    category: "Bacteria",
    taxon_id: 562,
    dataset: "RefSeq" as const,
    dataset_release: "September 2021",
    data_scope: "genomic species aggregate including descendant taxa",
    coding_sequences: 130008698,
    codon_count: 39718304295,
    gc_percent: 51.7,
    source_url: "https://dnahive.fda.gov/dna.cgi?cmd=cuts_main",
    metric_label: "Profile-relative CAI",
  },
  {
    key: "s_cerevisiae",
    name: "Saccharomyces cerevisiae",
    source: "FDA HIVE-CUTs/CoCoPUTs September 2021 snapshot",
    category: "Yeasts",
    taxon_id: 4932,
    dataset: "RefSeq" as const,
    dataset_release: "September 2021",
    data_scope: "genomic species aggregate including descendant taxa",
    coding_sequences: 5983,
    codon_count: 2929341,
    gc_percent: 39.63,
    source_url: "https://dnahive.fda.gov/dna.cgi?cmd=cuts_main",
    metric_label: "Profile-relative CAI",
  },
];

const DATASET = {
  name: "FDA HIVE-CUTs / CoCoPUTs",
  release: "September 2021",
  url: "https://dnahive.fda.gov/dna.cgi?cmd=cuts_main",
  sha256: "4c6d35b4c42dd449b9e441e6eef9cb17777211036f0cbc8e7c5f27616b9a8ab0",
};

function renderPage() {
  vi.spyOn(api, "codonHosts").mockResolvedValue({ hosts: HOSTS, default: "ecoli", dataset: DATASET });
  return render(
    <MemoryRouter>
      <WorkspaceStateProvider><Optimise /></WorkspaceStateProvider>
    </MemoryRouter>,
  );
}

beforeEach(() => {
  window.sessionStorage.clear();
  vi.restoreAllMocks();
});

describe("host-specific codon optimisation", () => {
  it("loads the host catalogue and exposes the selected table source", async () => {
    renderPage();

    const selector = await screen.findByLabelText("Expression host");
    expect(selector).toHaveValue("ecoli");
    fireEvent.change(selector, { target: { value: "s_cerevisiae" } });

    await waitFor(() => expect(selector).toHaveValue("s_cerevisiae"));
    expect(screen.getByText(/NCBI taxon 4932/)).toBeInTheDocument();
    expect(screen.getByText(/5,983 CDS/)).toBeInTheDocument();
    expect(screen.getByText(/does not predict expression yield/)).toBeInTheDocument();
  });

  it("exposes a documented custom reference set for context-specific CAI", async () => {
    renderPage();
    await waitFor(() => expect(screen.getByLabelText("Expression host")).not.toBeDisabled());
    fireEvent.click(screen.getByText(/strain-, cell- or tissue-specific/));
    const input = screen.getByLabelText("Highly expressed reference CDSs");
    fireEvent.change(input, { target: { value: "ATGCTGCTG\nATGCTGTTG" } });
    expect(input).toHaveValue("ATGCTGCTG\nATGCTGTTG");
    expect(screen.getByText(/overrides the species profile/i)).toBeInTheDocument();
  });

  it("makes peptide back-translation and initiator-methionine handling explicit", async () => {
    renderPage();
    await waitFor(() => expect(screen.getByLabelText("Expression host")).not.toBeDisabled());

    fireEvent.click(screen.getByRole("button", { name: "Peptide" }));
    expect(screen.getByLabelText("Peptide role")).toHaveValue("auto");
    expect(screen.getByText(/No N-terminal Met: preserve as a mature peptide/)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Back-translate & optimise" })).toBeInTheDocument();

    fireEvent.change(screen.getByLabelText("Peptide role"), {
      target: { value: "complete_orf" },
    });
    expect(screen.getByText(/One initiator Met \(ATG\) will be added/)).toBeInTheDocument();
  });
});
