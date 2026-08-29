import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { api, type SequenceAnalysis } from "../../api/client";
import { WorkspaceStateProvider } from "../../state/WorkspaceStateContext";
import ReverseComplement from "../ReverseComplement";
import Translate from "../Translate";

const RESULT: SequenceAnalysis = {
  sequence: "CCCATGAAATAAGGG",
  reverse_complement: "CCCTTATTTCATGGG",
  length: 15,
  gc: 46.7,
  minimum_codons: 2,
  frames: [
    { frame: 1, strand: "forward", offset: 0, dna: "CCCATGAAATAAGGG", protein: "PMK*G", first_atg: 3, protein_from_first_atg: "MK*G" },
    { frame: 2, strand: "forward", offset: 1, dna: "CCATGAAATAAGGG", protein: "P*N", first_atg: null, protein_from_first_atg: "" },
    { frame: 3, strand: "forward", offset: 2, dna: "CATGAAATAAGGG", protein: "HEIR", first_atg: null, protein_from_first_atg: "" },
    { frame: -1, strand: "reverse", offset: 0, dna: "CCCTTATTTCATGGG", protein: "PLFH", first_atg: null, protein_from_first_atg: "" },
    { frame: -2, strand: "reverse", offset: 1, dna: "CCTTATTTCATGGG", protein: "LYFM", first_atg: null, protein_from_first_atg: "" },
    { frame: -3, strand: "reverse", offset: 2, dna: "CTTATTTCATGGG", protein: "YFM", first_atg: null, protein_from_first_atg: "" },
  ],
  orfs: [{
    index: 1, frame: 1, strand: "forward", start: 3, end: 12,
    length_nt: 9, amino_acids: 2, dna: "ATGAAATAA", protein: "MK",
    start_codon: "ATG", stop_codon: "TAA",
  }],
};

describe("sequence tools", () => {
  it("renders six-frame results and ORF coordinates returned by the engine", async () => {
    vi.spyOn(api, "analyse").mockResolvedValue(RESULT);
    render(<WorkspaceStateProvider><Translate /></WorkspaceStateProvider>);

    fireEvent.click(screen.getByRole("button", { name: "Analyse" }));

    expect(await screen.findByRole("heading", { name: "Complete ORFs" })).toBeInTheDocument();
    expect(screen.getAllByText("MK").length).toBeGreaterThan(0);
    expect(screen.getByText("4–12")).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "-3" })).toBeInTheDocument();
  });

  it("shows the validated reverse complement returned by the engine", async () => {
    vi.spyOn(api, "analyse").mockResolvedValue(RESULT);
    render(<WorkspaceStateProvider><ReverseComplement /></WorkspaceStateProvider>);

    fireEvent.change(screen.getByLabelText("Sequence (A/C/G/T)"), {
      target: { value: RESULT.sequence },
    });
    fireEvent.click(screen.getByRole("button", { name: "Transform" }));

    expect(await screen.findByText(RESULT.reverse_complement)).toBeInTheDocument();
    expect(screen.getAllByText(/15 nt/).length).toBeGreaterThan(0);
  });
});
