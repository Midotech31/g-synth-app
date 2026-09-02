import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { describe, expect, it, vi } from "vitest";

import { api } from "../api/client";
import { WorkspaceStateProvider } from "../state/WorkspaceStateContext";
import Optimise from "./Optimise";

const HOSTS = [
  {
    key: "ecoli",
    name: "Escherichia coli (high-expression reference)",
    source: "Sharp & Li (1987), relative adaptiveness index",
  },
  {
    key: "s_cerevisiae",
    name: "Saccharomyces cerevisiae",
    source: "Kazusa Codon Usage Database, NCBI taxon 4932",
  },
];

function renderPage() {
  vi.spyOn(api, "codonHosts").mockResolvedValue({ hosts: HOSTS, default: "ecoli" });
  return render(
    <MemoryRouter>
      <WorkspaceStateProvider><Optimise /></WorkspaceStateProvider>
    </MemoryRouter>,
  );
}

describe("host-specific codon optimisation", () => {
  it("loads the host catalogue and exposes the selected table source", async () => {
    renderPage();

    const selector = await screen.findByLabelText("Expression host");
    expect(selector).toHaveValue("ecoli");
    fireEvent.change(selector, { target: { value: "s_cerevisiae" } });

    await waitFor(() => expect(selector).toHaveValue("s_cerevisiae"));
    expect(screen.getByText(/NCBI taxon 4932/)).toBeInTheDocument();
  });
});
