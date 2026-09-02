import { fireEvent, render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { describe, expect, it } from "vitest";

import Learn from "./Learn";

function renderLearn() {
  return render(<MemoryRouter><Learn /></MemoryRouter>);
}

describe("Learn deterministic knowledge library", () => {
  it("is useful without checking or calling an AI service", () => {
    renderLearn();

    expect(screen.getByText("Ask the library, not a black box")).toBeInTheDocument();
    expect(screen.getByText("10 reviewed topics")).toBeInTheDocument();
    expect(screen.getByText("Restriction-enzyme cloning")).toBeInTheDocument();
    expect(screen.getByText("Post-sequencing validation")).toBeInTheDocument();
  });

  it("answers a practical primer-tail search with the reviewed topic", () => {
    renderLearn();
    fireEvent.change(screen.getByLabelText("Search by question or keyword"), {
      target: { value: "5′ primer tail hybridise cycle 1" },
    });

    expect(screen.getByText("Cloning primers and 5′ tails")).toBeInTheDocument();
    expect(screen.getByText(/only the primer's 3′ template-complementary region anneals/i))
      .toBeInTheDocument();
  });

  it("filters to validation topics and can reset", () => {
    renderLearn();
    fireEvent.click(screen.getByRole("button", { name: "Validation" }));

    expect(screen.getByText("Diagnostic digests and agarose gels")).toBeInTheDocument();
    expect(screen.queryByText("Codon optimisation without changing protein")).not.toBeInTheDocument();

    fireEvent.click(screen.getByRole("button", { name: "Reset filters" }));
    expect(screen.getByText("Codon optimisation without changing protein")).toBeInTheDocument();
  });
});
