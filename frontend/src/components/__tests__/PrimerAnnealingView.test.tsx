import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import type { PcrPrimer } from "../../api/client";
import PrimerAnnealingView, { annealingRows } from "../PrimerAnnealingView";

const primer: PcrPrimer = {
  name: "insert_F",
  sequence: "GCGCGCAAGCTTATGCCGTA",
  tail: "GCGCGCAAGCTT",
  anneals: "ATGCCGTA",
  direction: 1,
  start: 10,
  end: 18,
  length: 20,
  anneal_length: 8,
  tm: 58,
  tm_full: 72,
  gc: 55,
  enzyme: "HindIII",
  restriction_site: "AAGCTT",
  has_gc_clamp: true,
  warnings: [],
};

describe("primer/template annealing view", () => {
  it("leaves the 5-prime cloning tail visibly unpaired", () => {
    const rows = annealingRows(primer);
    expect(rows.pairs.slice(0, primer.tail.length)).toBe(" ".repeat(primer.tail.length));
    expect(rows.pairs.slice(primer.tail.length)).toBe("|".repeat(primer.anneals.length));
    expect(rows.template.trimStart()).toBe("TACGGCAT");
  });

  it("states that only the 3-prime portion hybridizes in cycle one", () => {
    render(<PrimerAnnealingView forward={primer} reverse={{ ...primer, name: "insert_R" }} />);
    expect(screen.getAllByText(/intentionally unpaired in cycle 1/i)).toHaveLength(2);
    expect(screen.getAllByText(/only the 8-nt 3′ region hybridizes/i)).toHaveLength(2);
    expect(screen.getByText(/tail-derived restriction sites now have complements/i)).toBeInTheDocument();
    expect(screen.getAllByText("AAGCTT")).toHaveLength(2);
    expect(screen.getAllByText("AAGCTT")[0]).toHaveClass("annealing-site");
  });
});
