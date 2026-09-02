import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import type { JunctionView } from "../../api/client";
import LigationOutcome from "../LigationOutcome";

const junctions: JunctionView[] = [
  {
    name: "vector → insert",
    enzyme: "NdeI",
    overhang: "TA",
    kind: "5'",
    compatible: true,
    reason: "",
    left_top: "AAA",
    left_bottom: "TTT",
    right_top: "CCC",
    right_bottom: "GGG",
    joined_top: "AAATACCC",
    joined_bottom: "TTTATGGG",
    joined_pairs: "||||||||",
    seam: 3,
    overhang_span: [3, 5],
  },
  {
    name: "insert → vector",
    enzyme: "XhoI",
    overhang: "TCGA",
    kind: "5'",
    compatible: true,
    reason: "",
    left_top: "CCC",
    left_bottom: "GGG",
    right_top: "AAA",
    right_bottom: "TTT",
    joined_top: "CCCTCGAAAA",
    joined_bottom: "GGGAGCTTTT",
    joined_pairs: "||||||||||",
    seam: 3,
    overhang_span: [3, 7],
  },
];

describe("LigationOutcome", () => {
  it("announces a verified product and shows both closed sequence junctions", () => {
    render(
      <LigationOutcome
        vectorName="pET-21a(+)"
        backboneLength={5365}
        insertName="Glargine A"
        insertLength={122}
        productName="pGS-Glargine-A"
        productLength={5487}
        junctions={junctions}
      />,
    );

    const status = screen.getByRole("status", { name: "Ligation result" });
    expect(status).toHaveAttribute("aria-live", "polite");
    expect(screen.getByText("Ligation confirmed")).toBeInTheDocument();
    expect(screen.getByText("pGS-Glargine-A · 5,487 bp")).toBeInTheDocument();
    expect(screen.getByText("2 junctions verified")).toBeInTheDocument();
    expect(screen.getAllByText("Closed")).toHaveLength(2);
    expect(screen.getByLabelText("vector → insert ligated duplex")).toHaveTextContent("AAATACCC");
    expect(screen.getByLabelText("insert → vector ligated duplex")).toHaveTextContent("CCCTCGAAAA");
  });
});
