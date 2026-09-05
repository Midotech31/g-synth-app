import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import type { Enzyme } from "../../api/client";
import EnzymePicker from "../EnzymePicker";

const ENZYMES: Enzyme[] = [
  { name: "NdeI", recognition: "CATATG", overhang: "TA", overhang_type: "5′", supplies_start_codon: true, common: true },
  { name: "XhoI", recognition: "CTCGAG", overhang: "TCGA", overhang_type: "5′", supplies_start_codon: false, common: true },
  { name: "HindIII", aliases: ["Hind3"], recognition: "AAGCTT", overhang: "AGCT", overhang_type: "5′", supplies_start_codon: false, common: false },
];

describe("EnzymePicker", () => {
  it("filters by name, alias or recognition site and keeps the selected enzyme available", () => {
    const onChange = vi.fn();
    render(
      <EnzymePicker
        id="enzyme"
        label="5′ enzyme"
        enzymes={ENZYMES}
        value="NdeI"
        onChange={onChange}
      />,
    );

    fireEvent.change(screen.getByRole("searchbox"), { target: { value: "AAGCTT" } });

    expect(screen.getByRole("option", { name: /HindIII/ })).toBeInTheDocument();
    expect(screen.queryByRole("option", { name: /XhoI/ })).not.toBeInTheDocument();
    expect(screen.getByText("1 matching enzyme")).toBeInTheDocument();

    fireEvent.change(screen.getByRole("combobox"), { target: { value: "HindIII" } });
    expect(onChange).toHaveBeenCalledWith("HindIII");
  });
});
