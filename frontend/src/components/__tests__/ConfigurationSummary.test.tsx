import { render, screen } from "@testing-library/react";
import { expect, it } from "vitest";
import type { CloneResult, DesignParams } from "../../api/client";
import InsertSettingsSummary from "../InsertSettingsSummary";
import VectorConfiguration from "../VectorConfiguration";

it("reports the calculated frame and tags only for the simulated enzyme pair", () => {
  const result = { left_enzyme: "NdeI", right_enzyme: "XhoI", is_clonable: true,
    reading_frame: { summary: "Calculated frame assessment" },
    tags: [{ name: "His-tag", end: "C", present: false }],
  } as CloneResult;
  const props = { spec: null, leftEnzyme: "NdeI", rightEnzyme: "XhoI", result, bundled: false, loading: false };
  const { rerender } = render(<VectorConfiguration {...props} />);
  expect(screen.getByText("Calculated frame assessment")).toBeInTheDocument();
  expect(screen.getByText("not detected in predicted protein")).toBeInTheDocument();
  rerender(<VectorConfiguration {...props} leftEnzyme="BamHI" />);
  expect(screen.queryByText("Calculated frame assessment")).not.toBeInTheDocument();
  expect(screen.queryByText("not detected in predicted protein")).not.toBeInTheDocument();
  expect(screen.getByText(/Awaiting simulation/)).toBeInTheDocument();
});

it("summarizes actual insert options rather than assuming guided defaults", () => {
  const params = { include_his_tag: false, include_linkers: false, cleavage_site: "TEV", target_oligo_length: 120, overhang_length: 6 } as DesignParams;
  const { rerender } = render(<InsertSettingsSummary params={params} />);
  expect(screen.getByText(/no added His tag; no added linkers; TEV site/)).toHaveTextContent("120 nt target oligos with 6 nt assembly junctions");
  rerender(<InsertSettingsSummary params={{ ...params, include_his_tag: true, cleavage_site: null }} fragment={false} />);
  expect(screen.getByText(/6×His tag enabled/)).toHaveTextContent("no protease site");
  expect(screen.queryByText(/120 nt/)).not.toBeInTheDocument();
  expect(screen.getByText(/without fragmentation/)).toBeInTheDocument();
});
