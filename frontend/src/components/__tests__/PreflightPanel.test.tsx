import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import type { PreflightReport } from "../../api/client";
import PreflightPanel from "../PreflightPanel";

const report: PreflightReport = {
  workflow: "cloning",
  verdict: "blocked",
  can_export: false,
  checks: [{
    code: "CLONE_END_COMPATIBILITY",
    label: "Insert and vector ends are compatible",
    status: "block",
    detail: "The observed ends do not anneal.",
    remedy: "Redesign the insert ends.",
    evidence: { left: "TA", right: "TCGA" },
    passed: false,
  }],
  diagnostics: [],
};

describe("PreflightPanel", () => {
  it("states the verdict, stable code and blocking consequence in text", () => {
    render(<PreflightPanel report={report} />);
    expect(screen.getByText("Blocked")).toBeInTheDocument();
    expect(screen.getByText("CLONE_END_COMPATIBILITY")).toBeInTheDocument();
    expect(screen.getByText(/release stay disabled/i)).toBeInTheDocument();
  });
});
