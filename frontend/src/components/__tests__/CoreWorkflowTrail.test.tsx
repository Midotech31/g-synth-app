import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { describe, expect, it } from "vitest";

import CoreWorkflowTrail from "../CoreWorkflowTrail";

describe("CoreWorkflowTrail", () => {
  it("keeps Design, Hybridization and Restriction cloning in scientific order", () => {
    render(
      <MemoryRouter>
        <CoreWorkflowTrail active="hybridization" />
      </MemoryRouter>,
    );

    const links = screen.getAllByRole("link");
    expect(links.map((link) => link.textContent)).toEqual([
      "1. DesignSSD / ESD molecules",
      "2. HybridizationVerify both strands",
      "3. Restriction cloningDigest and ligate",
    ]);
    expect(links[1]).toHaveAttribute("aria-current", "step");
    expect(links.map((link) => link.getAttribute("href"))).toEqual([
      "/design",
      "/hybridize",
      "/clone",
    ]);
  });
});
