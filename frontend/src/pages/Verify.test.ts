import { describe, expect, it } from "vitest";

import type { Project } from "../api/client";
import { projectCapabilities } from "./Verify";

function project(data: Project["data"], sequence = "A".repeat(200)): Project {
  return {
    id: 1,
    name: "test construct",
    module: "extended_sequence_design",
    sequence,
    notes: "",
    data,
    provenance: {},
    created_at: "2026-08-30T00:00:00Z",
    updated_at: "2026-08-30T00:00:00Z",
  };
}

describe("Check workflow capabilities", () => {
  it("enables reads and primers, but not ligation, for a saved assembly", () => {
    expect(projectCapabilities(project({
      topology: "linear",
      insert_start: 40,
      insert_end: 140,
    }))).toMatchObject({
      hasRegion: true,
      hasLigationContext: false,
      circular: false,
    });
  });

  it("enables ligation only when an explicit backbone length is present", () => {
    expect(projectCapabilities(project({
      topology: "circular",
      insert_start: 40,
      insert_end: 140,
      backbone_length: 5_432,
    }))).toMatchObject({
      backboneLength: 5_432,
      hasRegion: true,
      hasLigationContext: true,
      circular: true,
    });
  });

  it("rejects malformed or out-of-bounds insert coordinates", () => {
    expect(projectCapabilities(project({ insert_start: 140, insert_end: 40 })))
      .toMatchObject({ hasRegion: false, hasLigationContext: false });
    expect(projectCapabilities(project({ insert_start: 40, insert_end: 240 })))
      .toMatchObject({ hasRegion: false, hasLigationContext: false });
  });
});
