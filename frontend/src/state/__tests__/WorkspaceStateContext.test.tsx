import { fireEvent, render, screen } from "@testing-library/react";
import { useState } from "react";
import { describe, expect, it } from "vitest";

import {
  WorkspaceStateProvider,
  useWorkspaceState,
} from "../WorkspaceStateContext";

function Workspace() {
  const [sequence, setSequence, clearSequence] = useWorkspaceState(
    "test.sequence",
    "sample",
  );
  return (
    <>
      <label htmlFor="sequence">Sequence</label>
      <input
        id="sequence"
        value={sequence}
        onChange={(event) => setSequence(event.target.value)}
      />
      <button onClick={clearSequence}>Clear workspace</button>
    </>
  );
}

function NavigationHarness() {
  const [open, setOpen] = useState(true);
  return (
    <WorkspaceStateProvider>
      <button onClick={() => setOpen((current) => !current)}>Change section</button>
      {open ? <Workspace /> : <div>Another section</div>}
    </WorkspaceStateProvider>
  );
}

describe("workspace state", () => {
  it("survives a section unmount and returns to its default only when cleared", () => {
    render(<NavigationHarness />);

    fireEvent.change(screen.getByLabelText("Sequence"), {
      target: { value: "my-result" },
    });
    fireEvent.click(screen.getByRole("button", { name: "Change section" }));
    fireEvent.click(screen.getByRole("button", { name: "Change section" }));

    expect(screen.getByLabelText("Sequence")).toHaveValue("my-result");

    fireEvent.click(screen.getByRole("button", { name: "Clear workspace" }));
    expect(screen.getByLabelText("Sequence")).toHaveValue("sample");
  });
});
