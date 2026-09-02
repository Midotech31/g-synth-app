import { fireEvent, render, screen } from "@testing-library/react";
import { useState } from "react";
import { describe, expect, it } from "vitest";

import PasswordInput from "../PasswordInput";

function PasswordFixture() {
  const [value, setValue] = useState("");
  return (
    <div>
      <label htmlFor="password">Password</label>
      <PasswordInput
        id="password"
        label="Password"
        autoComplete="current-password"
        value={value}
        onChange={setValue}
      />
    </div>
  );
}

describe("PasswordInput", () => {
  it("reveals the entered value and hides it again without clearing it", () => {
    render(<PasswordFixture />);

    const input = screen.getByLabelText("Password");
    const toggle = screen.getByRole("button", { name: "Show password" });

    expect(input).toHaveAttribute("type", "password");
    expect(toggle).toHaveAttribute("aria-pressed", "false");

    fireEvent.change(input, { target: { value: "DemoVisible123" } });
    fireEvent.click(toggle);

    expect(input).toHaveAttribute("type", "text");
    expect(input).toHaveValue("DemoVisible123");
    expect(screen.getByRole("button", { name: "Hide password" })).toHaveAttribute(
      "aria-pressed",
      "true",
    );

    fireEvent.click(screen.getByRole("button", { name: "Hide password" }));

    expect(input).toHaveAttribute("type", "password");
    expect(input).toHaveValue("DemoVisible123");
  });
});
