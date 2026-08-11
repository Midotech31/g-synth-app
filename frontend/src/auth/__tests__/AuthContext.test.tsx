/**
 * What the provider does on boot decides what the user sees before anything
 * else. A stored token can be expired or revoked, so it is checked against
 * the server rather than believed — and every way that check can end has to
 * leave `loading` false, because a `loading` that stays true is a workspace
 * that never appears and gives no reason.
 *
 * The real client is used here, over a mocked fetch: the boot path runs
 * through the refresh logic, and stubbing it out would test the stub.
 */
import { fireEvent, render, screen } from "@testing-library/react";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { AuthProvider, useAuth } from "../AuthContext";

const ACCESS_KEY = "gsynth.access";
const REFRESH_KEY = "gsynth.refresh";

const USER = { id: 1, email: "lab@example.org", name: "Lab", date_joined: "2026-01-01" };

const calls: string[] = [];

function reply(status: number, body: unknown): Response {
  const text = JSON.stringify(body);
  return {
    status,
    ok: status >= 200 && status < 300,
    text: async () => text,
    json: async () => body,
  } as unknown as Response;
}

function serve(handler: (url: string) => Response): void {
  vi.stubGlobal("fetch", (input: RequestInfo | URL) => {
    calls.push(String(input));
    return Promise.resolve(handler(String(input)));
  });
}

/** Says which of the three states the provider is in, in one line. */
function Probe() {
  const { user, loading, signIn } = useAuth();
  if (loading) return <p>checking</p>;
  return (
    <button onClick={() => void signIn("lab@example.org", "correct horse")}>
      {user ? `signed in as ${user.email}` : "signed out"}
    </button>
  );
}

beforeEach(() => {
  calls.length = 0;
  localStorage.clear();
});

describe("AuthProvider on boot", () => {
  it("settles as signed out without asking the server who a visitor is", async () => {
    // Nothing stored means nobody to identify; a request here would answer
    // 401 and greet a first-time visitor with an expired-session message.
    serve(() => reply(500, { detail: "should not be called" }));

    render(
      <AuthProvider>
        <Probe />
      </AuthProvider>,
    );

    expect(await screen.findByText("signed out")).toBeInTheDocument();
    expect(calls).toEqual([]);
  });

  it("shows the signed-in user when the stored token is still good", async () => {
    localStorage.setItem(ACCESS_KEY, "access-1");
    serve(() => reply(200, USER));

    render(
      <AuthProvider>
        <Probe />
      </AuthProvider>,
    );

    expect(await screen.findByText("signed in as lab@example.org")).toBeInTheDocument();
  });

  it("finishes loading and discards the tokens when the stored session is dead", async () => {
    // A revoked pair must not survive boot: kept, it attaches itself to every
    // later request and each one pays a refresh round trip to fail again.
    localStorage.setItem(ACCESS_KEY, "stale");
    localStorage.setItem(REFRESH_KEY, "revoked");
    serve(() => reply(401, { detail: "Token is invalid or expired" }));

    render(
      <AuthProvider>
        <Probe />
      </AuthProvider>,
    );

    expect(await screen.findByText("signed out")).toBeInTheDocument();
    expect(localStorage.getItem(ACCESS_KEY)).toBeNull();
    expect(localStorage.getItem(REFRESH_KEY)).toBeNull();
  });
});

describe("signing in", () => {
  it("stores the tokens and shows who signed in", async () => {
    serve((url) =>
      url.endsWith("/api/auth/login/")
        ? reply(200, { access: "access-1", refresh: "refresh-1" })
        : reply(200, USER),
    );
    render(
      <AuthProvider>
        <Probe />
      </AuthProvider>,
    );
    await screen.findByText("signed out");

    fireEvent.click(screen.getByRole("button"));

    expect(await screen.findByText("signed in as lab@example.org")).toBeInTheDocument();
    expect(localStorage.getItem(ACCESS_KEY)).toBe("access-1");
  });
});

describe("useAuth outside a provider", () => {
  it("names the mistake instead of failing on a null context", () => {
    // React re-throws through a DOM event, so the deliberate failure is
    // announced twice more — once by React's logger and once by jsdom's
    // uncaught-error reporting. Both are silenced: a test run that prints a
    // stack trace while passing teaches the reader to skim the output.
    vi.spyOn(console, "error").mockImplementation(() => {});
    const swallow = (event: ErrorEvent) => event.preventDefault();
    window.addEventListener("error", swallow);

    try {
      expect(() => render(<Probe />)).toThrow("useAuth must be used inside <AuthProvider>");
    } finally {
      window.removeEventListener("error", swallow);
    }
  });
});
