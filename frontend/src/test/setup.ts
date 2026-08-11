import "@testing-library/jest-dom/vitest";

import { cleanup } from "@testing-library/react";
import { afterEach, vi } from "vitest";

// Tokens live in localStorage, which jsdom keeps for the whole file. One test
// leaving a refresh token behind would make the next one's "signed out" case
// pass for the wrong reason, so the store is emptied between tests rather
// than trusted to be empty.
afterEach(() => {
  cleanup();
  localStorage.clear();
  vi.unstubAllEnvs();
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});
