import { defineConfig } from "vitest/config";

// Separate from vite.config.ts because Vitest reads this file *instead of*
// that one, and the dev proxy there has no meaning under a mocked fetch.
// No React plugin: JSX in tests is transformed from tsconfig's "react-jsx",
// and Fast Refresh — the rest of what the plugin does — has nothing to
// refresh in a test run.
export default defineConfig({
  test: {
    // localStorage, the DOM the download helper builds its anchor in, and
    // the FormData a trace upload is assembled with all come from here.
    environment: "jsdom",
    globals: true,
    setupFiles: ["./src/test/setup.ts"],
    include: ["src/**/*.{test,spec}.{ts,tsx}"],
  },
});
