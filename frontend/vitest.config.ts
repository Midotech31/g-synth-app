import react from "@vitejs/plugin-react";
import { defineConfig } from "vitest/config";

// Separate from vite.config.ts because Vitest reads this file *instead of*
// that one, and the dev proxy there has no meaning under a mocked fetch.
// The React plugin has to be repeated, otherwise .tsx tests get no JSX
// transform.
export default defineConfig({
  plugins: [react()],
  test: {
    // localStorage, the DOM the download helper builds its anchor in, and
    // the FormData a trace upload is assembled with all come from here.
    environment: "jsdom",
    globals: true,
    setupFiles: ["./src/test/setup.ts"],
    include: ["src/**/*.{test,spec}.{ts,tsx}"],
  },
});
