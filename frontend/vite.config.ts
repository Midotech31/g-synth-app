import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";

import packageMetadata from "./package.json" with { type: "json" };

export default defineConfig({
  // The rail displays the exact package version. Keeping this value sourced
  // from package.json prevents a decorative label drifting away from the
  // engine, citation and release metadata.
  define: {
    __APP_VERSION__: JSON.stringify(packageMetadata.version),
  },
  plugins: [react()],
  server: {
    port: 5173,
    // Proxy /api to Django in development so the browser sees one origin —
    // no CORS preflight, and the same relative URLs work in production
    // when the SPA is served behind the same host.
    proxy: {
      "/api": {
        target: process.env.VITE_API_TARGET ?? "http://127.0.0.1:8000",
        changeOrigin: true,
      },
    },
  },
});
