import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";

import packageMetadata from "./package.json" with { type: "json" };

export default defineConfig({


  define: {
    __APP_VERSION__: JSON.stringify(packageMetadata.version),
  },
  plugins: [react()],
  server: {
    port: 5173,


    proxy: {
      "/api": {
        target: process.env.VITE_API_TARGET ?? "http://127.0.0.1:8000",
        changeOrigin: true,
      },
    },
  },
});
