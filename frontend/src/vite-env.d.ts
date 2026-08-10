/// <reference types="vite/client" />

/**
 * Build-time configuration.
 *
 * Vite substitutes these at build time, so they are baked into the bundle
 * and cannot be changed without rebuilding. Declaring them here is what
 * makes a typo in the name a compile error rather than `undefined` at
 * runtime — which, for the API's address, would mean every request silently
 * going to the wrong origin.
 */
interface ImportMetaEnv {
  /** Absolute origin of the API. Empty in development: Vite proxies /api. */
  readonly VITE_API_BASE?: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}
