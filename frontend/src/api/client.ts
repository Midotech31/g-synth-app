/**
 * API client for the G-Synth Django backend.
 *
 * Access tokens are short-lived (30 min in production). Rather than making
 * every caller handle expiry, a 401 triggers one refresh attempt and the
 * original request is replayed. Concurrent 401s share a single refresh
 * promise, so a page issuing five requests doesn't fire five refreshes and
 * invalidate its own rotating refresh token.
 */

const ACCESS_KEY = "gsynth.access";
const REFRESH_KEY = "gsynth.refresh";

export type Tokens = { access: string; refresh: string };

export type User = {
  id: number;
  email: string;
  name: string;
  date_joined: string;
};

export type Annotation = {
  name: string;
  type: string;
  start: number;
  end: number;
  direction: number;
  color: string;
};

export type ParsedRecord = {
  name: string;
  description: string;
  sequence: string;
  length: number;
  topology: "circular" | "linear";
  gc_content: number;
  source_format: string;
  annotations: Annotation[];
};

export type Segment = {
  name: string;
  start: number;
  end: number;
  sequence: string;
};

export type SSDResult = {
  forward: string;
  reverse: string;
  forward_length: number;
  reverse_length: number;
  forward_gc: number;
  reverse_gc: number;
  forward_tm: number;
  reverse_tm: number;
  left_enzyme: string;
  right_enzyme: string;
  left_overhang: string;
  right_overhang: string;
  cleavage_site: string | null;
  orf_start: number;
  coding_region: string;
  segments: Segment[];
  warnings: string[];
  project_id?: number;
};

export type Fragment = {
  index: number;
  name: string;
  forward: string;
  reverse: string;
  forward_length: number;
  reverse_length: number;
  forward_tm: number;
  reverse_tm: number;
  top_start: number;
  top_end: number;
  left_overhang: string;
  right_overhang: string;
  /** Which strand the overhang sits on: "top", "bottom" or "blunt". */
  left_overhang_strand: string;
  right_overhang_strand: string;
  /** How far right the bottom strand's left end sits, in bases. */
  bottom_offset: number;
  is_first: boolean;
  is_last: boolean;
};

export type Oligo = Record<string, string | number>;

export type DuplexSpan = { name: string; start: number; end: number };

/**
 * Both strands in one coordinate frame. A space means that strand is absent
 * from the column — which is what a single-stranded overhang looks like.
 */
export type Duplex = {
  top: string;
  bottom: string;
  pairs: string;
  width: number;
  left_overhang: string;
  right_overhang: string;
  junctions: number[];
  mismatches: number[];
  segments: DuplexSpan[];
  top_fragments: DuplexSpan[];
  bottom_fragments: DuplexSpan[];
};

export type AssemblyResult = {
  construct_forward: string;
  construct_reverse: string;
  construct_length: number;
  construct_gc: number;
  fragment_count: number;
  oligo_count: number;
  overhang_length: number;
  longest_oligo: number;
  junction_overhangs: string[];
  fragments: Fragment[];
  oligos: Oligo[];
  ssd: SSDResult;
  duplex: Duplex;
  tm_conditions: { name: string; summary: string; model: string };
  warnings: string[];
  /** Empty means the oligos re-ligate to the design. Non-empty blocks ordering. */
  verification: string[];
  project_id?: number;
};

export type Junction = {
  name: string;
  enzyme: string;
  overhang: string;
  kind: string;
  position: number;
  context: string;
  site_regenerated: boolean;
};

export type Orf = {
  start: number;
  end: number;
  frame: number;
  codons: number;
  wraps: boolean;
  protein: string;
};

/** The recombinant plasmid: what you actually end up with. */
export type CloneResult = {
  plasmid: string;
  name: string;
  vector_name: string;
  length: number;
  gc: number;
  topology: string;
  insert_start: number;
  insert_end: number;
  insert_length: number;
  backbone_length: number;
  removed_length: number;
  left_enzyme: string;
  right_enzyme: string;
  protein: string;
  protein_length: number;
  /** True when the insert reads on the minus strand of the vector's numbering. */
  reversed_insert: boolean;
  tags: { name: string; end: string; present: boolean; position: number | null; note: string }[];
  vector: { recognised: boolean; spec: VectorSpec | null; check: VectorCheck | null };
  annotations: Annotation[];
  junctions: Junction[];
  orfs: Orf[];
  warnings: string[];
  /** Empty means these two molecules really do join. */
  problems: string[];
  is_clonable: boolean;
  insert: SSDResult;
  assembly: AssemblyResult | null;
  project_id?: number;
};

export type OptimiseParams = {
  sequence: string;
  is_protein?: boolean;
  keep_stop?: boolean;
  avoid_enzymes?: string[];
  avoid_motifs?: string[];
  max_homopolymer?: number;
  gc_min?: number;
  gc_max?: number;
  gc_window?: number;
  max_repeat?: number;
  avoid_rare?: boolean;
  reference_genes?: string[];
};

export type OptimiseResult = {
  sequence: string;
  protein: string;
  length: number;
  table: string;
  table_source: string;
  /** Null when the input was a protein: there was no gene to measure. */
  cai_before: number | null;
  cai_after: number;
  gc_before: number | null;
  gc_after: number;
  sites_removed: string[];
  rare_codons_before: number;
  rare_codons_after: number;
  changed_codons: number;
  /** Empty means the gene can be built and cut as asked. */
  problems: string[];
  warnings: string[];
  is_clean: boolean;
};

export type LigationReaction = {
  ratio: number;
  vector_ng: number;
  insert_ng: number;
  vector_fmol: number;
  insert_fmol: number;
  total_ng: number;
  rows: Record<string, string>[];
  warnings: string[];
};

export type SeqPrimer = {
  name: string;
  sequence: string;
  length: number;
  start: number;
  direction: number;
  tm: number;
  gc: number;
  reads_from: number;
  reads_to: number;
};

export type PrimerSet = {
  primers: SeqPrimer[];
  rows: Record<string, string | number>[];
  target_start: number;
  target_end: number;
  gaps: [number, number][];
  covers_target: boolean;
  warnings: string[];
};

export type Difference = {
  kind: string;
  position: number;
  expected: string;
  found: string;
  residue: number | null;
  from_residue: string;
  to_residue: string;
  silent: boolean | null;
  description: string;
};

export type VerifyReport = {
  design_length: number;
  coverage: number;
  gaps: [number, number][];
  fully_covered: boolean;
  /** Empty differences with at least one read means it is the design. */
  is_verified: boolean;
  differences: Difference[];
  reads: {
    name: string;
    length: number;
    start: number;
    end: number;
    covered: number;
    reverse_complemented: boolean;
    identity: number;
    matched: number;
    difference_count: number;
    is_clean: boolean;
    warnings: string[];
  }[];
  warnings: string[];
};

export type VectorTag = { name: string; end: string; note: string };

/** A backbone G-Synth knows about. `has_sequence` means it ships with one. */
export type VectorSpec = {
  key: string;
  name: string;
  length: number;
  resistance: string;
  promoter: string;
  host: string;
  supplier: string;
  summary: string;
  unique_sites: string[];
  recommended_pairs: string[];
  tags: VectorTag[];
  notes: string[];
  reference: string;
  has_sequence: boolean;
  supplies_translation_start: boolean;
  tag_summary: string;
};

export type VectorRecord = {
  key: string;
  name: string;
  length: number;
  topology: string;
  source: string;
  sequence: string;
  annotations: Annotation[];
  spec: VectorSpec;
};

export type VectorCheck = {
  matches: boolean;
  length: number;
  problems: string[];
  notes: string[];
  found_motifs: string[];
  missing_motifs: string[];
};

export type CloneParams = DesignParams & {
  vector_key?: string;
  vector?: string;
  vector_name?: string;
  vector_annotations?: Annotation[];
  vector_is_circular?: boolean;
  fragment?: boolean;
};

export type Enzyme = {
  name: string;
  recognition: string;
  overhang: string;
  overhang_type: string;
  supplies_start_codon: boolean;
};

export type Catalogue = {
  enzymes: Enzyme[];
  common_pairs: string[];
  cleavage_sites: { name: string; sequence: string }[];
};

export type DesignParams = {
  sequence: string;
  name?: string;
  left_enzyme: string;
  right_enzyme: string;
  is_coding: boolean;
  remove_stop: boolean;
  cleavage_site: string | null;
  include_his_tag: boolean;
  include_linkers: boolean;
  target_oligo_length?: number;
  overhang_length?: number;
  save_as_project?: boolean;
};

export type ProjectSummary = {
  id: number;
  name: string;
  module: string;
  updated_at: string;
};

export type Project = ProjectSummary & {
  sequence: string;
  notes: string;
  data: { annotations?: Annotation[]; topology?: string; gc_content?: number };
  created_at: string;
};

export class ApiError extends Error {
  status: number;
  /** Field-level errors from DRF, e.g. { email: ["already exists"] } */
  fields: Record<string, string[]>;

  constructor(status: number, message: string, fields: Record<string, string[]> = {}) {
    super(message);
    this.status = status;
    this.fields = fields;
  }
}

export const tokenStore = {
  get access() {
    return localStorage.getItem(ACCESS_KEY);
  },
  get refresh() {
    return localStorage.getItem(REFRESH_KEY);
  },
  save({ access, refresh }: Tokens) {
    localStorage.setItem(ACCESS_KEY, access);
    localStorage.setItem(REFRESH_KEY, refresh);
  },
  saveAccess(access: string) {
    localStorage.setItem(ACCESS_KEY, access);
  },
  clear() {
    localStorage.removeItem(ACCESS_KEY);
    localStorage.removeItem(REFRESH_KEY);
  },
};

/** Turn a DRF error body into one readable sentence plus per-field detail. */
function describe(status: number, body: unknown): ApiError {
  if (body && typeof body === "object") {
    const obj = body as Record<string, unknown>;
    if (typeof obj.detail === "string") return new ApiError(status, obj.detail);

    const fields: Record<string, string[]> = {};
    const parts: string[] = [];
    for (const [key, value] of Object.entries(obj)) {
      const messages = Array.isArray(value) ? value.map(String) : [String(value)];
      fields[key] = messages;
      parts.push(key === "non_field_errors" ? messages.join(" ") : `${key}: ${messages.join(" ")}`);
    }
    if (parts.length) return new ApiError(status, parts.join(" · "), fields);
  }
  if (status === 401) return new ApiError(status, "Your session has expired. Please sign in again.");
  if (status === 429) return new ApiError(status, "Too many attempts. Please wait a minute and try again.");
  return new ApiError(status, `Request failed (${status}).`);
}

/** Hand a blob to the browser as a download. */
function saveBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  anchor.click();
  URL.revokeObjectURL(url);
}

let refreshInFlight: Promise<string | null> | null = null;

async function refreshAccessToken(): Promise<string | null> {
  const refresh = tokenStore.refresh;
  if (!refresh) return null;

  const response = await fetch("/api/auth/refresh/", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ refresh }),
  });
  if (!response.ok) return null;

  const data = (await response.json()) as { access: string; refresh?: string };
  // ROTATE_REFRESH_TOKENS is on, so keep the new refresh token when sent.
  if (data.refresh) tokenStore.save({ access: data.access, refresh: data.refresh });
  else tokenStore.saveAccess(data.access);
  return data.access;
}

type RequestOptions = {
  method?: string;
  body?: unknown;
  /** Set for multipart uploads — the browser must write its own boundary. */
  formData?: FormData;
  auth?: boolean;
};

async function request<T>(path: string, options: RequestOptions = {}): Promise<T> {
  const { method = "GET", body, formData, auth = true } = options;

  const send = async (token: string | null): Promise<Response> => {
    const headers: Record<string, string> = {};
    if (!formData) headers["Content-Type"] = "application/json";
    if (auth && token) headers["Authorization"] = `Bearer ${token}`;
    return fetch(path, {
      method,
      headers,
      body: formData ?? (body === undefined ? undefined : JSON.stringify(body)),
    });
  };

  let response = await send(tokenStore.access);

  if (response.status === 401 && auth && tokenStore.refresh) {
    refreshInFlight = refreshInFlight ?? refreshAccessToken().finally(() => {
      refreshInFlight = null;
    });
    const fresh = await refreshInFlight;
    if (fresh) {
      response = await send(fresh);
    } else {
      tokenStore.clear();
      throw new ApiError(401, "Your session has expired. Please sign in again.");
    }
  }

  if (response.status === 204 || response.status === 205) return undefined as T;

  const text = await response.text();
  const payload = text ? JSON.parse(text) : null;

  if (!response.ok) throw describe(response.status, payload);
  return payload as T;
}

export const api = {
  register: (email: string, name: string, password: string, password2: string) =>
    request<User>("/api/auth/register/", {
      method: "POST",
      auth: false,
      body: { email, name, password, password2 },
    }),

  login: async (email: string, password: string) => {
    const tokens = await request<Tokens>("/api/auth/login/", {
      method: "POST",
      auth: false,
      body: { email, password },
    });
    tokenStore.save(tokens);
    return tokens;
  },

  logout: async () => {
    const refresh = tokenStore.refresh;
    try {
      if (refresh) await request<void>("/api/auth/logout/", { method: "POST", body: { refresh } });
    } finally {
      tokenStore.clear();
    }
  },

  me: () => request<User>("/api/auth/me/"),

  listProjects: () =>
    request<{ count: number; results: ProjectSummary[] }>("/api/projects/"),

  getProject: (id: number) => request<Project>(`/api/projects/${id}/`),

  createProject: (payload: Partial<Project>) =>
    request<Project>("/api/projects/", { method: "POST", body: payload }),

  deleteProject: (id: number) =>
    request<void>(`/api/projects/${id}/`, { method: "DELETE" }),

  catalogue: () => request<Catalogue>("/api/design/enzymes/", { auth: false }),

  designSSD: (params: DesignParams) =>
    request<SSDResult>("/api/design/ssd/", { method: "POST", body: params }),

  designAssembly: (params: DesignParams) =>
    request<AssemblyResult>("/api/design/assembly/", { method: "POST", body: params }),

  optimise: (params: OptimiseParams) =>
    request<OptimiseResult>("/api/design/optimise/", { method: "POST", body: params }),

  ligation: (params: {
    vector_length: number;
    insert_length: number;
    vector_ng?: number;
    ends?: string;
    ratios?: number[];
  }) =>
    request<{ reactions: LigationReaction[]; ends: string }>(
      "/api/design/ligation/", { method: "POST", body: params },
    ),

  primers: (params: {
    template: string;
    target_start: number;
    target_end: number;
    circular?: boolean;
    name?: string;
  }) => request<PrimerSet>("/api/design/primers/", { method: "POST", body: params }),

  verify: (params: {
    design: string;
    reads: Record<string, string>;
    circular?: boolean;
    trim?: number;
    coding_start?: number | null;
    coding_end?: number | null;
    region_start?: number | null;
    region_end?: number | null;
  }) => request<VerifyReport>("/api/design/verify/", { method: "POST", body: params }),

  clone: (params: CloneParams) =>
    request<CloneResult>("/api/design/clone/", { method: "POST", body: params }),

  vectors: () =>
    request<{ vectors: VectorSpec[]; default: string }>("/api/design/vectors/", {
      auth: false,
    }),

  vectorSequence: (key: string) =>
    request<VectorRecord>(`/api/design/vectors/${encodeURIComponent(key)}/`, {
      auth: false,
    }),

  /** GET a file the browser saves — for endpoints that need no body. */
  downloadUrl: async (path: string, filename: string) => {
    const response = await fetch(path, {
      headers: { Authorization: `Bearer ${tokenStore.access ?? ""}` },
    });
    if (!response.ok) throw new ApiError(response.status, "Download failed.");
    saveBlob(await response.blob(), filename);
  },

  /** Downloads stream as files, so they bypass the JSON request helper. */
  download: async (path: string, params: DesignParams | CloneParams, filename: string) => {
    const response = await fetch(path, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        Authorization: `Bearer ${tokenStore.access ?? ""}`,
      },
      body: JSON.stringify(params),
    });
    if (!response.ok) throw new ApiError(response.status, "Download failed.");
    saveBlob(await response.blob(), filename);
  },

  parseFile: (file: File) => {
    const formData = new FormData();
    formData.append("file", file);
    return request<ParsedRecord>("/api/sequences/parse/", { method: "POST", formData });
  },
};
