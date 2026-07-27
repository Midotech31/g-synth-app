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

  parseFile: (file: File) => {
    const formData = new FormData();
    formData.append("file", file);
    return request<ParsedRecord>("/api/sequences/parse/", { method: "POST", formData });
  },
};
