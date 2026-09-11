import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

import type { ApiError, DesignParams } from "./client";

type Client = typeof import("./client");
type Call = { url: string; init: RequestInit };
type Handler = (url: string, init: RequestInit) => Response | Promise<Response>;

const ACCESS_KEY = "gsynth.access";
const REFRESH_KEY = "gsynth.refresh";

const calls: Call[] = [];
let client: Client;


function reply(status: number, body: string | null = null): Response {
  return {
    status,
    ok: status >= 200 && status < 300,
    text: async () => body ?? "",
    json: async () => JSON.parse(body ?? "null") as unknown,
    blob: async () => new Blob([body ?? ""]),
  } as unknown as Response;
}

function json(status: number, body: unknown): Response {
  return reply(status, JSON.stringify(body));
}

function serve(handler: Handler): void {
  vi.stubGlobal("fetch", (input: RequestInfo | URL, init: RequestInit = {}) => {
    calls.push({ url: String(input), init });
    return Promise.resolve(handler(String(input), init));
  });
}


function header(init: RequestInit, name: string): string | null {
  const headers = init.headers;
  if (!headers) return null;
  if (typeof (headers as Headers).get === "function") return (headers as Headers).get(name);
  return (headers as Record<string, string>)[name] ?? null;
}

const isRefresh = (url: string) => url.endsWith("/api/auth/refresh/");
const refreshCalls = () => calls.filter((call) => isRefresh(call.url));


const flush = () => new Promise((resolve) => setTimeout(resolve, 0));

function deferred(): { promise: Promise<void>; release: () => void } {
  let release!: () => void;
  const promise = new Promise<void>((resolve) => {
    release = resolve;
  });
  return { promise, release };
}


async function rejection(promise: Promise<unknown>): Promise<ApiError> {
  const error = await promise.then(
    () => {
      throw new Error("expected the request to fail, but it resolved");
    },
    (reason: unknown) => reason,
  );
  expect(error).toBeInstanceOf(client.ApiError);
  return error as ApiError;
}


async function loadClient(base = ""): Promise<Client> {
  vi.stubEnv("VITE_API_BASE", base);
  vi.resetModules();
  return import("./client");
}

const USER = { id: 1, email: "lab@example.org", name: "Lab", date_joined: "2026-01-01" };

const DESIGN: DesignParams = {
  sequence: "ATGGGCAGCAGCCATCACCATCACCATCACTAA",
  left_enzyme: "NdeI",
  right_enzyme: "XhoI",
  is_coding: true,
  remove_stop: true,
  cleavage_site: null,
  include_his_tag: true,
  include_linkers: false,
};

beforeEach(async () => {
  calls.length = 0;
  localStorage.clear();
  client = await loadClient();
});

describe("tokenStore", () => {
  it("gives back the pair it was handed, so a reload does not ask for the password again", () => {
    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });

    expect(client.tokenStore.access).toBe("access-1");
    expect(client.tokenStore.refresh).toBe("refresh-1");
  });

  it("reports an empty store as null, which is what the 'is there a session' guards test", () => {
    expect(client.tokenStore.access).toBeNull();
    expect(client.tokenStore.refresh).toBeNull();
  });

  it("clears both halves: a refresh token left behind can mint a new session after sign-out", () => {
    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });

    client.tokenStore.clear();

    expect(localStorage.getItem(ACCESS_KEY)).toBeNull();
    expect(localStorage.getItem(REFRESH_KEY)).toBeNull();
  });

  it("keeps the refresh token when only the access token is replaced", () => {


    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });

    client.tokenStore.saveAccess("access-2");

    expect(client.tokenStore.access).toBe("access-2");
    expect(client.tokenStore.refresh).toBe("refresh-1");
  });
});

describe("error bodies", () => {
  it("passes DRF's detail through word for word", async () => {
    serve(() => json(400, { detail: "NdeI cuts inside the insert at 412." }));

    const error = await rejection(client.api.catalogue());

    expect(error.status).toBe(400);
    expect(error.message).toBe("NdeI cuts inside the insert at 412.");
  });

  it("names every field that failed, and keeps the per-field lists a form needs", async () => {
    serve(() =>
      json(400, {
        email: ["A user with that address already exists."],
        password: ["This password is too short."],
      }),
    );

    const error = await rejection(client.api.catalogue());

    expect(error.status).toBe(400);
    expect(error.message).toBe(
      "email: A user with that address already exists. · password: This password is too short.",
    );
    expect(error.fields).toEqual({
      email: ["A user with that address already exists."],
      password: ["This password is too short."],
    });
  });

  it("shows non_field_errors without its key, which means nothing to the reader", async () => {
    serve(() => json(400, { non_field_errors: ["The two passwords do not match."] }));

    const error = await rejection(client.api.catalogue());

    expect(error.message).toBe("The two passwords do not match.");
    expect(error.fields.non_field_errors).toEqual(["The two passwords do not match."]);
  });

  it("wraps a bare string field in a list, so the form renders it like any other", async () => {
    serve(() => json(400, { sequence: "Not a nucleotide sequence." }));

    const error = await rejection(client.api.catalogue());

    expect(error.fields.sequence).toEqual(["Not a nucleotide sequence."]);
    expect(error.message).toBe("sequence: Not a nucleotide sequence.");
  });

  it("turns a bodyless 401 into a sentence that says what to do about it", async () => {
    serve(() => json(401, {}));

    const error = await rejection(client.api.catalogue());

    expect(error.status).toBe(401);
    expect(error.message).toBe("Your session has expired. Please sign in again.");
  });

  it("turns a bodyless 429 into the wait, not into a bare status code", async () => {
    serve(() => json(429, {}));

    const error = await rejection(client.api.catalogue());

    expect(error.status).toBe(429);
    expect(error.message).toBe("Too many attempts. Please wait a minute and try again.");
  });

  it("prefers the throttle's own detail, because only it knows how long", async () => {
    serve(() => json(429, { detail: "Request was throttled. Expected available in 47 seconds." }));

    const error = await rejection(client.api.catalogue());

    expect(error.message).toBe("Request was throttled. Expected available in 47 seconds.");
  });

  it("survives a proxy's HTML error page instead of reporting a parse error", async () => {


    serve(() => reply(502, "<html><body>502 Bad Gateway</body></html>"));

    const error = await rejection(client.api.catalogue());

    expect(error.status).toBe(502);
    expect(error.message).toBe("Request failed (502).");
  });

  it("calls a 200 that is not JSON a server fault, not a design fault", async () => {
    serve(() => reply(200, "<html>login page</html>"));

    const error = await rejection(client.api.catalogue());

    expect(error.status).toBe(502);
    expect(error.message).toBe("The server returned an invalid response.");
  });

  it("accepts an empty 204 without trying to parse a body", async () => {
    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });
    serve(() => reply(204));

    await expect(client.api.deleteProject(4)).resolves.toBeUndefined();
  });
});

describe("refresh on 401", () => {
  it("refreshes once and replays the original request with the new token", async () => {
    client.tokenStore.save({ access: "stale", refresh: "refresh-1" });
    serve((url, init) => {
      if (isRefresh(url)) return json(200, { access: "fresh" });
      return header(init, "Authorization") === "Bearer fresh"
        ? json(200, USER)
        : json(401, { detail: "Given token not valid for any token type" });
    });

    const me = await client.api.me();

    expect(me.email).toBe("lab@example.org");
    expect(calls.map((call) => call.url)).toEqual([
      "/api/auth/me/",
      "/api/auth/refresh/",
      "/api/auth/me/",
    ]);
    expect(JSON.parse(String(calls[1].init.body))).toEqual({ refresh: "refresh-1" });
    expect(client.tokenStore.access).toBe("fresh");
  });

  it("answers three simultaneous expiries with one refresh, not three", async () => {


    client.tokenStore.save({ access: "stale", refresh: "refresh-1" });
    const gate = deferred();
    serve(async (url, init) => {
      if (isRefresh(url)) {
        await gate.promise;
        return json(200, { access: "fresh", refresh: "refresh-2" });
      }
      return header(init, "Authorization") === "Bearer fresh"
        ? json(200, USER)
        : json(401, { detail: "Given token not valid for any token type" });
    });

    const pending = Promise.all([
      client.api.me(),
      client.api.listProjects(),
      client.api.getProject(7),
    ]);


    await flush();
    expect(refreshCalls()).toHaveLength(1);
    expect(calls).toHaveLength(4);

    gate.release();
    await pending;

    expect(refreshCalls()).toHaveLength(1);
    expect(calls).toHaveLength(7);

    const replays = calls.slice(4);
    expect(replays.map((call) => call.url).sort()).toEqual(
      ["/api/auth/me/", "/api/projects/", "/api/projects/7/"].sort(),
    );
    expect(replays.every((call) => header(call.init, "Authorization") === "Bearer fresh")).toBe(true);
    expect(client.tokenStore.refresh).toBe("refresh-2");
  });

  it("refreshes again when a later token expires — the latch must not stay shut", async () => {
    client.tokenStore.save({ access: "stale", refresh: "refresh-1" });
    let issued = 0;
    serve((url, init) => {
      if (isRefresh(url)) {
        issued += 1;
        return json(200, { access: `fresh-${issued}` });
      }
      return header(init, "Authorization") === `Bearer fresh-${issued}` && issued > 0
        ? json(200, USER)
        : json(401, { detail: "Given token not valid for any token type" });
    });

    await client.api.me();


    localStorage.setItem(ACCESS_KEY, "stale-again");
    issued = 0;
    await client.api.me();

    expect(refreshCalls()).toHaveLength(2);
    expect(client.tokenStore.access).toBe("fresh-1");
  });

  it("clears both tokens and reports a 401 when the refresh itself is rejected", async () => {
    client.tokenStore.save({ access: "stale", refresh: "revoked" });
    serve((url) =>
      isRefresh(url)
        ? json(401, { detail: "Token is invalid or expired" })
        : json(401, { detail: "Given token not valid for any token type" }),
    );

    const error = await rejection(client.api.me());

    expect(error.status).toBe(401);
    expect(error.message).toBe("Your session has expired. Please sign in again.");
    expect(localStorage.getItem(ACCESS_KEY)).toBeNull();
    expect(localStorage.getItem(REFRESH_KEY)).toBeNull();


    expect(calls).toHaveLength(2);
  });

  it("keeps the old refresh token when the server sends no new one", async () => {
    client.tokenStore.save({ access: "stale", refresh: "refresh-1" });
    serve((url, init) =>
      isRefresh(url)
        ? json(200, { access: "fresh" })
        : header(init, "Authorization") === "Bearer fresh"
          ? json(200, USER)
          : json(401, {}),
    );

    await client.api.me();

    expect(client.tokenStore.access).toBe("fresh");
    expect(client.tokenStore.refresh).toBe("refresh-1");
  });

  it("does not attempt a refresh it has no token for", async () => {
    localStorage.setItem(ACCESS_KEY, "stale");
    serve(() => json(401, { detail: "Given token not valid for any token type" }));

    const error = await rejection(client.api.me());

    expect(error.status).toBe(401);
    expect(error.message).toBe("Given token not valid for any token type");
    expect(calls).toHaveLength(1);
  });

  it("does not refresh on failures that are not about the token", async () => {


    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });
    serve(() => json(403, { detail: "You do not have permission to perform this action." }));

    const error = await rejection(client.api.getProject(7));

    expect(error.status).toBe(403);
    expect(refreshCalls()).toHaveLength(0);
    expect(calls).toHaveLength(1);
  });

  it("sends no Authorization header on the endpoints that need no account", async () => {


    localStorage.setItem(ACCESS_KEY, "stale");
    serve(() => json(200, { enzymes: [], common_pairs: [], cleavage_sites: [] }));

    await client.api.catalogue();

    expect(header(calls[0].init, "Authorization")).toBeNull();
  });

  it("loads codon-host reference data without attaching account tokens", async () => {
    localStorage.setItem(ACCESS_KEY, "stale");
    serve(() => json(200, { hosts: [], default: "ecoli" }));

    await client.api.codonHosts();

    expect(calls[0].url.endsWith("/api/design/codon-hosts/")).toBe(true);
    expect(header(calls[0].init, "Authorization")).toBeNull();
  });
});

describe("API_BASE", () => {
  it("strips trailing slashes, so the path never doubles one", async () => {


    client = await loadClient("https://api.gsynth.test///");
    serve(() => json(200, { enzymes: [], common_pairs: [], cleavage_sites: [] }));

    await client.api.catalogue();

    expect(calls[0].url).toBe("https://api.gsynth.test/api/design/enzymes/");
  });

  it("leaves an absolute base without a trailing slash alone", async () => {
    client = await loadClient("https://api.gsynth.test");
    serve(() => json(200, { enzymes: [], common_pairs: [], cleavage_sites: [] }));

    await client.api.catalogue();

    expect(calls[0].url).toBe("https://api.gsynth.test/api/design/enzymes/");
  });

  it("turns Render's bare service host into an HTTPS origin", async () => {
    client = await loadClient("gsynth-api.onrender.com");
    serve(() => json(200, { enzymes: [], common_pairs: [], cleavage_sites: [] }));

    await client.api.catalogue();

    expect(calls[0].url).toBe("https://gsynth-api.onrender.com/api/design/enzymes/");
  });

  it("keeps the path relative when no base is set, which is what the dev proxy needs", async () => {
    client = await loadClient("");
    serve(() => json(200, { enzymes: [], common_pairs: [], cleavage_sites: [] }));

    await client.api.catalogue();

    expect(calls[0].url).toBe("/api/design/enzymes/");
  });
});

describe("session endpoints", () => {
  it("stores the pair a sign-in returns, so the next request carries it", async () => {
    serve((url, init) =>
      url.endsWith("/api/auth/login/")
        ? json(200, { access: "access-1", refresh: "refresh-1" })
        : header(init, "Authorization") === "Bearer access-1"
          ? json(200, USER)
          : json(401, {}),
    );

    await client.api.login("lab@example.org", "correct horse");
    const me = await client.api.me();

    expect(me.email).toBe("lab@example.org");
  });

  it("stores nothing when the credentials are refused", async () => {
    serve(() => json(401, { detail: "No active account found with the given credentials" }));

    const error = await rejection(client.api.login("lab@example.org", "wrong"));

    expect(error.message).toBe("No active account found with the given credentials");
    expect(localStorage.getItem(ACCESS_KEY)).toBeNull();
    expect(localStorage.getItem(REFRESH_KEY)).toBeNull();
  });

  it("signs out locally even when the server cannot be told", async () => {


    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });
    serve(() => json(500, { detail: "Blacklist unavailable." }));

    await rejection(client.api.logout());

    expect(localStorage.getItem(ACCESS_KEY)).toBeNull();
    expect(localStorage.getItem(REFRESH_KEY)).toBeNull();

    expect(JSON.parse(String(calls[0].init.body))).toEqual({ refresh: "refresh-1" });
  });
});

describe("multipart uploads", () => {
  it("sets no Content-Type on a trace upload, so the browser writes the boundary", async () => {


    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });
    serve(() => json(200, { differences: [] }));

    await client.api.verifyTraces({
      design: "ATGGGCAGC",
      files: [new File(["trace"], "forward.ab1"), new File(["trace"], "reverse.ab1")],
      circular: false,
      trim_quality: 20,
    });

    expect(header(calls[0].init, "Content-Type")).toBeNull();
    const form = calls[0].init.body as FormData;
    expect(form.get("design")).toBe("ATGGGCAGC");
    expect(form.getAll("traces")).toHaveLength(2);
    expect(form.get("trim_quality")).toBe("20");


    expect(form.get("circular")).toBe("false");

    expect(form.get("coding_start")).toBeNull();
  });

  it("sends a parsed file under the name the endpoint reads", async () => {
    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });
    serve(() => json(200, { name: "pET-21a", length: 5443 }));

    await client.api.parseFile(new File(["LOCUS"], "pET-21a.gb"));

    const form = calls[0].init.body as FormData;
    expect((form.get("file") as File).name).toBe("pET-21a.gb");
  });
});

describe("PCR design", () => {
  it("sends the explicit NdeI start-codon policy to the engine", async () => {
    serve(() => json(200, {}));

    await client.api.pcr({
      template: "ATGAAAGGTGAAGAATTGTTCACCGGTGTTGTTCCGATTCTG",
      left_enzyme: "NdeI",
      right_enzyme: "XhoI",
      start_codon_mode: "use_site",
    });

    expect(JSON.parse(String(calls[0].init.body))).toMatchObject({
      left_enzyme: "NdeI",
      right_enzyme: "XhoI",
      start_codon_mode: "use_site",
    });
  });

  it("sends edited primer sequences for server-side revalidation", async () => {
    serve(() => json(200, {}));

    await client.api.pcr({
      template: "ATGAAAGGTGAAGAATTGTTCACCGGTGTTGTTCCGATTCTG",
      forward_primer: "AACCGGCATATGAAAGGTGAAGAATTGTT",
      reverse_primer: "TTGGCCCTCGAGCAGAATCGGAACAACAC",
    });

    expect(JSON.parse(String(calls[0].init.body))).toMatchObject({
      forward_primer: "AACCGGCATATGAAAGGTGAAGAATTGTT",
      reverse_primer: "TTGGCCCTCGAGCAGAATCGGAACAACAC",
    });
  });
});

describe("hybridization analysis", () => {
  it("sends both strands in 5-prime order with explicit thermodynamic conditions", async () => {
    serve(() => json(200, {}));

    await client.api.hybridize({
      first: "AATTATGCGT",
      second: "ACGCAT",
      analysis_temperature_c: 25,
      oligo_nM: 50_000,
      na_mM: 50,
      mg_mM: 1.5,
      dntp_mM: 0.2,
    });

    expect(calls[0].url).toContain("/api/design/hybridize/");
    expect(JSON.parse(String(calls[0].init.body))).toEqual({
      first: "AATTATGCGT",
      second: "ACGCAT",
      analysis_temperature_c: 25,
      oligo_nM: 50_000,
      na_mM: 50,
      mg_mM: 1.5,
      dntp_mM: 0.2,
    });
  });
});

describe("file downloads", () => {
  const revoke = vi.fn();

  beforeEach(() => {


    URL.createObjectURL = vi.fn(() => "blob:gsynth/1");
    URL.revokeObjectURL = revoke;
    revoke.mockClear();
  });

  afterEach(() => {
    expect(document.querySelectorAll("a")).toHaveLength(0);
  });

  it("hands the blob to the browser under the name asked for, then tidies up", async () => {
    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });
    serve(() => reply(200, "LOCUS pET-21a-insert 6031 bp DNA circular"));
    const clicked: { href: string; download: string; attached: boolean }[] = [];
    vi.spyOn(HTMLAnchorElement.prototype, "click").mockImplementation(function (
      this: HTMLAnchorElement,
    ) {


      clicked.push({ href: this.href, download: this.download, attached: this.isConnected });
    });

    await client.api.download("/api/design/assembly/genbank/", DESIGN, "assembly.gb");

    expect(clicked).toEqual([
      { href: "blob:gsynth/1", download: "assembly.gb", attached: true },
    ]);
    await vi.waitFor(() => expect(revoke).toHaveBeenCalledWith("blob:gsynth/1"));
  });

  it("throws rather than saving the server's error page as a GenBank file", async () => {


    client.tokenStore.save({ access: "access-1", refresh: "refresh-1" });
    serve(() => reply(500, "<html><body>Server Error</body></html>"));

    const error = await rejection(
      client.api.download("/api/design/assembly/genbank/", DESIGN, "assembly.gb"),
    );

    expect(error.status).toBe(500);
    expect(error.message).toBe("Download failed.");
    expect(URL.createObjectURL).not.toHaveBeenCalled();
  });

  it("refreshes an expired token before a download rather than failing it", async () => {
    client.tokenStore.save({ access: "stale", refresh: "refresh-1" });
    serve((url, init) => {
      if (isRefresh(url)) return json(200, { access: "fresh" });
      return header(init, "Authorization") === "Bearer fresh"
        ? reply(200, "LOCUS pET-21a-insert 6031 bp DNA circular")
        : json(401, { detail: "Given token not valid for any token type" });
    });
    vi.spyOn(HTMLAnchorElement.prototype, "click").mockImplementation(() => {});

    await client.api.downloadUrl("/api/projects/7/genbank/", "project.gb");

    expect(refreshCalls()).toHaveLength(1);
    expect(URL.createObjectURL).toHaveBeenCalledOnce();
  });
});
