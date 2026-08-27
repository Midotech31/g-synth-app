/**
 * API client for the G-Synth Django backend.
 *
 * Access tokens are short-lived (30 min in production). Rather than making
 * every caller handle expiry, a 401 triggers one refresh attempt and the
 * original request is replayed. Concurrent 401s share a single refresh
 * promise, so a page issuing five requests doesn't fire five refreshes and
 * invalidate its own rotating refresh token.
 */

/**
 * Where the API lives.
 *
 * Empty in development: Vite proxies `/api` to :8000, so a relative path is
 * same-origin and no CORS is involved. In production the workspace is a
 * static site and the API is a separate service, so the two are on different
 * origins and the path has to be absolute — `VITE_API_BASE` supplies it at
 * build time. Trailing slashes are stripped so `${BASE}/api/...` never
 * doubles one.
 */
function normaliseApiBase(value: string): string {
  const base = value.trim().replace(/\/+$/, "");
  if (!base || /^https?:\/\//i.test(base)) return base;

  // Render's Blueprint service references expose `host` but not `url`.
  // Production services are HTTPS, so turn the supplied hostname into the
  // absolute origin fetch() needs. Local development keeps the empty value.
  return `https://${base}`;
}

const API_BASE = normaliseApiBase(import.meta.env.VITE_API_BASE ?? "");

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
  /** Clipped where a vector feature met the insert junction. Set only on a
   *  cloned plasmid's own features — a truncated promoter is worth seeing,
   *  not silently keeping its pre-cut length. */
  truncated?: boolean;
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

export type Diagnostic = {
  code: string;
  severity: "info" | "warning" | "error";
  message: string;
  remedy: string;
  positions: number[];
};

export type PreflightCheck = {
  code: string;
  label: string;
  status: "pass" | "review" | "block";
  detail: string;
  remedy: string;
  evidence: Record<string, unknown>;
  passed: boolean;
};

export type PreflightReport = {
  workflow: string;
  verdict: "ready" | "review" | "blocked";
  can_export: boolean;
  checks: PreflightCheck[];
  diagnostics: Diagnostic[];
};

export type Provenance = {
  schema: "gsynth.provenance/v1";
  engine_version: string;
  workflow: string;
  generated_at: string | null;
  parameters: Record<string, unknown>;
  parameters_sha256: string;
  output_sha256: string;
  vector_sha256: string | null;
  enzyme_table: { sha256: string; source: string; entries: number };
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
  preflight?: PreflightReport;
  provenance?: Provenance;
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

export type TerminalEnd = {
  side: "left" | "right";
  enzyme: string;
  overhang: string;
  /** "5'", "3'" or "blunt" — polarity follows the side, not the strand. */
  kind: string;
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
  /** The outer ends as the assembled fragments present them, not as designed. */
  terminal_ends: TerminalEnd[];
  fragments: Fragment[];
  oligos: Oligo[];
  ssd: SSDResult;
  duplex: Duplex;
  tm_conditions: { name: string; summary: string; model: string };
  warnings: string[];
  /** Empty means the oligos re-ligate to the design. Non-empty blocks ordering. */
  verification: string[];
  preflight?: PreflightReport;
  provenance?: Provenance;
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

export type JunctionView = {
  name: string;
  enzyme: string;
  overhang: string;
  kind: string;
  compatible: boolean;
  reason: string;
  left_top: string;
  left_bottom: string;
  right_top: string;
  right_bottom: string;
  joined_top: string;
  joined_bottom: string;
  joined_pairs: string;
  seam: number;
  overhang_span: [number, number];
};

export type RestrictionSite = Annotation & {
  cuts: number;
  used: boolean;
  recognition: string;
  wraps: boolean;
};

export type ValidationCheck = { check: string; passed: boolean; detail: string };

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
  junction_views: JunctionView[];
  restriction_sites: RestrictionSite[];
  validation: ValidationCheck[];
  warnings: string[];
  /** Empty means these two molecules really do join. */
  problems: string[];
  is_clonable: boolean;
  preflight?: PreflightReport;
  provenance?: Provenance;
  /** Null when the insert was supplied already cut — there was no SSD design. */
  insert: SSDResult | null;
  assembly: AssemblyResult | null;
  project_id?: number;
};

/** One PCR primer. `tail` is the 5' addition; `anneals` binds template. */
export type PcrPrimer = {
  name: string;
  sequence: string;
  tail: string;
  anneals: string;
  direction: number;
  start: number;
  end: number;
  length: number;
  anneal_length: number;
  /** Tm of the annealing part — what the annealing temperature comes from. */
  tm: number;
  /** Tm of the whole oligo, which only applies once the tail is copied. */
  tm_full: number;
  gc: number;
  enzyme: string | null;
  has_gc_clamp: boolean;
  warnings: string[];
};

export type DigestEnd = {
  sequence: string;
  strand: string;
  side: string;
  /** "5'", "3'" or "blunt" — polarity follows the side, not the strand. */
  kind: string;
};

/** The product after both ends are cut: the insert that goes into a vector. */
export type PcrDigest = {
  top: string;
  bottom: string;
  length: number;
  left_end: DigestEnd;
  right_end: DigestEnd;
  trimmed_left: number;
  trimmed_right: number;
};

export type PcrResult = {
  forward: PcrPrimer;
  reverse: PcrPrimer;
  product: string;
  product_length: number;
  amplified_region: string;
  template_start: number;
  template_end: number;
  annealing_temperature: number;
  left_enzyme: string | null;
  right_enzyme: string | null;
  insert_orf_start: number | null;
  /** Empty means the product can be cut into the insert as designed. */
  problems: string[];
  warnings: string[];
  is_clean: boolean;
  /** Null for conventional PCR, and when a problem blocks the digest. */
  digest: PcrDigest | null;
  preflight?: PreflightReport;
  provenance?: Provenance;
};

export type PcrParams = {
  template: string;
  target_start?: number;
  target_end?: number | null;
  left_enzyme?: string | null;
  right_enzyme?: string | null;
  clamp?: number;
  keep_frame?: boolean;
  start_codon_mode?: "use_site" | "keep_both";
  name?: string;
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
  preflight?: PreflightReport;
  provenance?: Provenance;
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
  /** Null when the read came as letters — "unknown", not "fine". */
  quality?: number | null;
  confident?: boolean | null;
  read_index?: number | null;
};

export type VerifyReport = {
  design_length: number;
  coverage: number;
  gaps: [number, number][];
  fully_covered: boolean;
  /** True only when the requested region is completely covered and agrees. */
  is_verified: boolean;
  verification_state?:
    | "fully_verified"
    | "differences_detected"
    | "partial_match"
    | "reads_unplaced"
    | "not_checked";
  region_start?: number;
  region_end?: number;
  preflight?: PreflightReport;
  provenance?: Provenance;
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
    /** Null when the read arrived as letters rather than as a trace. */
    mean_quality?: number | null;
    trimmed_start?: number;
    trimmed_end?: number;
  }[];
  warnings: string[];
  /** Present only on the trace endpoint. */
  traces?: TraceSummary[];
  trace_windows?: TraceWindow[];
};

export type TraceSummary = {
  name: string;
  length: number;
  mean_quality: number;
  trim_start: number;
  trim_stop: number;
  trimmed_length: number;
  high_quality_bases: number;
  sample_count: number;
  /** Enough good sequence to be worth comparing to a design at all. */
  usable: boolean;
};

/** The peaks around one difference — never a whole trace, which is megabytes. */
export type TraceWindow = {
  read: string;
  position: number;
  samples: [number, number];
  traces: Record<string, number[]>;
  bases: { index: number; base: string; quality: number; at: number }[];
  centre: number;
};

export type AlignRow = {
  top: string;
  marks: string;
  bottom: string;
  top_start: number | null;
  top_end: number;
  bottom_start: number | null;
  bottom_end: number;
};

export type AlignResult = {
  top: string;
  marks: string;
  bottom: string;
  rows: AlignRow[];
  text: string;
  score: number;
  mode: string;
  length: number;
  identity: number;
  similarity: number;
  identities: number;
  similarities: number;
  gaps: number;
  start_a: number;
  end_a: number;
  start_b: number;
  end_b: number;
  reverse_complemented: boolean;
  is_protein: boolean;
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
  /** Ligate the fragment as supplied instead of designing an insert around
   *  it. Requires `insert_reverse`: the stagger between the two strands is
   *  the overhang, so one strand alone cannot show both ends. */
  pre_digested?: boolean;
  insert_reverse?: string;
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
  /** In this lab's freezer. The pickers offer these first; the rest are
   *  still selectable, because an enzyme in your vector's polylinker
   *  should not need a code change to be usable. */
  common: boolean;
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
  data: {
    annotations?: Annotation[];
    topology?: string;
    gc_content?: number;
    construct_gc?: number;
    gc?: number;
  };
  pro