import { useEffect, useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";

import {
  ApiError,
  api,
  type AlignResult,
  type Annotation,
  type HybridizationResult,
} from "../api/client";
import CoreWorkflowTrail from "../components/CoreWorkflowTrail";
import HybridizationView from "../components/HybridizationView";
import Icon from "../components/Icon";
import LiveStatus from "../components/LiveStatus";
import { useWorkspaceState } from "../state/WorkspaceStateContext";

const MODES = [
  { key: "global", label: "Whole of both", hint: "Two variants of one gene" },
  { key: "local", label: "Best stretch", hint: "The one region they share" },
  { key: "semi-global", label: "First in second", hint: "Align all of the first sequence within the second" },
] as const;

const SAMPLE_A = "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGAAATAATGGAGCACATATGGGA";
const SAMPLE_B = "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGTAATAATGGAGCACATATGGGA";
const STICKY_FIRST = "AATTATGCCGTAGCTAGCTA";
const STICKY_SECOND = "GGCCTAGCTAGCTACGGCAT";

type Tool = "alignment" | "hybridization";
type Handoff = {
  tool?: Tool;
  first?: string;
  second?: string;
  name?: string;
  leftEnzyme?: string;
  rightEnzyme?: string;
  orfStart?: number;
  insertAnnotations?: Annotation[];
  autoRun?: boolean;
};

const STATE_LABELS: Record<HybridizationResult["predicted_state"], string> = {
  favourable_at_temperature: "Favourable at this temperature",
  temperature_above_tm: "Temperature is above predicted Tm",
  mismatches_not_thermodynamically_scored: "Mismatches require review",
  insufficient_complementarity: "Insufficient complementarity",
  not_scored: "Not thermodynamically scored",
};

export default function Align({ initialTool = "alignment" }: { initialTool?: Tool }) {
  const location = useLocation() as { pathname: string; state?: Handoff | null };
  const navigate = useNavigate();
  const [tool, setTool, clearTool] = useWorkspaceState<Tool>(
    "align.tool",
    initialTool,
  );
  const [first, setFirst, clearFirst] = useWorkspaceState("align.first", SAMPLE_A);
  const [second, setSecond, clearSecond] = useWorkspaceState("align.second", SAMPLE_B);
  const [mode, setMode, clearMode] = useWorkspaceState<
    (typeof MODES)[number]["key"]
  >("align.mode", "global");
  const [isProtein, setIsProtein, clearIsProtein] = useWorkspaceState(
    "align.isProtein",
    false,
  );
  const [tryReverse, setTryReverse, clearTryReverse] = useWorkspaceState(
    "align.tryReverse",
    true,
  );
  const [alignment, setAlignment, clearAlignment] = useWorkspaceState<
    AlignResult | null
  >("align.result", null);
  const [hybridization, setHybridization, clearHybridization] = useWorkspaceState<
    HybridizationResult | null
  >("align.hybridization", null);
  const [temperature, setTemperature, clearTemperature] = useWorkspaceState(
    "align.temperature",
    25,
  );
  const [oligoMicromolar, setOligoMicromolar, clearOligoMicromolar] =
    useWorkspaceState("align.oligoMicromolar", 50);
  const [naMm, setNaMm, clearNaMm] = useWorkspaceState("align.naMm", 50);
  const [mgMm, setMgMm, clearMgMm] = useWorkspaceState("align.mgMm", 0);
  const [source, setSource, clearSource] = useWorkspaceState<
    Omit<Handoff, "tool" | "first" | "second" | "autoRun"> | null
  >("align.source", null);
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const [pendingHybridization, setPendingHybridization] = useState<{
    first: string;
    second: string;
  } | null>(null);

  useEffect(() => {
    const handed = location.state;
    if (!handed?.first || !handed.second) {
      setTool(initialTool);
      return;
    }
    setFirst(handed.first);
    setSecond(handed.second);
    setTool(handed.tool ?? "hybridization");
    setSource({
      name: handed.name,
      leftEnzyme: handed.leftEnzyme,
      rightEnzyme: handed.rightEnzyme,
      orfStart: handed.orfStart,
      insertAnnotations: handed.insertAnnotations,
    });
    setAlignment(null);
    setHybridization(null);
    if (handed.autoRun) {
      setPendingHybridization({ first: handed.first, second: handed.second });
    }
    navigate(location.pathname, { replace: true, state: null });
  }, [
    initialTool,
    location.pathname,
    location.state,
    navigate,
    setAlignment,
    setFirst,
    setHybridization,
    setSecond,
    setSource,
    setTool,
  ]);

  useEffect(() => {
    if (!pendingHybridization) return;
    setPendingHybridization(null);
    setBusy(true);
    setError("");
    api.hybridize({
      first: pendingHybridization.first,
      second: pendingHybridization.second,
      analysis_temperature_c: temperature,
      oligo_nM: oligoMicromolar * 1_000,
      na_mM: naMm,
      mg_mM: mgMm,
      dntp_mM: 0,
    }).then((data) => {
      setHybridization(data);
      setAlignment(null);
    }).catch((err) => {
      setError(err instanceof ApiError ? err.message : "The hybridization analysis failed.");
      setHybridization(null);
    }).finally(() => setBusy(false));
  }, [
    mgMm,
    naMm,
    oligoMicromolar,
    pendingHybridization,
    setAlignment,
    setHybridization,
    temperature,
  ]);

  const clean = (text: string) => text.replace(/[^A-Za-z]/g, "").length;

  function invalidate() {
    setAlignment(null);
    setHybridization(null);
    setError("");
  }

  async function run() {
    setBusy(true);
    setError("");
    try {
      if (tool === "alignment") {
        setAlignment(await api.align({
          first,
          second,
          mode,
          is_protein: isProtein,
          try_reverse: tryReverse,
        }));
        setHybridization(null);
      } else {
        setHybridization(
          await api.hybridize({
            first,
            second,
            analysis_temperature_c: temperature,
            oligo_nM: oligoMicromolar * 1_000,
            na_mM: naMm,
            mg_mM: mgMm,
            dntp_mM: 0,
          }),
        );
        setAlignment(null);
      }
    } catch (err) {
      setError(err instanceof ApiError ? err.message : `The ${tool} analysis failed.`);
      if (tool === "alignment") setAlignment(null);
      else setHybridization(null);
    } finally {
      setBusy(false);
    }
  }

  function clearWorkspace() {
    clearTool();
    clearFirst();
    clearSecond();
    clearMode();
    clearIsProtein();
    clearTryReverse();
    clearAlignment();
    clearHybridization();
    clearTemperature();
    clearOligoMicromolar();
    clearNaMm();
    clearMgMm();
    clearSource();
    setError("");
  }

  function useHybridizationExample() {
    setTool("hybridization");
    setFirst(STICKY_FIRST);
    setSecond(STICKY_SECOND);
    setSource(null);
    invalidate();
  }

  function continueToClone() {
    if (!hybridization || hybridization.complementarity !== "exact") return;
    navigate("/clone", {
      state: {
        preDigested: {
          top: hybridization.first,
          bottom: hybridization.second,
          leftEnzyme: source?.leftEnzyme ?? null,
          rightEnzyme: source?.rightEnzyme ?? null,
          orfStart: source?.orfStart ?? null,
          insertAnnotations: source?.insertAnnotations,
          name: source?.name,
          origin: "hybridization",
        },
      },
    });
  }

  function switchTool(next: Tool) {
    setTool(next);
    if (next === "hybridization") setIsProtein(false);
    invalidate();
    navigate(next === "hybridization" ? "/hybridize" : "/align", { replace: true });
  }

  const result = tool === "alignment" ? alignment : hybridization;
  const status = busy
    ? tool === "alignment"
      ? "Aligning…"
      : "Testing antiparallel hybridization…"
    : alignment && tool === "alignment"
      ? `${alignment.identity}% identity over ${alignment.length} ${alignment.is_protein ? "residues" : "bases"}, ${alignment.gaps} gaps.`
      : hybridization && tool === "hybridization"
        ? `${hybridization.paired_bases} paired bases, ${hybridization.mismatches} mismatches, ${hybridization.overhangs.length} cohesive ends.`
        : "";

  return (
    <>
      <LiveStatus message={status} />
      <div className="topbar">
        <div className="grow">
          <h1>{tool === "hybridization" ? "Verify hybridization" : "Align sequences"}</h1>
          <p className="sub">
            {tool === "hybridization"
              ? "Confirm antiparallel pairing, mismatches and exposed cohesive ends before restriction cloning."
              : "Compare biological similarity, identity and gaps without confusing alignment with physical annealing."}
          </p>
          {tool === "hybridization" && <CoreWorkflowTrail active="hybridization" />}
        </div>
        <button className="btn btn-outline" onClick={clearWorkspace} disabled={busy}>
          Clear
        </button>
        <button
          className="btn btn-primary"
          onClick={() => void run()}
          disabled={busy || !clean(first) || !clean(second)}
        >
          {busy && <span className="spinner" />}
          {busy ? "Analysing…" : tool === "alignment" ? "Align" : "Simulate hybridization"}
        </button>
      </div>

      <div className="content align-workspace" aria-busy={busy}>
        <div className="analysis-switch" role="tablist" aria-label="Sequence analysis">
          <button
            type="button"
            role="tab"
            aria-selected={tool === "alignment"}
            className={tool === "alignment" ? "active" : ""}
            onClick={() => switchTool("alignment")}
          >
            <Icon name="scales" size={20} />
            <span>
              <strong>Alignment</strong>
              <small>identity, homology and gaps</small>
            </span>
          </button>
          <button
            type="button"
            role="tab"
            aria-selected={tool === "hybridization"}
            className={tool === "hybridization" ? "active" : ""}
            onClick={() => switchTool("hybridization")}
          >
            <Icon name="helix" size={20} />
            <span>
              <strong>Hybridization</strong>
              <small>antiparallel pairs and sticky ends</small>
            </span>
          </button>
        </div>

        {source?.name && tool === "hybridization" && (
          <div className="notice notice-info hybrid-source">
            <strong>Transferred from {source.name}.</strong>{" "}
            The designed {source.leftEnzyme ?? "left"} and{" "}
            {source.rightEnzyme ?? "right"} ends remain attached to this duplex.
          </div>
        )}
        {error && <div className="notice notice-error" role="alert">{error}</div>}

        <div className="design-layout align-layout">
          <div className="card">
            <div className="card-head">
              <h2 style={{ flex: 1 }}>{tool === "alignment" ? "Sequences" : "Strands"}</h2>
              {tool === "hybridization" && (
                <button
                  className="btn btn-ghost"
                  type="button"
                  onClick={useHybridizationExample}
                >
                  Load sticky-end example
                </button>
              )}
            </div>
            <div className="card-body align-inputs">
              {tool === "hybridization" && (
                <div className="notice notice-info compact">
                  Enter both ordered molecules 5′→3′. The second is reversed only
                  in the physical drawing so the strands are antiparallel.
                </div>
              )}
              <div className="field">
                <label htmlFor="a">{tool === "alignment" ? "First" : "Top strand · 5′→3′"}</label>
                <textarea
                  id="a"
                  value={first}
                  onChange={(event) => {
                    setFirst(event.target.value);
                    setSource(null);
                    invalidate();
                  }}
                  rows={6}
                  className="mono"
                  style={{ fontSize: "0.76rem" }}
                  aria-describedby="a-count"
                />
                <span className="label" id="a-count">
                  {clean(first)} {isProtein ? "residues" : "nt"}
                </span>
              </div>
              <div className="field">
                <label htmlFor="b">{tool === "alignment" ? "Second" : "Partner strand · 5′→3′"}</label>
                <textarea
                  id="b"
                  value={second}
                  onChange={(event) => {
                    setSecond(event.target.value);
                    setSource(null);
                    invalidate();
                  }}
                  rows={6}
                  className="mono"
                  style={{ fontSize: "0.76rem" }}
                  aria-describedby="b-count"
                />
                <span className="label" id="b-count">
                  {clean(second)} {isProtein ? "residues" : "nt"}
                </span>
              </div>

              {tool === "alignment" ? (
                <>
                  <div className="field">
                    <span className="field-label" id="mode-label">What are you asking?</span>
                    <div className="mode-list" role="group" aria-labelledby="mode-label">
                      {MODES.map((option) => (
                        <button
                          key={option.key}
                          type="button"
                          className={mode === option.key ? "mode on" : "mode"}
                          aria-pressed={mode === option.key}
                          onClick={() => {
                            setMode(option.key);
                            invalidate();
                          }}
                        >
                          <strong>{option.label}</strong>
                          <span>{option.hint}</span>
                        </button>
                      ))}
                    </div>
                  </div>
                  <div className="checks">
                    <label>
                      <input
                        type="checkbox"
                        checked={isProtein}
                        onChange={(event) => {
                          setIsProtein(event.target.checked);
                          invalidate();
                        }}
                      />
                      These are proteins (BLOSUM62)
                    </label>
                    {!isProtein && (
                      <label>
                        <input
                          type="checkbox"
                          checked={tryReverse}
                          onChange={(event) => {
                            setTryReverse(event.target.checked);
                            invalidate();
                          }}
                        />
                        Also try the reverse complement
                      </label>
                    )}
                  </div>
                </>
              ) : (
                <details className="advanced-control" open>
                  <summary>Hybridization conditions</summary>
                  <div className="hybrid-condition-grid advanced-control-body">
                    <ConditionInput
                      id="hybrid-temperature"
                      label="Analysis temperature"
                      value={temperature}
                      unit="°C"
                      min={-20}
                      max={120}
                      step={0.5}
                      set={(value) => {
                        setTemperature(value);
                        invalidate();
                      }}
                    />
                    <ConditionInput
                      id="hybrid-concentration"
                      label="Total strand"
                      value={oligoMicromolar}
                      unit="µM"
                      min={0.001}
                      max={100_000}
                      step={0.1}
                      set={(value) => {
                        setOligoMicromolar(value);
                        invalidate();
                      }}
                    />
                    <ConditionInput
                      id="hybrid-na"
                      label="Na⁺"
                      value={naMm}
                      unit="mM"
                      min={0}
                      max={5_000}
                      step={1}
                      set={(value) => {
                        setNaMm(value);
                        invalidate();
                      }}
                    />
                    <ConditionInput
                      id="hybrid-mg"
                      label="Mg²⁺"
                      value={mgMm}
                      unit="mM"
                      min={0}
                      max={1_000}
                      step={0.1}
                      set={(value) => {
                        setMgMm(value);
                        invalidate();
                      }}
                    />
                  </div>
                </details>
              )}
            </div>
          </div>

          <div className="align-results">
            {!result ? (
              <div className="card">
                <div className="empty">
                  <Icon
                    name={tool === "alignment" ? "scales" : "helix"}
                    size={38}
                    className="glyph"
                  />
                  <strong>
                    {tool === "alignment"
                      ? "Nothing aligned yet"
                      : "No hybridization simulated yet"}
                  </strong>
                  <span>
                    {tool === "alignment"
                      ? "Paste two sequences and choose the biological question."
                      : "Add both 5′→3′ strands to inspect their pairing and exposed ends."}
                  </span>
                </div>
              </div>
            ) : tool === "alignment" && alignment ? (
              <AlignmentPanel result={alignment} />
            ) : hybridization ? (
              <HybridizationPanel
                result={hybridization}
                source={source}
                continueToClone={continueToClone}
              />
            ) : null}
          </div>
        </div>
      </div>
    </>
  );
}

type ConditionInputProps = {
  id: string;
  label: string;
  value: number;
  unit: string;
  min: number;
  max: number;
  step: number;
  set: (value: number) => void;
};

function ConditionInput({
  id,
  label,
  value,
  unit,
  min,
  max,
  step,
  set,
}: ConditionInputProps) {
  return (
    <div className="field">
      <label htmlFor={id}>{label}</label>
      <div className="input-unit">
        <input
          id={id}
          type="number"
          value={value}
          min={min}
          max={max}
          step={step}
          onChange={(event) => set(Number(event.target.value))}
        />
        <span>{unit}</span>
      </div>
    </div>
  );
}

type HybridizationPanelProps = {
  result: HybridizationResult;
  source: Omit<Handoff, "tool" | "first" | "second" | "autoRun"> | null;
  continueToClone: () => void;
};

function HybridizationPanel({
  result,
  source,
  continueToClone,
}: HybridizationPanelProps) {
  const verdict = result.complementarity === "exact"
    ? "Exact antiparallel complementarity"
    : result.complementarity === "partial"
      ? "Partial complementarity"
      : "Insufficient complementarity";

  return (
    <>
      {result.warnings.map((note) => (
        <div key={note} className="notice notice-info">{note}</div>
      ))}

      <div className={`hybrid-verdict ${result.complementarity}`}>
        <Icon
          name={result.complementarity === "exact"
            ? "check"
            : result.complementarity === "partial"
              ? "target"
              : "cross"}
          size={24}
        />
        <div>
          <strong>{verdict}</strong>
          <span>{STATE_LABELS[result.predicted_state]}</span>
        </div>
        <span className="pill">
          {result.overhangs.length
            ? `${result.overhangs.length} cohesive ends`
            : "flush duplex"}
        </span>
      </div>

      <div className="card hybrid-stats">
        <div className="card-body stat-row">
          <HybridStat label="Paired" value={result.paired_bases} unit="bp" />
          <HybridStat label="Overlap" value={result.paired_percent} unit="% paired" />
          <HybridStat label="Mismatches" value={result.mismatches} />
          <HybridStat
            label="Longest exact run"
            value={result.longest_perfect_run}
            unit="bp"
          />
          <HybridStat label="Tm" value={result.tm_c ?? "—"} unit={result.tm_c === null ? undefined : "°C"} />
          <HybridStat
            label="Tm margin"
            value={result.tm_margin_c === null
              ? "—"
              : `${result.tm_margin_c > 0 ? "+" : ""}${result.tm_margin_c}`}
            unit={result.tm_margin_c === null ? "not scored" : "°C"}
          />
        </div>
      </div>

      <div className="card">
        <div className="card-head">
          <div>
            <h2>Hybridization evidence</h2>
            <span className="label">Summary geometry and nucleotide-level double strand</span>
          </div>
        </div>
        <div className="card-body hybrid-evidence-stack">
          <section className="hybrid-evidence-section" aria-labelledby="hybrid-overview-heading">
            <div className="hybrid-evidence-heading">
              <span className="hybrid-evidence-number">1</span>
              <div>
                <h3 id="hybrid-overview-heading">Pairing overview</h3>
                <p>Read the paired core and each exposed end at a glance.</p>
              </div>
            </div>
            <HybridizationView result={result} detail="simple" />
          </section>
          <section className="hybrid-evidence-section" aria-labelledby="hybrid-double-strand-heading">
            <div className="hybrid-evidence-heading">
              <span className="hybrid-evidence-number">2</span>
              <div>
                <h3 id="hybrid-double-strand-heading">Antiparallel double-strand verification</h3>
                <p>Inspect every paired, mismatched and unpaired nucleotide before cloning.</p>
              </div>
            </div>
            <HybridizationView result={result} detail="detailed" />
          </section>
        </div>
      </div>

      <div className="card hybrid-thermo-card">
        <div className="card-head"><h2>Thermodynamic scope</h2></div>
        <div className="card-body">
          <dl className="hybrid-thermo-list">
            <div><dt>Conditions</dt><dd>{result.conditions.summary}</dd></div>
            <div>
              <dt>Model</dt>
              <dd>SantaLucia nearest-neighbour with Owczarzy salt correction</dd>
            </div>
            <div>
              <dt>ΔH</dt>
              <dd>
                {result.delta_h_kcal_mol === null
                  ? "not scored"
                  : `${result.delta_h_kcal_mol} kcal/mol`}
              </dd>
            </div>
            <div>
              <dt>ΔS</dt>
              <dd>
                {result.delta_s_cal_mol_k === null
                  ? "not scored"
                  : `${result.delta_s_cal_mol_k} cal/mol/K`}
              </dd>
            </div>
          </dl>
          <p className="note">
            Tm is intentionally withheld when the overlap contains mismatches;
            a perfect-match parameter set must not be presented as a mismatch
            prediction.
          </p>
        </div>
      </div>

      <div className="workflow-next-card">
        <div>
          <span className="label">Next core step</span>
          <strong>Simulate restriction enzyme cloning</strong>
          <small>
            {source?.leftEnzyme && source.rightEnzyme
              ? `${source.leftEnzyme} / ${source.rightEnzyme} are carried forward automatically.`
              : "Choose the vector and restriction enzymes in Cloning; the verified duplex is carried forward."}
          </small>
        </div>
        <button
          className="btn btn-primary"
          onClick={continueToClone}
          disabled={result.complementarity !== "exact"}
        >
          Simulate restriction cloning <Icon name="arrowRight" size={17} />
        </button>
      </div>
    </>
  );
}

function HybridStat({ label, value, unit }: {
  label: string;
  value: string | number;
  unit?: string;
}) {
  return (
    <div className="stat">
      <div className="k">{label}</div>
      <div className="v">{value}{unit && <small>{unit}</small>}</div>
    </div>
  );
}

function AlignmentPanel({ result }: { result: AlignResult }) {
  return (
    <>
      {result.warnings.map((note) => (
        <div key={note} className="notice notice-info">{note}</div>
      ))}

      <div className="card">
        <div className="card-body stat-row">
          <HybridStat label="Identity" value={result.identity} unit="%" />
          {result.is_protein && (
            <HybridStat label="Similarity" value={result.similarity} unit="%" />
          )}
          <HybridStat
            label="Aligned"
            value={result.length}
            unit={result.is_protein ? "aa" : "nt"}
          />
          <HybridStat label="Gaps" value={result.gaps} />
          <HybridStat label="Score" value={result.score} />
        </div>
      </div>

      <div className="card">
        <div className="card-head">
          <h2 style={{ flex: 1 }}>Alignment</h2>
          <span className="label">
            {result.start_a + 1}–{result.end_a} vs {result.start_b + 1}–{result.end_b}
          </span>
          <button
            className="btn btn-outline"
            onClick={() => void navigator.clipboard?.writeText(result.text)}
          >
            Copy
          </button>
        </div>
        <div className="card-body">
          <div className="duplex-scroll">
            {result.rows.map((row, index) => (
              <div className="duplex-row align-row" key={index}>
                <span className="dx-num">{row.top_start ?? ""}</span>
                <span className="dx-seq">{row.top}</span>
                <span className="dx-num">{row.top_end}</span>
                <span className="dx-num" />
                <span className="dx-seq dx-ticks">{row.marks}</span>
                <span className="dx-num" />
                <span className="dx-num">{row.bottom_start ?? ""}</span>
                <span className="dx-seq">{row.bottom}</span>
                <span className="dx-num">{row.bottom_end}</span>
              </div>
            ))}
          </div>
          <p className="note" style={{ marginTop: "0.7rem" }}>
            <b>|</b> identical
            {result.is_protein && <> · <b>:</b> conservative substitution</>}
            {" · "}<b>.</b> different · <b>-</b> gap.
          </p>
        </div>
      </div>
    </>
  );
}
