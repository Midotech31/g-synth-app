import { useState } from "react";

import { ApiError, api, type AlignResult } from "../api/client";

/**
 * Comparing two sequences that are not assumed to be the same thing.
 *
 * Distinct from the Check page, which places a read against the construct it
 * is supposed to be. This one is for two genes from different strains, a
 * design against what a supplier returned, a protein against its homologue.
 */

const MODES = [
  { key: "global", label: "Whole of both", hint: "Two variants of one gene" },
  { key: "local", label: "Best stretch", hint: "The one region they share" },
  { key: "semi-global", label: "Shorter in longer", hint: "Where a gene sits in a plasmid" },
] as const;

const SAMPLE_A = "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGAAATAATGGAGCACATATGGGA";
const SAMPLE_B = "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGTAATAATGGAGCACATATGGGA";

export default function Align() {
  const [first, setFirst] = useState(SAMPLE_A);
  const [second, setSecond] = useState(SAMPLE_B);
  const [mode, setMode] = useState<(typeof MODES)[number]["key"]>("global");
  const [isProtein, setIsProtein] = useState(false);
  const [tryReverse, setTryReverse] = useState(true);
  const [result, setResult] = useState<AlignResult | null>(null);
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);

  async function run() {
    setBusy(true);
    setError("");
    try {
      setResult(await api.align({
        first, second, mode, is_protein: isProtein, try_reverse: tryReverse,
      }));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "The alignment failed.");
      setResult(null);
    } finally {
      setBusy(false);
    }
  }

  const clean = (text: string) => text.replace(/[^A-Za-z]/g, "").length;

  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>Compare two sequences</h1>
          <p className="sub">
            Two strains, a design against what came back, a protein against
            its homologue. Gaps are scored the way biology makes them.
          </p>
        </div>
        <button
          className="btn btn-primary"
          onClick={() => void run()}
          disabled={busy || !clean(first) || !clean(second)}
        >
          {busy && <span className="spinner" />}
          {busy ? "Aligning…" : "Align"}
        </button>
      </div>

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
        {error && <div className="notice notice-error">{error}</div>}

        <div className="design-layout">
          <div className="card">
            <div className="card-head"><h2>Sequences</h2></div>
            <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.85rem" }}>
              <div className="field">
                <label htmlFor="a">First</label>
                <textarea id="a" value={first} onChange={(e) => { setFirst(e.target.value); setResult(null); }}
                          rows={6} className="mono" style={{ fontSize: "0.76rem" }} />
                <span className="label">{clean(first)} {isProtein ? "residues" : "nt"}</span>
              </div>
              <div className="field">
                <label htmlFor="b">Second</label>
                <textarea id="b" value={second} onChange={(e) => { setSecond(e.target.value); setResult(null); }}
                          rows={6} className="mono" style={{ fontSize: "0.76rem" }} />
                <span className="label">{clean(second)} {isProtein ? "residues" : "nt"}</span>
              </div>

              <div className="field">
                <label>What are you asking?</label>
                <div className="mode-list">
                  {MODES.map((option) => (
                    <button
                      key={option.key}
                      type="button"
                      className={mode === option.key ? "mode on" : "mode"}
                      onClick={() => { setMode(option.key); setResult(null); }}
                    >
                      <strong>{option.label}</strong>
                      <span>{option.hint}</span>
                    </button>
                  ))}
                </div>
              </div>

              <div className="checks">
                <label>
                  <input type="checkbox" checked={isProtein}
                         onChange={(e) => { setIsProtein(e.target.checked); setResult(null); }} />
                  These are proteins (BLOSUM62)
                </label>
                {!isProtein && (
                  <label>
                    <input type="checkbox" checked={tryReverse}
                           onChange={(e) => { setTryReverse(e.target.checked); setResult(null); }} />
                    Also try the reverse complement
                  </label>
                )}
              </div>
            </div>
          </div>

          <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
            {!result ? (
              <div className="card">
                <div className="empty">
                  <div className="glyph">⚖️</div>
                  <strong>Nothing aligned yet</strong>
                  <span>Paste two sequences and choose what you are asking.</span>
                </div>
              </div>
            ) : (
              <>
                {result.warnings.map((note) => (
                  <div key={note} className="notice notice-info">{note}</div>
                ))}

                <div className="card">
                  <div className="card-body stat-row">
                    <div className="stat">
                      <div className="k">Identity</div>
                      <div className="v">{result.identity}<small>%</small></div>
                    </div>
                    {result.is_protein && (
                      <div className="stat">
                        <div className="k">Similarity</div>
                        <div className="v">{result.similarity}<small>%</small></div>
                      </div>
                    )}
                    <div className="stat">
                      <div className="k">Aligned</div>
                      <div className="v">{result.length}<small>{result.is_protein ? "aa" : "nt"}</small></div>
                    </div>
                    <div className="stat">
                      <div className="k">Gaps</div>
                      <div className="v">{result.gaps}</div>
                    </div>
                    <div className="stat">
                      <div className="k">Score</div>
                      <div className="v">{result.score}</div>
                    </div>
                  </div>
                </div>

                <div className="card">
                  <div className="card-head">
                    <h2 style={{ flex: 1 }}>Alignment</h2>
                    <span className="label">
                      {result.start_a + 1}–{result.end_a} vs {result.start_b + 1}–{result.end_b}
                    </span>
                    <button className="btn btn-outline"
                            onClick={() => void navigator.clipboard?.writeText(result.text)}>
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
                      {" "}· <b>.</b> different · <b>-</b> gap. One long gap is one
                      event, so opening one is expensive and extending it is cheap.
                    </p>
                  </div>
                </div>
              </>
            )}
          </div>
        </div>
      </div>
    </>
  );
}
