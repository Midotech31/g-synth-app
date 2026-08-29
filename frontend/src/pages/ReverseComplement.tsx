import { useMemo, useState } from "react";

import { ApiError, api, type SequenceAnalysis } from "../api/client";
import Icon from "../components/Icon";
import LiveStatus from "../components/LiveStatus";
import { useWorkspaceState } from "../state/WorkspaceStateContext";

export default function ReverseComplement() {
  const [sequence, setSequence, clearSequence] = useWorkspaceState("reverse.sequence", "");
  const [result, setResult, clearResult] = useWorkspaceState<SequenceAnalysis | null>("reverse.result", null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");
  const [copied, setCopied] = useState(false);
  const cleanLength = useMemo(() => sequence.replace(/[^A-Za-z]/g, "").length, [sequence]);

  async function run() {
    setBusy(true); setError(""); setCopied(false);
    try { setResult(await api.analyse({ sequence, minimum_codons: 1 })); }
    catch (err) {
      setResult(null);
      setError(err instanceof ApiError ? err.message : "The reverse complement could not be generated.");
    } finally { setBusy(false); }
  }

  function clearWorkspace() { clearSequence(); clearResult(); setError(""); setCopied(false); }
  async function copy() {
    if (!result) return;
    await navigator.clipboard.writeText(result.reverse_complement);
    setCopied(true);
  }

  return <>
    <LiveStatus message={busy ? "Calculating reverse complement…" : copied ? "Reverse complement copied." : ""} />
    <div className="topbar"><div className="grow">
      <h1>Reverse complement</h1><p className="sub">Transform DNA without changing or silently accepting ambiguous bases.</p>
    </div>
      <button className="btn btn-outline" onClick={clearWorkspace} disabled={busy}>Clear</button>
      <button className="btn btn-primary" onClick={() => void run()} disabled={busy || !cleanLength}>
        {busy && <span className="spinner" />}{busy ? "Calculating…" : "Transform"}
      </button>
    </div>
    <div className="content design-layout" aria-busy={busy}>
      <div className="card"><div className="card-head"><h2>Input DNA</h2></div><div className="card-body">
        <div className="field"><label htmlFor="reverse-input">Sequence (A/C/G/T)</label>
          <textarea id="reverse-input" className="mono" rows={12} value={sequence}
            onChange={(event) => { setSequence(event.target.value); setResult(null); setCopied(false); }}
            aria-describedby="reverse-count" />
          <span id="reverse-count" className="label">{cleanLength} nt</span>
        </div>
      </div></div>
      <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
        {error && <div className="notice notice-error" role="alert">{error}</div>}
        {!result ? <div className="card"><div className="empty">
          <Icon name="helix" size={38} className="glyph" /><strong>No result yet</strong>
          <span>Paste DNA and transform it.</span>
        </div></div> : <div className="card">
          <div className="card-head"><h2 style={{ flex: 1 }}>Reverse complement · 5&prime;→3&prime;</h2>
            <button className="btn btn-outline" onClick={() => void copy()}>{copied ? "Copied" : "Copy"}</button>
          </div>
          <div className="card-body"><div className="seq-block">{result.reverse_complement}</div>
            <p className="note" style={{ marginTop: "0.7rem" }}>{result.length} nt · {result.gc}% GC · validated A/C/G/T only</p>
          </div>
        </div>}
      </div>
    </div>
  </>;
}
