import { useMemo, useState } from "react";

import { ApiError, api, type ORFRecord, type SequenceAnalysis } from "../api/client";
import Icon from "../components/Icon";
import LiveStatus from "../components/LiveStatus";
import { useWorkspaceState } from "../state/WorkspaceStateContext";

const SAMPLE = "CCCATGAAACCCGGGTTTAAACCCGGGTAAGGG";
const FRAMES = [1, 2, 3, -1, -2, -3];

function downloadText(filename: string, text: string, type: string) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  anchor.style.display = "none";
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  window.setTimeout(() => URL.revokeObjectURL(url), 0);
}

function orfsAsFasta(orfs: ORFRecord[]): string {
  return orfs.map((orf) =>
    `>ORF_${orf.index} frame=${orf.frame > 0 ? "+" : ""}${orf.frame} `
    + `range=${orf.start + 1}-${orf.end} aa=${orf.amino_acids}\n${orf.protein}`,
  ).join("\n");
}

function orfsAsCsv(orfs: ORFRecord[]): string {
  const rows = ["index,strand,frame,start_1based,end_1based,length_nt,amino_acids,stop_codon,dna,protein"];
  orfs.forEach((orf) => rows.push([
    orf.index, orf.strand, orf.frame, orf.start + 1, orf.end, orf.length_nt,
    orf.amino_acids, orf.stop_codon, orf.dna, orf.protein,
  ].join(",")));
  return `${rows.join("\n")}\n`;
}

export default function Translate() {
  const [sequence, setSequence, clearSequence] = useWorkspaceState("translate.sequence", SAMPLE);
  const [minimumCodons, setMinimumCodons, clearMinimum] = useWorkspaceState("translate.minimumCodons", 5);
  const [selectedFrame, setSelectedFrame, clearFrame] = useWorkspaceState("translate.frame", 1);
  const [startAtAtg, setStartAtAtg, clearStartAtAtg] = useWorkspaceState("translate.startAtAtg", true);
  const [result, setResult, clearResult] = useWorkspaceState<SequenceAnalysis | null>("translate.result", null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");
  const cleanLength = useMemo(() => sequence.replace(/[^A-Za-z]/g, "").length, [sequence]);

  async function run() {
    setBusy(true);
    setError("");
    try {
      setResult(await api.analyse({ sequence, minimum_codons: minimumCodons }));
    } catch (err) {
      setResult(null);
      setError(err instanceof ApiError ? err.message : "The sequence could not be analysed.");
    } finally {
      setBusy(false);
    }
  }

  function clearWorkspace() {
    clearSequence(); clearMinimum(); clearFrame(); clearStartAtAtg(); clearResult(); setError("");
  }

  const frame = result?.frames.find((item) => item.frame === selectedFrame) ?? null;
  const protein = frame ? (startAtAtg ? frame.protein_from_first_atg : frame.protein) : "";

  return (
    <>
      <LiveStatus message={busy ? "Analysing six frames…" : result ? `${result.orfs.length} complete ORFs found.` : ""} />
      <div className="topbar">
        <div className="grow">
          <h1>Translate &amp; find ORFs</h1>
          <p className="sub">Six reading frames, both strands, exact coordinates and exportable candidates.</p>
        </div>
        <button className="btn btn-outline" onClick={clearWorkspace} disabled={busy}>Clear</button>
        <button className="btn btn-primary" onClick={() => void run()} disabled={busy || !cleanLength}>
          {busy && <span className="spinner" />}{busy ? "Analysing…" : "Analyse"}
        </button>
      </div>

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }} aria-busy={busy}>
        {error && <div className="notice notice-error" role="alert">{error}</div>}
        <div className="design-layout">
          <div className="card">
            <div className="card-head"><h2>DNA</h2></div>
            <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.9rem" }}>
              <div className="field">
                <label htmlFor="translation-sequence">Sequence (A/C/G/T)</label>
                <textarea id="translation-sequence" className="mono" rows={9} value={sequence}
                  onChange={(event) => { setSequence(event.target.value); setResult(null); }}
                  aria-describedby="translation-count" />
                <span id="translation-count" className="label">{cleanLength} nt</span>
              </div>
              <div className="field">
                <label htmlFor="minimum-codons">Minimum translated amino acids</label>
                <input id="minimum-codons" type="number" min={1} max={100000} value={minimumCodons}
                  onChange={(event) => { setMinimumCodons(Number(event.target.value)); setResult(null); }} />
                <span className="note">Complete ORFs begin with ATG and end at TAA, TAG or TGA.</span>
              </div>
            </div>
          </div>

          {!result ? (
            <div className="card"><div className="empty">
              <Icon name="helix" size={38} className="glyph" />
              <strong>No translation yet</strong><span>Analyse DNA to inspect all six reading frames.</span>
            </div></div>
          ) : (
            <div style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
              <div className="card"><div className="card-body stat-row">
                <div className="stat"><div className="k">Length</div><div className="v">{result.length}<small>nt</small></div></div>
                <div className="stat"><div className="k">GC</div><div className="v">{result.gc}<small>%</small></div></div>
                <div className="stat"><div className="k">Complete ORFs</div><div className="v">{result.orfs.length}</div></div>
              </div></div>
              <div className="card">
                <div className="card-head"><h2 style={{ flex: 1 }}>Translation</h2>
                  <button className="btn btn-outline" disabled={!protein}
                    onClick={() => void navigator.clipboard?.writeText(protein)}>Copy protein</button>
                </div>
                <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.85rem" }}>
                  <div className="seg-toggle" role="group" aria-label="Reading frame">
                    {FRAMES.map((value) => <button key={value} className={selectedFrame === value ? "on" : ""}
                      onClick={() => setSelectedFrame(value)} aria-pressed={selectedFrame === value}>
                      {value > 0 ? `+${value}` : value}
                    </button>)}
                  </div>
                  <div className="checks"><label><input type="checkbox" checked={startAtAtg}
                    onChange={(event) => setStartAtAtg(event.target.checked)} />Start at the first in-frame ATG</label></div>
                  {startAtAtg && frame?.first_atg === null && <div className="notice notice-info compact">No in-frame ATG in this frame.</div>}
                  <div className="seq-block" aria-label={`Protein translation frame ${selectedFrame}`}>{protein || "—"}</div>
                </div>
              </div>
            </div>
          )}
        </div>

        {result && <div className="card">
          <div className="card-head"><h2 style={{ flex: 1 }}>Complete ORFs</h2>
            <button className="btn btn-outline" disabled={!result.orfs.length}
              onClick={() => downloadText("gsynth_orfs.fasta", orfsAsFasta(result.orfs), "text/plain")}>FASTA</button>
            <button className="btn btn-outline" disabled={!result.orfs.length}
              onClick={() => downloadText("gsynth_orfs.csv", orfsAsCsv(result.orfs), "text/csv")}>CSV</button>
          </div>
          {result.orfs.length ? <div className="table-scroll"><table className="data">
            <thead><tr><th>ORF</th><th>Strand</th><th>Frame</th><th>Range</th><th>Length</th><th>Stop</th><th>Protein</th></tr></thead>
            <tbody>{result.orfs.map((orf) => <tr key={orf.index}>
              <td className="num">{orf.index}</td><td>{orf.strand}</td>
              <td className="num">{orf.frame > 0 ? `+${orf.frame}` : orf.frame}</td>
              <td className="num">{orf.start + 1}–{orf.end}</td><td className="num">{orf.amino_acids} aa</td>
              <td className="mono">{orf.stop_codon}</td><td className="mono seq-cell" title={orf.protein}>{orf.protein}</td>
            </tr>)}</tbody>
          </table></div> : <div className="card-body"><p className="note">No complete ORF meets the selected minimum length.</p></div>}
        </div>}
      </div>
    </>
  );
}
