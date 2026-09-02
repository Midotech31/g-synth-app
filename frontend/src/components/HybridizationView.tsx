import type {
  HybridizationEnd,
  HybridizationResult,
  HybridizationRow,
} from "../api/client";

type Props = {
  result: HybridizationResult;
  detail?: "simple" | "detailed";
};

function SequenceLine({
  row,
  text,
  result,
}: {
  row: HybridizationRow;
  text: string;
  result: HybridizationResult;
}) {
  return (
    <span className="dx-seq">
      {[...text].map((base, index) => {
        const position = row.start + index;
        const mark = result.marks[position];
        const overhang = position < result.overlap_start || position >= result.overlap_end;
        const classes = ["dx-base"];
        if (base === " ") classes.push("dx-gap");
        else if (overhang) classes.push("hybrid-overhang-base");
        else if (mark === "×") classes.push("hybrid-mismatch-base");
        else classes.push("hybrid-paired-base");
        return (
          <span key={position} className={classes.join(" ")}>
            {base === " " ? "\u00A0" : base}
          </span>
        );
      })}
    </span>
  );
}

function EndCard({ end }: { end: HybridizationEnd }) {
  if (end.kind === "blunt") {
    return (
      <div className="hybrid-end-card blunt">
        <span className="hybrid-end-label">{end.end} end</span>
        <strong>Blunt</strong>
        <small>Both strands terminate at the same column.</small>
      </div>
    );
  }

  return (
    <div className="hybrid-end-card sticky">
      <span className="hybrid-end-label">{end.end} cohesive end</span>
      <div className="hybrid-end-title">
        <strong>{end.polarity} overhang</strong>
        <span>{end.length} nt</span>
      </div>
      <code>5′-{end.sequence}-3′</code>
      <small>Exposed on the {end.strand} input strand.</small>
    </div>
  );
}

export default function HybridizationView({ result, detail = "detailed" }: Props) {
  if (detail === "simple") {
    return (
      <div className="hybridization-view hybridization-simple">
        <div className="hybrid-simple-flow">
          <EndCard end={result.left_end} />
          <div className="hybrid-core-card">
            <span className="hybrid-end-label">Antiparallel core</span>
            <strong>{result.paired_bases} complementary base pairs</strong>
            <div className="hybrid-core-bar" aria-hidden="true">
              <span style={{ width: `${result.paired_percent}%` }} />
            </div>
            <small>
              {result.mismatches === 0
                ? "No mismatches in the overlapping region."
                : `${result.mismatches} mismatch${result.mismatches === 1 ? "" : "es"}; inspect the detailed view.`}
            </small>
          </div>
          <EndCard end={result.right_end} />
        </div>
        <p className="note hybrid-explanation">
          Ends are reported in the protruding strand&rsquo;s own 5′→3′ sequence.
          Use the nucleotide-level double-strand view to inspect every base and
          the physical 5′/3′ geometry.
        </p>
      </div>
    );
  }

  return (
    <div className="hybridization-view">
      <div className="hybrid-legend" aria-label="Hybridization legend">
        <span><i className="hybrid-key-paired" /> complementary pair</span>
        <span><i className="hybrid-key-overhang" /> unpaired overhang</span>
        <span><i className="hybrid-key-mismatch" /> mismatch</span>
      </div>

      <div
        className="duplex-scroll hybrid-duplex"
        tabIndex={0}
        aria-label="Nucleotide-level antiparallel double-strand visualization"
      >
        {result.rows.map((row) => (
          <div className="hybrid-row-block" key={row.start}>
            <div className="hybrid-strand-labels" aria-hidden="true">
              <span>Input 1</span>
              <span>Input 2 · reverse-complement orientation</span>
            </div>
            <div className="duplex-row hybrid-row">
              <span className="dx-num">{row.top_start ?? ""}</span>
              <span className="dx-end">5′</span>
              <SequenceLine row={row} text={row.top} result={result} />
              <span className="dx-end">3′</span>
            </div>
            <div className="duplex-row hybrid-row">
              <span className="dx-num" />
              <span className="dx-end" />
              <span className="dx-seq hybrid-pair-marks" aria-hidden="true">
                {[...row.marks].map((mark, index) => (
                  <span key={row.start + index}>{mark === "×" ? "×" : mark}</span>
                ))}
              </span>
              <span className="dx-end" />
            </div>
            <div className="duplex-row hybrid-row">
              <span className="dx-num">{row.bottom_start ?? ""}</span>
              <span className="dx-end">3′</span>
              <SequenceLine row={row} text={row.bottom} result={result} />
              <span className="dx-end">5′</span>
            </div>
          </div>
        ))}
      </div>

      <p className="note hybrid-explanation">
        Both inputs were entered 5′→3′. G-Synth reversed the second strand for this
        physical, antiparallel view. Vertical rules are complementary bases; × marks
        a mismatch. Amber bases are single-stranded and therefore remain exposed as
        cohesive ends rather than being counted as failed hybridization.
      </p>

      <div className="hybrid-end-grid" aria-label="Detected duplex ends">
        <EndCard end={result.left_end} />
        <EndCard end={result.right_end} />
      </div>
    </div>
  );
}
