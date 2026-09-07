import type { CloneResult, VectorSpec } from "../api/client";

/** Catalogue attributes describe the source vector. Only a simulation can
 * determine which features survive the selected cuts in the recombinant. */
export default function VectorConfiguration({ spec, leftEnzyme, rightEnzyme, result, bundled, loading }: {
  spec: VectorSpec | null;
  leftEnzyme: string;
  rightEnzyme: string;
  result: CloneResult | null;
  bundled: boolean;
  loading: boolean;
}) {
  const currentResult = result?.left_enzyme === leftEnzyme && result.right_enzyme === rightEnzyme ? result : null;
  return (
    <section className="vector-brief" aria-label="Cloning configuration">
      {spec && <>
        <span className="label">Catalogue reference · {spec.name}</span>
        <div className="vector-facts">
          <span>Reference promoter: <b>{spec.promoter}</b></span>
          <span>Reference resistance: <b>{spec.resistance}</b></span>
        </div>
        <p className="note">Reference tags: {spec.tag_summary}. Their retention and translation depend on the recombinant sequence.</p>
        {!bundled && <p className="note">Imported or edited sequence: catalogue identity must be checked by simulation.</p>}
      </>}
      <div className="vector-facts" aria-live="polite">
        <span>Selected enzymes: <b>{leftEnzyme} / {rightEnzyme}</b></span>
      </div>
      <p className="note vector-note">
        {loading ? "Loading vector sequence…" : currentResult
          ? currentResult.is_clonable
            ? currentResult.reading_frame.summary
            : "The selected configuration did not pass cloning checks. Review the simulation results."
          : "Awaiting simulation for these inputs. RBS, start codon, reading frame and retained tags are assessed on the recombinant sequence."}
      </p>
      {currentResult?.is_clonable && currentResult.tags.length > 0 && <div className="vector-facts">
        {currentResult.tags.map((tag, index) => <span key={`${tag.name}-${tag.end}-${index}`}>
          {tag.end}-terminal {tag.name}: <b>{tag.present ? "detected in predicted protein" : "not detected in predicted protein"}</b>
        </span>)}
      </div>}
    </section>
  );
}
