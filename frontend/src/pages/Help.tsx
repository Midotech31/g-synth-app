import { Link } from "react-router-dom";

import Icon, { type IconName } from "../components/Icon";

/** Concise bench-facing reference for each workflow and its evidence. */

type Section = {
  to: string;
  name: string;
  icon: IconName;
  summary: string;
  points: string[];
};

const SECTIONS: Section[] = [
  {
    to: "/optimise", name: "Optimise", icon: "helix",
    summary: "Rewrite a gene for the organism that will express it. The protein never changes.",
    points: [
      "Codon usage is matched to the host, and the translated protein is checked to be identical before and after.",
      "List the enzymes the construct will be cloned with here — a gene that translates beautifully but carries an internal NdeI site cannot be cloned NdeI/XhoI, and this is where that gets caught.",
      "“Send to Design” hands the optimised sequence straight to the next stage.",
    ],
  },
  {
    to: "/design", name: "Design", icon: "helix",
    summary: "Use Small Sequence Design (SSD) for one oligo pair or Extended Sequence Design (ESD) for tiled pairs.",
    points: [
      "Choose the tag, linker, protease site, and the enzyme pair the cassette will carry at each end.",
      "The page hands back the oligos to order, a bench protocol, and the hybridisation view — both strands drawn aligned, with the overhangs showing.",
      "Nothing can be downloaded until re-ligating the fragments in silico reproduces the construct base for base, on both strands, with the sticky ends the chosen enzymes actually leave.",
    ],
  },
  {
    to: "/hybridize", name: "Hybridization", icon: "check",
    summary: "Verify the designed duplex before it enters a cloning simulation.",
    points: [
      "Both ordered molecules remain entered 5′→3′; G-Synth reverses the partner only in the physical drawing so the duplex is antiparallel.",
      "The result always includes both the compact pairing overview and a nucleotide-level double-strand view, with mismatches and exposed 5′/3′ cohesive ends distinguished.",
      "A design transfer runs the hybridization automatically and carries its verified strands and enzyme pair into restriction cloning in one click.",
    ],
  },
  {
    to: "/clone", name: "Restriction cloning", icon: "plate",
    summary: "Simulate restriction digestion, end compatibility and ligation into a vector.",
    points: [
      "pET-21a(+) and pET-21(+) ship with their sequences; any other backbone is imported — SnapGene .dna, GenBank or FASTA — and checked against the catalogue entry, so pasting the wrong one is caught rather than cloned into.",
      "Click a feature or a restriction site on the plasmid map to see exactly what it is and the bases it spans.",
      "Restriction sites default to useful single-cutters plus the cloning pair. Enable multi-cutters to display every occurrence in the full catalogue; origin-crossing sites are split across both map ends rather than hidden.",
      "Each seam is drawn as the two ends that made it, so “the overhangs match” can be checked instead of believed.",
      "The diagnostic-gel simulation digests the final recombinant sequence with the cloning pair and plots every calculated fragment beside a selectable ladder.",
    ],
  },
  {
    to: "/verify", name: "Validate", icon: "microscope",
    summary: "Close the loop: ligation amounts, sequencing primers, and what the reads say.",
    points: [
      "Ligation is worked out in fmol, not nanograms — at equal mass a 5.4 kb vector outnumbers a 150 bp insert thirty-six to one.",
      "Sequencing primers sit back from the insert rather than at it.",
      "Upload the ABIF (.ab1) or SCF traces the facility sends back and compare them to the design directly — differences below Q20 confidence are marked unconfident rather than reported as mutations.",
      "In Projects, the Annotated view expands a locus into coordinates, overlapping feature tracks and translation. Add, rename, edit or delete any feature; an exact-motif scan proposes common elements for review before saving them to GenBank.",
    ],
  },
  {
    to: "/pcr", name: "Primer design & PCR", icon: "target",
    summary: "A supporting workflow for conventional or cloning PCR and primer annealing.",
    points: [
      "The primer–template view aligns both oligos base by base. A cloning primer’s 5′ clamp and restriction site are shown unpaired in cycle 1; only its 3′ region hybridizes and sets the annealing temperature.",
      "The product is assembled from the primer and template sequences, then cut in silico with the selected enzymes. Internal or junction-created sites block the design.",
      "The predicted gel shows the calculated amplicon beside an expected no-template control and a selectable generic DNA ladder. It is a size prediction, not experimental evidence.",
    ],
  },
  {
    to: "/align", name: "Sequence alignment", icon: "scales",
    summary: "Compare sequence similarity as a supporting analysis, separately from physical hybridization.",
    points: [
      "Two strains, a design against what a supplier returned, a protein against its homologue — anything the rest of the workflow does not already cover.",
      "Whole-of-both, best-stretch, and shorter-in-longer alignment modes, on nucleotide or protein sequences.",
      "Alignment remains available in the same analysis workspace but has its own route and language, so sequence similarity is never presented as annealing evidence.",
    ],
  },
];

export default function Help() {
  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>Help</h1>
          <p className="sub">What each page does, the method behind Design, and what the messages on screen mean.</p>
        </div>
        <Link to="/" className="btn btn-outline">
          <Icon name="arrowLeft" size={16} /> Home
        </Link>
      </div>

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
        <div className="card">
          <div className="card-head"><h2>The workflow</h2></div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "1.3rem" }}>
            {SECTIONS.map((s) => (
              <div key={s.to} className="help-section">
                <Link to={s.to} className="help-stage-head">
                  <Icon name={s.icon} size={22} />
                  <h3>{s.name}</h3>
                </Link>
                <p className="note" style={{ margin: "0.35rem 0 0.55rem" }}>{s.summary}</p>
                <ul className="help-points">
                  {s.points.map((point) => (
                    <li key={point}>{point}</li>
                  ))}
                </ul>
              </div>
            ))}
          </div>
        </div>

        <div className="card">
          <div className="card-head"><h2>SSD and ESD, in short</h2></div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.8rem" }}>
            <p className="note">
              Small Sequence Design (SSD) emits one complementary synthesis pair. When the construct is
              longer, Extended Sequence Design (ESD) cuts it in silico into short fragments and orders each as two oligos &mdash; a forward strand and a reverse
              strand &mdash; which anneal into a short double-stranded piece with a single-stranded
              overhang left sticking out at each end.
            </p>
            <p className="note">
              Every overhang in the assembly is chosen so it only pairs with the one overhang it is meant
              to join &mdash; never with any other junction in the same design, and never with itself.
              Mixed together and ligated in one tube, the fragments can therefore only come together in
              the order they were designed in; there is no other way for the sticky ends to fit.
            </p>
            <p className="note">
              The two outer ends of the whole cassette are the sticky ends the chosen cloning enzymes
              actually leave, so the finished piece drops straight into a vector cut with the same pair.
              Before a design can be downloaded, the software reassembles the fragments itself, base by
              base on both strands, and checks the result against the construct and against what the
              enzymes should leave &mdash; a check run on the molecule, not on what the design intended it
              to be.
            </p>
          </div>
        </div>

        <div className="card">
          <div className="card-head"><h2>Reading a result</h2></div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.8rem" }}>
            <p className="note">
              <strong>Something that stops you</strong> &mdash; a leftover restriction site, a construct
              that will not clone, a verification that fails &mdash; is written out in full in place of
              the result: the message on screen is the exact reason, not a generic failure.
            </p>
            <p className="note">
              <strong>Something worth knowing but not blocking</strong> &mdash; a rare codon left in to
              satisfy a GC window, for instance &mdash; shows up as a note alongside the result instead. It
              costs a little translation speed, not the strategy.
            </p>
            <p className="note">
              On the Validate page, a difference between a sequencing read and the design is marked
              <strong> unconfident</strong> when the trace&rsquo;s quality at that position falls below the
              threshold a base call is trusted at &mdash; the same letters can mean a real change or a bad
              peak, and only the trace tells them apart.
            </p>
            <p className="note">
              Sequencing has five explicit outcomes: <strong>fully verified</strong>, <strong>differences detected</strong>,
              <strong> partial match</strong>, <strong>reads unplaced</strong>, and <strong>not checked</strong>.
              Only the first means the entire requested region is covered and agrees with the design.
            </p>
          </div>
        </div>

        <div className="card">
          <div className="card-head"><h2>Preflight &amp; reproducibility</h2></div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.8rem" }}>
            <p className="note">
              Every workflow uses the same three release words. <strong>Ready</strong> has no unresolved checks.
              <strong> Review required</strong> is scientifically usable but needs a documented choice.
              <strong> Blocked</strong> disables release because the proposed molecule or evidence fails a required check.
            </p>
            <p className="note">
              Expand any preflight row to see its stable code, evidence and next action. Codes stay fixed even when
              explanatory wording improves, so saved projects and audit records remain interpretable.
            </p>
            <p className="note">
              Saved projects carry the engine version and SHA-256 checksums of the parameters, output sequence,
              vector and enzyme table. GenBank and FASTA exports include the output identity, allowing a later file
              to be traced back to the exact calculated molecule without duplicating raw sequences in the manifest.
            </p>
          </div>
        </div>

        <div className="card">
          <div className="card-head"><h2>Glossary</h2></div>
          <div className="card-body">
            <dl className="help-glossary">
              <div><dt>Recognition site</dt><dd>The DNA motif an enzyme recognises; it is not necessarily the bases retained at the cut end.</dd></div>
              <div><dt>Overhang</dt><dd>The single-stranded bases exposed after cleavage. Sequence, strand and side determine compatibility and polarity.</dd></div>
              <div><dt>Orthogonal junction</dt><dd>An assembly overhang that does not pair with itself, another junction, or another junction&rsquo;s reverse complement.</dd></div>
              <div><dt>Reading frame</dt><dd>The grouping of coding bases into triplets from the selected translation start through both junctions.</dd></div>
              <div><dt>Coverage gap</dt><dd>A requested interval supported by no placed sequencing read; agreement elsewhere cannot verify it.</dd></div>
              <div><dt>Provenance</dt><dd>The versioned record and checksums that identify the inputs, reference table and exact calculated output.</dd></div>
            </dl>
          </div>
        </div>

        <div className="card">
          <div className="card-head"><h2>Limits of the prediction</h2></div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.8rem" }}>
            <p className="note">
              G-Synth checks sequence logic. It does not predict star activity, methylation sensitivity, partial
              digestion, enzyme-lot performance, secondary structure, transformation efficiency, toxicity,
              expression yield or protein solubility. Confirm buffers, units and incubation conditions against the
              current manufacturer datasheets and retain appropriate experimental controls.
            </p>
            <p className="note">
              A calculated Tm, compatible end or translated protein is a design prediction, not an experimental
              result. Sequence-verify the complete requested region before expression and adjudicate every confident
              difference from the chromatogram.
            </p>
          </div>
        </div>

        <div className="card">
          <div className="card-head"><h2>About &amp; contact</h2></div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.25rem" }}>
            <strong>Designed by Prof. Merzoug Mohamed</strong>
            <span className="note">Full Professor</span>
            <span className="note">Genomics Technology Platform</span>
            <span className="note">Higher School of Biological Sciences of Oran</span>
            <span className="note" style={{ marginTop: "0.35rem" }}>
              <a href="mailto:mohamed.merzoug.essbo@gmail.com">mohamed.merzoug.essbo@gmail.com</a>
            </span>
          </div>
        </div>
      </div>
    </>
  );
}
