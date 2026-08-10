import { Link } from "react-router-dom";

import Icon, { type IconName } from "../components/Icon";

/**
 * What each page does, the method behind Design, and what the messages on
 * screen mean. Reference, not a tutorial — the wording below is the same
 * claim the README and CLAUDE.md make, said for someone at the bench rather
 * than someone reading the source.
 */

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
    summary: "Build the cassette and cut it into oligo pairs for Merzoug assembly — see below.",
    points: [
      "Choose the tag, linker, protease site, and the enzyme pair the cassette will carry at each end.",
      "The page hands back the oligos to order, a bench protocol, and the hybridisation view — both strands drawn aligned, with the overhangs showing.",
      "Nothing can be downloaded until re-ligating the fragments in silico reproduces the construct base for base, on both strands, with the sticky ends the chosen enzymes actually leave.",
    ],
  },
  {
    to: "/clone", name: "Clone", icon: "plate",
    summary: "Cut a vector and put the construct in.",
    points: [
      "pET-21a(+) and pET-21(+) ship with their sequences; any other backbone is imported — SnapGene .dna, GenBank or FASTA — and checked against the catalogue entry, so pasting the wrong one is caught rather than cloned into.",
      "Click a feature or a restriction site on the plasmid map to see exactly what it is and the bases it spans.",
      "Each seam is drawn as the two ends that made it, so “the overhangs match” can be checked instead of believed.",
    ],
  },
  {
    to: "/verify", name: "Check", icon: "microscope",
    summary: "Close the loop: ligation amounts, sequencing primers, and what the reads say.",
    points: [
      "Ligation is worked out in fmol, not nanograms — at equal mass a 5.4 kb vector outnumbers a 150 bp insert thirty-six to one.",
      "Sequencing primers sit back from the insert rather than at it.",
      "Upload the .ab1 traces the facility sends back and compare them to the design directly — differences below Q20 confidence are marked unconfident rather than reported as mutations.",
    ],
  },
  {
    to: "/align", name: "Compare", icon: "scales",
    summary: "Align two sequences that are not assumed to be the same thing.",
    points: [
      "Two strains, a design against what a supplier returned, a protein against its homologue — anything the rest of the workflow does not already cover.",
      "Whole-of-both, best-stretch, and shorter-in-longer alignment modes, on nucleotide or protein sequences.",
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
          <div className="card-head"><h2>The Merzoug method, in short</h2></div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.8rem" }}>
            <p className="note">
              Design does not build the cassette with PCR. It cuts the finished construct, on paper, into
              short fragments and orders each one as two oligos &mdash; a forward strand and a reverse
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
              On the Check page, a difference between a sequencing read and the design is marked
              <strong> unconfident</strong> when the trace&rsquo;s quality at that position falls below the
              threshold a base call is trusted at &mdash; the same letters can mean a real change or a bad
              peak, and only the trace tells them apart.
            </p>
          </div>
        </div>
      </div>
    </>
  );
}
