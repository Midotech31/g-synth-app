import { Fragment, useEffect, useState } from "react";
import { Link } from "react-router-dom";

import { ApiError, api, type ProjectSummary } from "../api/client";
import { useAuth } from "../auth/AuthContext";
import Icon, { type IconName } from "../components/Icon";
import { Logo } from "../components/Logo";

/** One pipeline stage: its own icon (reused from that page's own empty
 *  state, so a returning user already knows what it means) and the exact
 *  claim the README makes for it — this page draws on that wording rather
 *  than inventing separate marketing copy that could drift out of truth
 *  with what the tool actually does. */
const STAGES: { to: string; name: string; icon: IconName; blurb: string }[] = [
  {
    to: "/optimise", name: "Optimise", icon: "helix",
    blurb: "Rewrite a gene for the host that will express it. The protein never changes.",
  },
  {
    to: "/design", name: "Design", icon: "helix",
    blurb: "SSD for one synthesis pair; ESD for tiled oligo pairs — exact sticky ends, no PCR.",
  },
  {
    to: "/hybridize", name: "Hybridize", icon: "check",
    blurb: "Verify the intended antiparallel duplex and inspect every exposed cohesive end.",
  },
  {
    to: "/clone", name: "Restriction clone", icon: "plate",
    blurb: "Digest vector and insert, verify both junctions, then ligate the compatible product.",
  },
  {
    to: "/verify", name: "Validate", icon: "microscope",
    blurb: "Ligation amounts, sequencing primers, and what ABIF/SCF traces say when they come back.",
  },
];

const SUPPORTING_TOOLS: { to: string; name: string; icon: IconName; blurb: string }[] = [
  {
    to: "/pcr", name: "Primer design & PCR", icon: "microscope",
    blurb: "Design conventional or cloning primers, inspect annealing, and predict an amplicon gel.",
  },
  {
    to: "/align", name: "Sequence alignment", icon: "scales",
    blurb: "Compare nucleotide or protein sequences by identity, similarity and gaps.",
  },
  {
    to: "/learn", name: "Learn", icon: "book",
    blurb: "Review concise molecular-biology references behind the workflows.",
  },
];

/** Concrete safeguards expressed in language a bench scientist can act on. */
const GUARANTEES = [
  {
    title: "Nothing ships unverified",
    body: "A design re-ligates in silico before it can be downloaded — on both "
        + "strands, sticky ends included. A tool that skips this step can hand "
        + "back a construct it never actually simulated.",
  },
  {
    title: "A trace is read, not just displayed",
    body: "Upload the ABIF or SCF trace the facility sends back. Differences below Q20 are "
        + "marked as unconfident rather than reported as mutations — the same "
        + "letters can mean a real change or a bad peak, and only the trace "
        + "tells them apart.",
  },
  {
    title: "Amounts are molar, not by mass",
    body: "At equal mass a 5.4 kb vector outnumbers a 150 bp insert thirty-six "
        + "to one. The ligation table is worked out in fmol so the ratio on the "
        + "bench is the ratio that was meant.",
  },
];

function formatDate(iso: string): string {
  return new Date(iso).toLocaleDateString(undefined, { year: "numeric", month: "short", day: "numeric" });
}

export default function Home() {
  const { user } = useAuth();
  const [projects, setProjects] = useState<ProjectSummary[] | null>(null);
  const [error, setError] = useState("");

  useEffect(() => {
    let cancelled = false;
    api.listProjects()
      .then((page) => { if (!cancelled) setProjects(page.results); })
      .catch((err) => {
        if (!cancelled) setError(err instanceof ApiError ? err.message : "");
      });
    return () => { cancelled = true; };
  }, []);

  const firstName = (user?.name || "").split(" ").find((w) => !/^(Dr|Pr|Mr|Mrs|Ms)\.?$/i.test(w));
  const recent = (projects ?? []).slice(0, 3);

  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>{firstName ? `Welcome back, ${firstName}` : "Welcome"}</h1>
          <p className="sub">
            Gene → oligos to order, a plasmid to build, and the checks that say it is the
            construct you designed.
          </p>
        </div>
        <Link to="/help" className="btn btn-outline">
          <Icon name="target" size={16} /> Help
        </Link>
      </div>

      <div className="content" style={{ display: "flex", flexDirection: "column", gap: "1.1rem" }}>
        <div className="card home-hero">
          <div className="card-body home-hero-body">
            <Logo size={64} />
            <div className="home-hero-copy">
              <span className="label">Design → hybridize → restriction clone → validate</span>
              <h2>One auditable record from sequence intent to verified construct.</h2>
              <p>
                Create synthesis-ready oligos, simulate the cloning geometry and bring the
                resulting sequencing evidence back to the exact design that produced it.
              </p>
              <div className="home-hero-actions">
                <Link to="/design" className="btn btn-primary">
                  <Icon name="helix" size={17} /> Design a construct
                </Link>
                <Link to="/verify" className="btn btn-outline">
                  <Icon name="microscope" size={17} /> Validate sequencing
                </Link>
              </div>
            </div>
          </div>
        </div>

        <div className="card">
          <div className="card-head"><h2>The workflow</h2></div>
          <div className="pipeline">
            {STAGES.map((stage, i) => (
              <Fragment key={stage.to}>
                <Link to={stage.to} className="pipeline-stage">
                  <Icon name={stage.icon} size={26} />
                  <strong>{stage.name}</strong>
                  <span>{stage.blurb}</span>
                </Link>
                {i < STAGES.length - 1 && (
                  <div className="pipeline-arrow" aria-hidden="true">
                    <Icon name="arrowRight" size={18} />
                  </div>
                )}
              </Fragment>
            ))}
          </div>
        </div>

        <div className="card">
          <div className="card-head">
            <div>
              <h2>Supporting tools</h2>
              <span className="label">Useful after or alongside the core design-to-cloning path</span>
            </div>
          </div>
          <div className="supporting-tools-grid">
            {SUPPORTING_TOOLS.map((tool) => (
              <Link to={tool.to} className="pipeline-stage" key={tool.to}>
                <Icon name={tool.icon} size={24} />
                <strong>{tool.name}</strong>
                <span>{tool.blurb}</span>
              </Link>
            ))}
          </div>
        </div>

        <div className="home-grid">
          <div className="card">
            <div className="card-head">
              <h2 style={{ flex: 1 }}>Recent projects</h2>
              <Link to="/projects" className="label">all projects →</Link>
            </div>
            {error ? (
              <div
                className="card-body"
                style={{ color: "var(--muted)", fontSize: "0.88rem" }}
                role="alert"
              >
                Could not load your projects.
              </div>
            ) : projects === null ? (
              <div
                className="card-body"
                style={{ color: "var(--muted)", fontSize: "0.88rem" }}
                role="status"
                aria-busy="true"
              >
                <span className="spinner" /> Loading…
              </div>
            ) : recent.length === 0 ? (
              <div className="card-body" style={{ color: "var(--muted)", fontSize: "0.88rem" }}>
                Nothing saved yet. Start with <Link to="/design">Design</Link> or{" "}
                <Link to="/projects">import a sequence</Link>.
              </div>
            ) : (
              <div className="feature-list">
                {recent.map((p) => (
                  <Link key={p.id} to={`/projects/${p.id}`} className="feature-row" style={{ textDecoration: "none" }}>
                    <span className="dot" style={{ background: "var(--accent)" }} />
                    <span style={{ minWidth: 0 }}>
                      <span className="nm" style={{ display: "block" }}>{p.name}</span>
                      <span className="ty">{p.module.replace(/_/g, " ")}</span>
                    </span>
                    <span className="rg">{formatDate(p.updated_at)}</span>
                  </Link>
                ))}
              </div>
            )}
          </div>

          <div className="card">
            <div className="card-head"><h2>What makes this different</h2></div>
            <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.9rem" }}>
              {GUARANTEES.map((g) => (
                <div key={g.title}>
                  <strong style={{ display: "block", fontSize: "0.9rem", marginBottom: "0.15rem" }}>
                    {g.title}
                  </strong>
                  <p className="note">{g.body}</p>
                </div>
              ))}
            </div>
          </div>
        </div>

        <footer className="home-footer">
          Designed by <strong>Prof. Merzoug Mohamed</strong> &middot; Full Professor, Genomics
          Technology Platform, Higher School of Biological Sciences of Oran &middot;{" "}
          <a href="mailto:mohamed.merzoug.essbo@gmail.com">mohamed.merzoug.essbo@gmail.com</a>
        </footer>
      </div>
    </>
  );
}
