import type { ReactNode } from "react";

import Icon, { type IconName } from "./Icon";
import { Logo } from "./Logo";

type TrustPoint = {
  icon: IconName;
  title: string;
  body: string;
};

const TRUST_POINTS: TrustPoint[] = [
  {
    icon: "helix",
    title: "Synthesis-ready sequence design",
    body: "SSD and ESD turn a target sequence into explicit oligos and assembly junctions.",
  },
  {
    icon: "plate",
    title: "Build logic you can inspect",
    body: "PCR, digestion, ligation and cloning are simulated from the actual sequences.",
  },
  {
    icon: "microscope",
    title: "Trace-aware validation",
    body: "AB1 and SCF evidence is aligned back to the saved reference before a clone is accepted.",
  },
];

type Props = {
  eyebrow: string;
  title: string;
  description: string;
  children: ReactNode;
  footer: ReactNode;
};

/** Shared, evidence-led entry point for sign-in and registration. */
export default function AuthShell({ eyebrow, title, description, children, footer }: Props) {
  return (
    <div className="auth">
      <section className="auth-story" aria-label="About G-Synth">
        <Logo size={68} />
        <div className="auth-story-copy">
          <span className="auth-eyebrow">End-to-end synthetic biology workspace</span>
          <h1>Move from sequence intent to a construct you can defend.</h1>
          <p>
            G-Synth keeps design, cloning and post-sequencing evidence in one auditable
            scientific record.
          </p>
        </div>

        <ul className="auth-trust-list">
          {TRUST_POINTS.map((point) => (
            <li key={point.title}>
              <span className="auth-trust-icon" aria-hidden="true">
                <Icon name={point.icon} size={20} />
              </span>
              <span>
                <strong>{point.title}</strong>
                <small>{point.body}</small>
              </span>
            </li>
          ))}
        </ul>

        <p className="auth-privacy">
          Your sequence records remain private to your authenticated workspace.
        </p>
      </section>

      <main className="auth-panel" aria-labelledby="auth-title">
        <div className="auth-form-card">
          <div className="auth-form-head">
            <span className="auth-eyebrow">{eyebrow}</span>
            <h1 id="auth-title">{title}</h1>
            <p>{description}</p>
          </div>
          {children}
        </div>
        <p className="auth-alt">{footer}</p>
        <p className="auth-signature">
          Designed by <strong>Prof. Merzoug Mohamed</strong>
        </p>
      </main>
    </div>
  );
}
