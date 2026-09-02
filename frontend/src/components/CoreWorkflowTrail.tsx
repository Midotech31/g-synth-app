import { Link } from "react-router-dom";

import Icon, { type IconName } from "./Icon";

type StageKey = "design" | "hybridization" | "cloning";

const STAGES: Array<{
  key: StageKey;
  label: string;
  detail: string;
  to: string;
  icon: IconName;
}> = [
  {
    key: "design",
    label: "Design",
    detail: "SSD / ESD molecules",
    to: "/design",
    icon: "helix",
  },
  {
    key: "hybridization",
    label: "Hybridization",
    detail: "Verify both strands",
    to: "/hybridize",
    icon: "check",
  },
  {
    key: "cloning",
    label: "Restriction cloning",
    detail: "Digest and ligate",
    to: "/clone",
    icon: "plate",
  },
];

export default function CoreWorkflowTrail({ active }: { active: StageKey }) {
  return (
    <nav className="core-workflow-trail" aria-label="Core G-Synth workflow">
      {STAGES.map((stage, index) => (
        <div className="core-workflow-stage" key={stage.key}>
          {index > 0 && <Icon name="arrowRight" size={15} className="core-workflow-arrow" />}
          <Link
            to={stage.to}
            className={stage.key === active ? "active" : ""}
            aria-current={stage.key === active ? "step" : undefined}
          >
            <Icon name={stage.icon} size={16} />
            <span>
              <strong>{index + 1}. {stage.label}</strong>
              <small>{stage.detail}</small>
            </span>
          </Link>
        </div>
      ))}
    </nav>
  );
}
