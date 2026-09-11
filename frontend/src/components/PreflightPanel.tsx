import type { PreflightReport } from "../api/client";
import Icon from "./Icon";

const VERDICT = {
  ready: { label: "Ready", icon: "check" as const },
  review: { label: "Review required", icon: "target" as const },
  blocked: { label: "Blocked", icon: "cross" as const },
};


export default function PreflightPanel({ report }: { report?: PreflightReport }) {
  if (!report) return null;
  const verdict = VERDICT[report.verdict];

  return (
    <section className={`card preflight preflight-${report.verdict}`} aria-labelledby="preflight-title">
      <div className="card-head">
        <h2 id="preflight-title" style={{ flex: 1 }}>Preflight</h2>
        <span className={`preflight-verdict ${report.verdict}`}>
          <Icon name={verdict.icon} size={16} />
          {verdict.label}
        </span>
      </div>
      <div className="card-body preflight-list">
        {report.checks.map((check) => (
          <details className={`preflight-check ${check.status}`} key={check.code}>
            <summary>
              <Icon name={check.status === "pass" ? "check" : check.status === "block" ? "cross" : "target"} size={16} />
              <span className="preflight-check-copy">
                <strong>{check.label}</strong>
                <span>{check.detail}</span>
              </span>
              <code>{check.code}</code>
            </summary>
            <div className="preflight-detail">
              {check.remedy && <p><strong>Next action:</strong> {check.remedy}</p>}
              {Object.keys(check.evidence).length > 0 && (
                <dl>
                  {Object.entries(check.evidence).map(([key, value]) => (
                    <div key={key}>
                      <dt>{key.replace(/_/g, " ")}</dt>
                      <dd>{typeof value === "string" ? value : JSON.stringify(value)}</dd>
                    </div>
                  ))}
                </dl>
              )}
            </div>
          </details>
        ))}
        <p className="preflight-foot">
          {report.can_export
            ? "All blocking checks passed. Review items remain visible in exported records."
            : "Export and project release stay disabled until every blocking check passes."}
        </p>
      </div>
    </section>
  );
}
