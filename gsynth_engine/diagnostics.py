"""Stable, machine-readable verdicts shared by every design workflow.

Human prose changes as the interface improves. Codes and statuses do not:
they are the contract saved projects, exports and clients can rely on.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Literal

CheckStatus = Literal["pass", "review", "block"]
Severity = Literal["info", "warning", "error"]
Verdict = Literal["ready", "review", "blocked"]


@dataclass(frozen=True)
class Diagnostic:
    code: str
    severity: Severity
    message: str
    remedy: str = ""
    positions: tuple[int, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class PreflightCheck:
    code: str
    label: str
    status: CheckStatus
    detail: str
    remedy: str = ""
    evidence: dict[str, object] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return self.status == "pass"

    def to_dict(self) -> dict[str, object]:
        return {**asdict(self), "passed": self.passed}


@dataclass(frozen=True)
class PreflightReport:
    workflow: str
    checks: tuple[PreflightCheck, ...]
    diagnostics: tuple[Diagnostic, ...] = ()

    @property
    def verdict(self) -> Verdict:
        statuses = {check.status for check in self.checks}
        if "block" in statuses:
            return "blocked"
        if "review" in statuses:
            return "review"
        return "ready"

    @property
    def can_export(self) -> bool:
        return self.verdict != "blocked"

    def to_dict(self) -> dict[str, object]:
        return {
            "workflow": self.workflow,
            "verdict": self.verdict,
            "can_export": self.can_export,
            "checks": [check.to_dict() for check in self.checks],
            "diagnostics": [diagnostic.to_dict() for diagnostic in self.diagnostics],
        }


def diagnostic_for(check: PreflightCheck) -> Diagnostic | None:
    """Convert a non-passing check into its stable diagnostic."""
    if check.status == "pass":
        return None
    return Diagnostic(
        code=check.code,
        severity="error" if check.status == "block" else "warning",
        message=check.detail,
        remedy=check.remedy,
    )


def make_report(workflow: str, checks: list[PreflightCheck]) -> PreflightReport:
    diagnostics = tuple(filter(None, (diagnostic_for(check) for check in checks)))
    return PreflightReport(workflow, tuple(checks), diagnostics)  # type: ignore[arg-type]
