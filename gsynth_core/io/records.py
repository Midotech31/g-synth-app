"""Framework-agnostic record types (Biopython-compatible, no dep)."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Literal


@dataclass(frozen=True, slots=True)
class FeatureLocation:
    """0-based half-open location; optionally compound for split features."""
    start: int
    end: int
    strand: Literal[1, -1, 0] = 1
    compound_parts: tuple[tuple[int, int], ...] = ()

    @property
    def length(self) -> int:
        if self.compound_parts:
            return sum(e - s for s, e in self.compound_parts)
        return self.end - self.start


@dataclass(slots=True)
class Feature:
    type: str
    location: FeatureLocation
    qualifiers: dict[str, Any] = field(default_factory=dict)

    @property
    def name(self) -> str:
        return (self.qualifiers.get("label")
                or self.qualifiers.get("gene")
                or self.qualifiers.get("product")
                or self.qualifiers.get("note")
                or self.type)

    @property
    def color(self) -> str:
        return self.qualifiers.get("color", self.qualifiers.get("ApEinfo_fwdcolor", ""))


@dataclass(slots=True)
class SeqRecord:
    sequence: str
    name: str = ""
    description: str = ""
    features: list[Feature] = field(default_factory=list)
    topology: Literal["linear", "circular"] = "linear"
    molecule_type: str = "DNA"
    annotations: dict[str, Any] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.sequence)

    def add_feature(self, feature: Feature) -> None:
        self.features.append(feature)

    def features_in(self, start: int, end: int) -> list[Feature]:
        return [f for f in self.features
                if f.location.start < end and f.location.end > start]
