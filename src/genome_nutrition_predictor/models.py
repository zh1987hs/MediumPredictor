from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class AnnotationRecord:
    """Gene-level functional annotation used as evidence."""

    gene_id: str
    kos: list[str] = field(default_factory=list)
    ecs: list[str] = field(default_factory=list)
    source: str = "unknown"
    confidence: float = 0.7


@dataclass
class PathwayResult:
    """Scored pathway/module evaluation result."""

    pathway_id: str
    name: str
    category: str
    status: str
    score_0_100: float
    essential_ratio: float
    optional_ratio: float
    missing_essential_groups: list[str]
    evidence_genes: list[str]


@dataclass
class CapabilityResult:
    """Result for substrate utilization capability and preference."""

    substrate: str
    kind: str
    capability: bool
    preference_score: float
    transport_score: float
    catabolism_score: float
    regulation_score: float
    evidence_gene_ids: list[str]
    evidence_kos: list[str]


JSONDict = dict[str, Any]
