"""Immutable Form — the frozen geodetic frame.

κ is a CONVENTION (InfoNCE gauge), not (h ln 2)².
dim = 2 is inhabit H² (embeddability floor), not a fitted discovery.
Radius is advisory; θ is candidate until the freeze-gate passes.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from typing import Literal

from triangleccs.ledger import Tag


@dataclass(frozen=True)
class FreezeGate:
    """Conformance thresholds. Candidate until all pass."""

    angular_median_deg_max: float = 10.0
    determinate_quartet_agreement_min: float = 0.90


@dataclass(frozen=True)
class Form:
    version: str = "triangleccs-1.0"
    kappa: float = 1.25  # CONVENTION — frozen gauge, not a theorem
    dim: Literal[2] = 2  # CONVENTION — inhabit H²
    tokenizer_id: str = "BPE-4096"
    prime_meridian: str = "GCF_000005845.2"  # E. coli K-12 → θ = 0
    chirality_anchor: str = "GCF_000091665.1"  # M. jannaschii → handedness
    epsilon: float = 1e-3
    radial_proxy: str = "ssu+cog+kmer_entropy"
    certified_axis: Literal["theta"] = "theta"
    advisory_axis: Literal["radius"] = "radius"
    freeze_gate: FreezeGate = field(default_factory=FreezeGate)

    def __post_init__(self) -> None:
        if self.kappa <= 0:
            raise ValueError("kappa must be positive")
        if self.dim != 2:
            raise ValueError("v1 Form inhabits H^2 only (dim must be 2)")
        if self.epsilon <= 0:
            raise ValueError("epsilon must be positive")

    @property
    def ball_radius(self) -> float:
        return 1.0 / (self.kappa**0.5)

    @property
    def form_hash(self) -> str:
        payload = {
            "version": self.version,
            "kappa": self.kappa,
            "dim": self.dim,
            "tokenizer_id": self.tokenizer_id,
            "prime_meridian": self.prime_meridian,
            "chirality_anchor": self.chirality_anchor,
            "epsilon": self.epsilon,
            "radial_proxy": self.radial_proxy,
            "certified_axis": self.certified_axis,
            "advisory_axis": self.advisory_axis,
            "freeze_gate": asdict(self.freeze_gate),
        }
        blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]

    def summary(self) -> dict[str, object]:
        return {
            "version": self.version,
            "form_hash": self.form_hash,
            "kappa": self.kappa,
            "kappa_status": Tag.CONVENTION.value,
            "dim": self.dim,
            "dim_status": Tag.CONVENTION.value,
            "ball_radius": round(self.ball_radius, 6),
            "tokenizer_id": self.tokenizer_id,
            "prime_meridian": self.prime_meridian,
            "chirality_anchor": self.chirality_anchor,
            "epsilon": self.epsilon,
            "radial_proxy": self.radial_proxy,
            "certified_axis": self.certified_axis,
            "advisory_axis": self.advisory_axis,
            "freeze_gate": asdict(self.freeze_gate),
            "canonical_primitive": "distance",
            "note": (
                "Decoder datum for a tree source on an overwriting tape. "
                "κ is CONVENTION; origin is chart origin, not LUCA; "
                "r is ADVISORY."
            ),
        }


DEFAULT_FORM = Form()
