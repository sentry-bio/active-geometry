"""Freeze-gate conformance suite.

Candidate until angular agreement and determinate-quartet agreement pass.
Warm-start inheritance alone cannot certify.
"""

from __future__ import annotations

from typing import Any, Mapping

from triangleccs.datum.form import Form


def evaluate_freeze_gate(
    *,
    angular_median_deg: float,
    determinate_quartet_agreement: float,
    inheritance: str = "independent",
    form: Form | None = None,
    balloon_reported: bool = False,
    sextant_reported: bool = False,
) -> dict[str, Any]:
    form = form or Form()
    gate = form.freeze_gate
    angular_ok = angular_median_deg <= gate.angular_median_deg_max
    quartet_ok = (
        determinate_quartet_agreement >= gate.determinate_quartet_agreement_min
    )
    inheritance_ok = inheritance != "warm_start"
    passed = bool(angular_ok and quartet_ok and inheritance_ok)
    return {
        "passed": passed,
        "theta_status": "certified" if passed else "candidate",
        "checks": {
            "angular_median_deg": {
                "value": angular_median_deg,
                "max": gate.angular_median_deg_max,
                "ok": angular_ok,
            },
            "determinate_quartet_agreement": {
                "value": determinate_quartet_agreement,
                "min": gate.determinate_quartet_agreement_min,
                "ok": quartet_ok,
            },
            "inheritance_independent": {
                "value": inheritance,
                "ok": inheritance_ok,
                "note": "warm_start is CIRCULAR and cannot alone satisfy the freeze-gate",
            },
            "balloon_witnesses_reported": {
                "value": balloon_reported,
                "ok": balloon_reported,
                "note": "reported for honesty; not a topology supervisor",
            },
            "independent_sextant_reported": {
                "value": sextant_reported,
                "ok": sextant_reported,
                "note": (
                    "load-bearing witness when available: distances onto the "
                    "chart, not a genomic LM; not required to flip passed"
                ),
            },
        },
    }


def summarize_gate(result: Mapping[str, Any]) -> str:
    status = result["theta_status"]
    return f"freeze-gate {status}: passed={result['passed']}"
