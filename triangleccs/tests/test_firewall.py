"""Constitution firewall: overlays cannot write Form.kappa; no (h,κ)→n API."""

from __future__ import annotations

import inspect

import triangleccs
from triangleccs.datum.form import Form
import triangleccs.overlays as overlays


def test_form_is_frozen():
    f = Form()
    try:
        f.kappa = 2.0  # type: ignore[misc]
        raised = False
    except Exception:
        raised = True
    assert raised


def test_no_public_h_kappa_to_n():
    """Public package must not expose a back-solve n from (h, kappa)."""
    public_names = set(triangleccs.__all__)
    forbidden = {"backsolve_n", "n_from_kappa", "infer_n", "state_equation_n"}
    assert public_names.isdisjoint(forbidden)
    # Scan exported callables' parameter names for a suspicious pair API
    for name in public_names:
        obj = getattr(triangleccs, name)
        if not callable(obj):
            continue
        try:
            params = set(inspect.signature(obj).parameters)
        except (TypeError, ValueError):
            continue
        assert not ({"h", "kappa", "n"} <= params), name


def test_overlays_package_is_empty_of_writers():
    assert not hasattr(overlays, "set_kappa")
    assert not hasattr(overlays, "apply_state_equation")
