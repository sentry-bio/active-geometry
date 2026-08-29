"""Offline math for the map-query probe (no live API)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments" / "map_query"))

from panel import canonical_kmers, mash_distance  # noqa: E402
from analyze import drop_radial_head, interpret_vectors  # noqa: E402
from triangleccs.datum.form import Form  # noqa: E402


def test_identical_sequences_zero_mash():
    s = "ACGT" * 40
    d = mash_distance(canonical_kmers(s), canonical_kmers(s))
    assert d == 0.0


def test_random_vs_poly_at_far():
    a = canonical_kmers("ACGT" * 80)
    b = canonical_kmers("AT" * 160)
    assert mash_distance(a, b) > 0.1


def test_drop_radial_head_finds_correlated_axis():
    rng = np.random.default_rng(0)
    n, d = 20, 8
    ang = rng.normal(size=(n, d))
    r = np.linspace(0.2, 0.8, n)
    ang[:, 3] = r * 2.0 + 0.01 * rng.normal(size=n)
    out = drop_radial_head(ang, r)
    assert out["dropped"] is True
    assert out["radial_index"] == 3
    assert out["kept_dim"] == d - 1


def test_interpret_small_norm_as_poincare():
    form = Form()
    x = np.zeros((4, 8))
    x[:, 0] = 0.1
    out = interpret_vectors(x, form)
    assert out["kind"] == "poincare_coords"
