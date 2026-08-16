"""Sextant v0: distances onto the chart, Address instruments filled."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from triangleccs.datum.form import Form
from triangleccs.datum.registration import addresses_from_registration, load_transform
from triangleccs.sextant import place_sequences
from triangleccs.tape.balloon import simulate_jc69_on_tree, yule_tree


ROOT = Path(__file__).resolve().parents[1]


def test_sextant_fills_address_instruments():
    form = Form()
    rng = np.random.default_rng(3)
    tree = yule_tree(8, 0.4, rng)
    tips = simulate_jc69_on_tree(tree, 120, rng)
    report = place_sequences(tips, form=form, meridian_index=0, chirality_index=1)
    meridian_index, chirality_index = 0, 1
    assert len(report.addresses) == 8
    a0 = report.addresses[0]
    assert a0.delta is not None
    assert a0.resolvable is not None
    assert a0.block_sep is not None
    assert a0.residual is not None
    assert a0.tags["delta"] == "INSTRUMENT"
    assert a0.tags["r"] == "ADVISORY"
    assert abs(report.addresses[meridian_index].theta) < 1e-6
    assert report.addresses[chirality_index].theta >= -1e-12
    assert report.residual >= 0.0
    assert 0.0 <= report.block_sep <= 1.0


def test_sextant_four_strings():
    form = Form()
    seqs = [
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "ACGTACGTACGTACGTACGTTCGTACGTACGTACGT",
        "TTTTACGTACGTACGTACGTACGTACGTACGTACGT",
        "TTTTACGTACGTACGTACGTTCGTACGTACCCACGT",
    ]
    report = place_sequences(seqs, form=form)
    assert len(report.addresses) == 4
    assert all(a.form_hash == form.form_hash for a in report.addresses)


def test_registration_fills_instruments_and_uses_fixture_transform():
    form = Form()
    path = ROOT / "examples" / "atlas_transform.v1.json"
    transform = load_transform(path)
    assert transform["form_hash"] == form.form_hash
    assert transform["certified"] is False
    assert transform["packing_metric"] == "poincare"
    rng = np.random.default_rng(0)
    D = 8
    raw = rng.normal(size=(12, D))
    raw = raw / np.linalg.norm(raw, axis=1, keepdims=True) * 0.25
    coords = np.tanh(raw) * (0.4 / np.sqrt(form.kappa))
    addrs = addresses_from_registration(
        coords, transform, coords[0], coords[1], form=form
    )
    assert len(addrs) == 12
    assert addrs[0].delta is not None
    assert addrs[0].resolvable is not None
    assert addrs[0].block_sep is not None
    assert addrs[0].residual is not None
    assert addrs[0].r >= 0.0
    # Meridian θ = 0
    assert abs(addrs[0].theta) < 1e-8
