"""Balloon witnesses and phase0 certify."""

from __future__ import annotations

import math

from triangleccs.tape.balloon import (
    fit_dstar_logL,
    phase0_certify,
    run_balloon_cell,
)


def test_phase0_certify():
    out = phase0_certify(seed=0)
    assert out["passed"] is True
    assert out["jc_formula_max_dev"] < 1e-9


def test_balloon_jc_resolvability_drops_with_depth():
    shallow = run_balloon_cell(n_tips=16, depth=0.3, seq_len=400, seed=1, max_quartets=400)
    deep = run_balloon_cell(n_tips=16, depth=2.5, seq_len=400, seed=1, max_quartets=400)
    assert deep["jc69"].block_fraction >= 0.9
    assert deep["jc69"].resolvable_fraction < shallow["jc69"].resolvable_fraction
    assert deep["jc69"].resolvable_fraction < 0.2
    # Infinite-sites recovers topology; JC among still-resolvable need not be perfect.
    assert deep["infinite_sites"].exact_among_resolvable >= 0.95
    assert deep["infinite_sites"].resolvable_fraction == 1.0


def test_longer_sequences_delay_channel_exhaustion():
    short = run_balloon_cell(n_tips=12, depth=1.8, seq_len=80, seed=2, max_quartets=300)
    long = run_balloon_cell(n_tips=12, depth=1.8, seq_len=500, seed=2, max_quartets=300)
    assert long["jc69"].resolvable_fraction >= short["jc69"].resolvable_fraction - 1e-9


def test_dstar_logL_fitter_recovers_reported_slope():
    lengths = (200, 500, 2000, 20000)
    dstar = tuple(-0.043 + 0.244 * math.log(L) for L in lengths)
    fit = fit_dstar_logL(lengths, dstar)
    assert abs(fit["slope_lnL"] - 0.244) < 1e-9
    assert abs(fit["intercept"] + 0.043) < 1e-9
    assert fit["r_squared"] > 0.999


def test_freeze_gate_rejects_warm_start():
    from triangleccs.datum.freeze import evaluate_freeze_gate

    r = evaluate_freeze_gate(
        angular_median_deg=5.0,
        determinate_quartet_agreement=0.95,
        inheritance="warm_start",
        balloon_reported=True,
        sextant_reported=True,
    )
    assert r["passed"] is False
    assert r["theta_status"] == "candidate"
