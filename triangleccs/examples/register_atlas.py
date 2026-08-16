#!/usr/bin/env python3
"""Consumer: register fixture Atlas-like coords onto the Form; fill Address."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from triangleccs.datum.form import Form
from triangleccs.datum.registration import addresses_from_registration, load_transform
from triangleccs.sextant import place_sequences


def main() -> None:
    form = Form()
    transform = load_transform(ROOT / "examples" / "atlas_transform.v1.json")
    assert transform["form_hash"] == form.form_hash
    assert transform["certified"] is False
    rng = np.random.default_rng(0)
    D = 8
    raw = rng.normal(size=(20, D))
    raw = raw / np.linalg.norm(raw, axis=1, keepdims=True) * 0.3
    coords = np.tanh(raw) * (0.5 / np.sqrt(form.kappa))

    addrs = addresses_from_registration(
        coords, transform, coords[0], coords[1], form=form
    )
    print(json.dumps(form.summary(), indent=2))
    a0 = addrs[0]
    print(f"registered {len(addrs)} addresses; theta_status={a0.theta_status}")
    print(
        f"example: theta={a0.theta:.4f} rad, r={a0.r:.4f} "
        f"(proxy={a0.r_proxy}, ADVISORY) delta={a0.delta:.4f} "
        f"resolvable={a0.resolvable:.3f} block_sep={a0.block_sep:.3f} "
        f"residual={a0.residual:.4f}"
    )

    # Independent sextant lineage on a four-taxon JC quartet (not an LM).
    seqs = [
        "ACGTACGTACGTACGTACGTACGTACGTACGT",
        "ACGTACGTACGTACGTACGTTCGTACGTACGT",
        "GGGGACGTACGTACGTACGTACGTACGTACGT",
        "GGGGACGTACGTACGTACGTTCGTACGTACCC",
    ]
    report = place_sequences(seqs, form=form)
    print(
        f"sextant: n={len(report.addresses)} delta={report.delta:.4f} "
        f"resolvable={report.resolvable:.3f} block_sep={report.block_sep:.3f} "
        f"residual={report.residual:.4f}"
    )
    print(report.note)


if __name__ == "__main__":
    main()
