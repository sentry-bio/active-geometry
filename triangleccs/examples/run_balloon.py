#!/usr/bin/env python3
"""Run a small Yule/JC69 balloon cell and the D* = a + b ln L fitter on a fixture."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from triangleccs.tape.balloon import fit_dstar_logL, phase0_certify, run_balloon_cell


def main() -> None:
    print(json.dumps(phase0_certify(seed=0), indent=2, default=float))
    shallow = run_balloon_cell(n_tips=16, depth=0.3, seq_len=400, seed=1, max_quartets=400)
    deep = run_balloon_cell(n_tips=16, depth=2.5, seq_len=400, seed=1, max_quartets=400)
    for name, cell in ("shallow", shallow), ("deep", deep):
        jc, inf = cell["jc69"], cell["infinite_sites"]
        print(
            f"{name}: JC block={jc.block_fraction:.3f} "
            f"resolvable={jc.resolvable_fraction:.3f} "
            f"topo={jc.exact_among_resolvable:.3f} | "
            f"inf-sites resolvable={inf.resolvable_fraction:.3f} "
            f"topo={inf.exact_among_resolvable:.3f}"
        )
    # Published-scale coefficients as a fitter sanity check, not a new measurement.
    lengths = (200, 500, 2000, 20000)
    # D* = -0.043 + 0.244 ln L  (reported Yule/JC69 sweep, THROUGHLINE)
    dstar = tuple(-0.043 + 0.244 * __import__("math").log(L) for L in lengths)
    print(json.dumps(fit_dstar_logL(lengths, dstar), indent=2))


if __name__ == "__main__":
    main()
