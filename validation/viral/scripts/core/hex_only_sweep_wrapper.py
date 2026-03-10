#!/usr/bin/env python3
"""Wrapper for hex_only_sweep with extended argument parsing."""

import argparse
import sys
from sweep_fixed_kappa import run_sweep

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--data_path", type=str, required=True)
    p.add_argument("--log_dir", type=str, required=True)
    p.add_argument("--kappa_grid", type=str, default="1.25,1.30,1.35,1.40,1.45,1.50,1.55")
    p.add_argument("--seeds", type=str, default="13,17,23")
    p.add_argument("--min_length", type=int, default=1000)
    p.add_argument("--max_length", type=int, default=50000)
    p.add_argument("--max_sequences", type=int, default=10000)
    p.add_argument("--hex_weight", type=float, default=0.8)
    p.add_argument("--batch_size", type=int, default=4)
    p.add_argument("--steps_per_kappa", type=int, default=150)
    args = p.parse_args()
    
    kappas = [float(x) for x in args.kappa_grid.split(",")]
    seeds = [int(x) for x in args.seeds.split(",")]
    
    print(f"Running sweep: {args.log_dir}")
    print(f"  Data: {args.data_path}")
    print(f"  κ grid: {kappas}")
    print(f"  Seeds: {seeds}")
    print(f"  Max seq: {args.max_sequences}, Max len: {args.max_length}")
    print()
    
    run_sweep(
        data_path=args.data_path,
        name=args.log_dir,
        max_seq=args.max_sequences,
        max_len=args.max_length,
        batch=args.batch_size,
        kappas=kappas,
        seeds=seeds,
    )


