#!/usr/bin/env python3
import re
import csv
from pathlib import Path


def parse_log(fp: Path):
    text = fp.read_text(errors="ignore")
    # Sample size
    m_n = re.search(r"Loaded (\d+) sequences", text)
    n = int(m_n.group(1)) if m_n else None
    # Kappa lines
    per = []
    for line in text.splitlines():
        line = line.strip()
        m = re.match(r"kappa=(\d+\.\d+)\s+HEX\(mean\xc2\xb1sd\)=(\d+\.\d+).*", line)
        if not m:
            m = re.match(r"kappa=(\d+\.\d+)\s+HEX=(\d+\.\d+).*", line)
        if m:
            k = float(m.group(1))
            hx = float(m.group(2))
            per.append((k, hx))
    if not per:
        return n, None, []
    k_best, hx_best = sorted(per, key=lambda x: x[1])[0]
    return n, (k_best, hx_best), per


def main() -> None:
    root = Path.home() / "RNA-Validation" / "Biosphere_codec"
    logs = {
        "SARS-CoV-2": root / "sweep_sars_public.log",
        "Influenza_A": root / "sweep_infa_public.log",
        "HCV": root / "sweep_hcv_public.log",
        "DENV": root / "sweep_denv_public.log",
        "HIV1_hash": root / "sweep_hiv1_public.log",
        "HIV1_subtype": root / "sweep_hiv_subtype.log",
    }
    rows = []
    for name, fp in logs.items():
        if not fp.exists():
            continue
        n, best, per = parse_log(fp)
        if best is None:
            continue
        k_best, hx_best = best
        rows.append({
            "dataset": name,
            "n": n or "",
            "kappa_best": f"{k_best:.2f}",
            "hex_best": f"{hx_best:.2f}",
        })

    out_dir = root / "publication_outputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "sweep_summary.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dataset", "n", "kappa_best", "hex_best"])
        w.writeheader()
        w.writerows(rows)

    md_lines = ["| Dataset | N | Best κ | HEX |", "|---|---:|---:|---:|"]
    for r in rows:
        md_lines.append(f"| {r['dataset']} | {r['n']} | {r['kappa_best']} | {r['hex_best']} |")
    md_path = out_dir / "sweep_summary.md"
    md_path.write_text("\n".join(md_lines))

    print("Wrote", csv_path)
    print("Wrote", md_path)


if __name__ == "__main__":
    main()







