#!/usr/bin/env python3
import os
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterator, Tuple


ENA_API = "https://www.ebi.ac.uk/ena/browser/api/fasta"


def ena_url(query: str) -> str:
    params = {
        "download": "true",
        "result": "sequence",
        "query": query,
    }
    return ENA_API + "?" + urllib.parse.urlencode(params)


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as r, open(dest, "wb") as f:
        f.write(r.read())


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    parts = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None and parts:
                    yield header, "".join(parts)
                header = line.strip()
                parts = []
            else:
                s = line.strip().upper()
                if s:
                    parts.append(s)
        if header is not None and parts:
            yield header, "".join(parts)


def filter_len(inp: Path, outp: Path, min_len: int, max_len: int, max_keep: int) -> int:
    kept = 0
    with open(outp, "w") as out:
        for h, s in iter_fasta(inp):
            L = len(s)
            if L < min_len or L > max_len:
                continue
            out.write(h + "\n")
            for i in range(0, L, 80):
                out.write(s[i : i + 80] + "\n")
            kept += 1
            if kept >= max_keep:
                break
    return kept


def main() -> None:
    base = Path("./data/viruses")
    base.mkdir(parents=True, exist_ok=True)
    hcv_raw = base / "HCV_public_raw.fasta"
    denv_raw = base / "DENV_public_raw.fasta"

    # Queries target complete genomes; ENA free-text query
    hcv_q = '"Hepatitis C virus" AND "complete genome"'
    denv_q = '"Dengue virus" AND "complete genome"'

    print("Downloading from ENA...")
    download(ena_url(hcv_q), hcv_raw)
    download(ena_url(denv_q), denv_raw)
    print("Downloads complete")

    # Filter to plausible genome lengths (rough windows)
    hcv_out = Path("./data/viruses/HCV_multi.fasta")
    denv_out = Path("./data/viruses/DENV_multi.fasta")
    hcv_n = filter_len(hcv_raw, hcv_out, min_len=9000, max_len=13000, max_keep=10000)
    denv_n = filter_len(denv_raw, denv_out, min_len=9000, max_len=13000, max_keep=10000)
    print(f"HCV written: {hcv_out} count={hcv_n}")
    print(f"DENV written: {denv_out} count={denv_n}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)







