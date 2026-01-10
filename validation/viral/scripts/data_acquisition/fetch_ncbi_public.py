#!/usr/bin/env python3
import sys
import time
import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterator, Tuple, List


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def http_get(url: str) -> bytes:
    with urllib.request.urlopen(url) as r:
        return r.read()


def esearch_ids(db: str, query: str, retmax: int = 50000) -> List[str]:
    params = {
        "db": db,
        "term": query,
        "retmax": str(retmax),
        "retmode": "json",
    }
    url = f"{EUTILS}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    data = json.loads(http_get(url).decode("utf-8"))
    ids = data.get("esearchresult", {}).get("idlist", [])
    return ids


def efetch_fasta(db: str, ids: List[str], chunk: int = 200) -> str:
    fasta_parts: List[str] = []
    for i in range(0, len(ids), chunk):
        sub = ids[i : i + chunk]
        params = {
            "db": db,
            "id": ",".join(sub),
            "rettype": "fasta",
            "retmode": "text",
        }
        url = f"{EUTILS}/efetch.fcgi?{urllib.parse.urlencode(params)}"
        fasta_parts.append(http_get(url).decode("utf-8", errors="ignore"))
        time.sleep(0.34)  # be gentle with NCBI
    return "".join(fasta_parts)


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    seq = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None and seq:
                    yield header, "".join(seq)
                header = line.strip()
                seq = []
            else:
                s = line.strip().upper()
                if s:
                    seq.append(s)
        if header is not None and seq:
            yield header, "".join(seq)


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


def fetch_and_write(query: str, raw_out: Path) -> int:
    ids = esearch_ids("nucleotide", query)
    if not ids:
        raw_out.write_text("")
        return 0
    fasta = efetch_fasta("nucleotide", ids)
    raw_out.write_text(fasta)
    return fasta.count("\n>") + (1 if fasta.startswith(">") else 0)


def main() -> None:
    base = Path("/home/rohit/public_viruses")
    base.mkdir(parents=True, exist_ok=True)

    # Queries tuned for complete genomes; broaden if needed
    hcv_q = "Hepatitis C virus[Organism] AND (complete genome)"
    denv_q = "Dengue virus[Organism] AND (complete genome)"

    hcv_raw = base / "HCV_public_raw.fasta"
    denv_raw = base / "DENV_public_raw.fasta"

    print("Querying NCBI E-utilities ...", flush=True)
    hcv_n = fetch_and_write(hcv_q, hcv_raw)
    denv_n = fetch_and_write(denv_q, denv_raw)
    print(f"HCV raw headers: {hcv_n}")
    print(f"DENV raw headers: {denv_n}")

    hcv_out = Path("/home/rohit/HCV_multi.fasta")
    denv_out = Path("/home/rohit/DENV_multi.fasta")
    hcv_keep = filter_len(hcv_raw, hcv_out, 9000, 13000, 10000)
    denv_keep = filter_len(denv_raw, denv_out, 9000, 13000, 10000)
    print(f"HCV kept: {hcv_keep} -> {hcv_out}")
    print(f"DENV kept: {denv_keep} -> {denv_out}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)







