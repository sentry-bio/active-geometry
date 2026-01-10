#!/usr/bin/env python3
import sys
import time
import json
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def http_get(url: str) -> bytes:
    with urllib.request.urlopen(url) as r:
        return r.read()


def esearch_ids(db: str, query: str, retmax: int = 50000):
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode({
        "db": db,
        "term": query,
        "retmax": str(retmax),
        "retmode": "json",
    })
    data = json.loads(http_get(url).decode("utf-8"))
    return data.get("esearchresult", {}).get("idlist", [])


def efetch_fasta(db: str, ids, chunk: int = 200) -> str:
    parts = []
    for i in range(0, len(ids), chunk):
        sub = ids[i : i + chunk]
        url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode({
            "db": db,
            "id": ",".join(sub),
            "rettype": "fasta",
            "retmode": "text",
        })
        parts.append(http_get(url).decode("utf-8", errors="ignore"))
        time.sleep(0.34)
    return "".join(parts)


def write_capped_fasta_text(fasta_text: str, out_path: Path, max_keep: int = 10000) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with open(out_path, "w") as out:
        h = None
        seq_parts = []
        for line in fasta_text.splitlines():
            if line.startswith(">"):
                if h is not None and seq_parts:
                    out.write(h + "\n")
                    seq = "".join(seq_parts).upper()
                    for i in range(0, len(seq), 80):
                        out.write(seq[i : i + 80] + "\n")
                    kept += 1
                    if kept >= max_keep:
                        break
                h = line.strip()
                seq_parts = []
            else:
                t = line.strip()
                if t:
                    seq_parts.append(t)
        if kept < max_keep and h is not None and seq_parts:
            out.write(h + "\n")
            seq = "".join(seq_parts).upper()
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")
            kept += 1
    return kept


def main() -> None:
    # Capture both complete genomes and segments to ensure robust sample
    query = "Influenza A virus[Organism] AND (complete genome OR segment)"
    ids = esearch_ids("nucleotide", query, retmax=50000)
    print("Influenza A: ids=", len(ids))
    fasta = efetch_fasta("nucleotide", ids)
    approx_headers = fasta.count("\n>") + (1 if fasta.startswith(">") else 0)
    print("Fetched records ~=", approx_headers)
    out = Path("/home/rohit/INFA_multi.fasta")
    kept = write_capped_fasta_text(fasta, out, max_keep=10000)
    print(f"Wrote {out} kept {kept}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ERROR:", e, file=sys.stderr)
        sys.exit(1)







