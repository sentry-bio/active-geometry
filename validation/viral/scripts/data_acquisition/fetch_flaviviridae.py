#!/usr/bin/env python3
"""Fetch Flaviviridae sequences from NCBI for κ validation sweeps."""

import sys
import time
import json
import hashlib
from pathlib import Path
from urllib import request, parse

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = "biosphere-codec/validation"

def http_get(url, params=None):
    if params:
        url = url + ("?" + parse.urlencode(params))
    req = request.Request(url, headers={"User-Agent": UA})
    with request.urlopen(req, timeout=120) as resp:
        return resp.read()

def esearch(term):
    data = http_get(BASE + "/esearch.fcgi", {
        "db": "nuccore",
        "term": term,
        "retmax": 100000,
        "usehistory": "y",
        "retmode": "json",
    })
    js = json.loads(data.decode("utf-8"))
    res = js["esearchresult"]
    return int(res["count"]), res["webenv"], res["querykey"]

def efetch_chunks(webenv, query_key, chunk=10000):
    retstart = 0
    while True:
        txt = http_get(BASE + "/efetch.fcgi", {
            "db": "nuccore",
            "query_key": query_key,
            "WebEnv": webenv,
            "rettype": "fasta",
            "retmode": "text",
            "retstart": retstart,
            "retmax": chunk,
        }).decode("utf-8", errors="ignore")
        if not txt.strip():
            break
        yield txt
        retstart += chunk
        time.sleep(0.34)

def parse_fasta(text):
    header = None
    seq_lines = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_lines)
            header = line.strip()
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        yield header, "".join(seq_lines)

def stream_filter_and_write(name, term, out_path, min_len, max_len):
    print(f"\n{name}: Starting fetch...")
    count, webenv, qk = esearch(term)
    print(f"{name}: esearch count={count}")
    
    kept = 0
    fetched = 0
    seen = set()
    
    with open(out_path, "w") as out_f:
        for chunk in efetch_chunks(webenv, qk):
            for h, s in parse_fasta(chunk):
                fetched += 1
                L = len(s)
                if L < min_len or L > max_len:
                    continue
                key = hashlib.sha256(s.encode()).hexdigest()
                if key in seen:
                    continue
                seen.add(key)
                kept += 1
                out_f.write(f"{h}\n")
                for i in range(0, L, 80):
                    out_f.write(s[i:i+80] + "\n")
            print(f"{name}: fetched={fetched}, kept={kept}", end="\r")
    
    print(f"\n{name}: DONE - fetched={fetched}, kept={kept}")
    return count, fetched, kept

def main():
    root = Path.cwd()
    root.mkdir(parents=True, exist_ok=True)
    
    jobs = [
        ("Zika", '"Zika virus"[Organism] AND biomol_genomic[PROP] AND 9000:13000[SLEN]', 
         root/"zika_multi.fasta", 9000, 13000),
        ("WestNile", '"West Nile virus"[Organism] AND biomol_genomic[PROP] AND 9000:13000[SLEN]', 
         root/"west_nile_multi.fasta", 9000, 13000),
        ("YellowFever", '"Yellow fever virus"[Organism] AND biomol_genomic[PROP] AND 9000:13000[SLEN]', 
         root/"yellow_fever_multi.fasta", 9000, 13000),
    ]
    
    summary = {}
    for name, term, path, mn, mx in jobs:
        try:
            cnt, fet, kep = stream_filter_and_write(name, term, str(path), mn, mx)
            summary[name] = {
                "esearch_count": cnt, 
                "fetched": fet, 
                "kept": kep, 
                "path": str(path)
            }
        except Exception as e:
            print(f"ERROR {name}: {e}")
            summary[name] = {"error": str(e)}
    
    summary_file = Path("fetch_summary_flaviviridae.json")
    summary_file.write_text(json.dumps(summary, indent=2))
    
    print("\n" + "="*60)
    print("SUMMARY:")
    print(json.dumps(summary, indent=2))
    print("="*60)

if __name__ == "__main__":
    main()


