#!/usr/bin/env python3
"""Fetch all validation datasets for κ sweeps with retry logic."""

import sys
import time
import json
import hashlib
from pathlib import Path
from urllib import request, parse, error

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = "biosphere-codec/validation"

def http_get(url, params=None, retries=3):
    """HTTP GET with exponential backoff retry."""
    for attempt in range(retries):
        try:
            if params:
                url_final = url + ("?" + parse.urlencode(params))
            else:
                url_final = url
            req = request.Request(url_final, headers={"User-Agent": UA})
            with request.urlopen(req, timeout=120) as resp:
                return resp.read()
        except error.HTTPError as e:
            if attempt == retries - 1:
                raise
            wait = 2 ** attempt
            print(f"  HTTP {e.code}, retrying in {wait}s...")
            time.sleep(wait)
        except Exception as e:
            if attempt == retries - 1:
                raise
            wait = 2 ** attempt
            print(f"  Error {e}, retrying in {wait}s...")
            time.sleep(wait)

def esearch(term):
    """Search NCBI and return count + WebEnv."""
    data = http_get(BASE + "/esearch.fcgi", {
        "db": "nuccore",
        "term": term,
        "retmax": 100000,
        "usehistory": "y",
        "retmode": "json",
    })
    js = json.loads(data.decode("utf-8"))
    res = js["esearchresult"]
    return int(res["count"]), res.get("webenv"), res.get("querykey")

def efetch_batch(id_list, batch_size=200):
    """Fetch sequences by ID list in small batches using POST."""
    all_text = []
    for i in range(0, len(id_list), batch_size):
        batch_ids = id_list[i:i+batch_size]
        id_str = ",".join(batch_ids)
        
        # Use POST to avoid URL length limits
        post_data = parse.urlencode({
            "db": "nuccore",
            "id": id_str,
            "rettype": "fasta",
            "retmode": "text",
        }).encode("utf-8")
        
        for attempt in range(3):
            try:
                req = request.Request(BASE + "/efetch.fcgi", data=post_data, headers={"User-Agent": UA})
                with request.urlopen(req, timeout=120) as resp:
                    txt = resp.read().decode("utf-8", errors="ignore")
                    all_text.append(txt)
                    break
            except Exception as e:
                if attempt == 2:
                    raise
                wait = 2 ** attempt
                print(f"  Retry in {wait}s...", end="\r")
                time.sleep(wait)
        
        print(f"  fetched {i+len(batch_ids)}/{len(id_list)}", end="\r")
        time.sleep(0.34)
    return "".join(all_text)

def parse_fasta(text):
    """Parse FASTA text into (header, sequence) tuples."""
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

def fetch_and_filter(name, term, out_path, min_len, max_len):
    """Fetch sequences via ID list method (more reliable than WebEnv)."""
    print(f"\n{name}: Starting esearch...")
    count, webenv, qk = esearch(term)
    print(f"{name}: esearch count={count}")
    
    if count == 0:
        return 0, 0, 0
    
    # Get ID list
    print(f"{name}: Fetching ID list...")
    id_data = http_get(BASE + "/esearch.fcgi", {
        "db": "nuccore",
        "term": term,
        "retmax": min(count, 50000),  # Cap at 50k
        "retmode": "json",
    })
    id_list = json.loads(id_data.decode("utf-8"))["esearchresult"]["idlist"]
    print(f"{name}: Got {len(id_list)} IDs")
    
    # Fetch in batches
    print(f"{name}: Fetching sequences...")
    fasta_text = efetch_batch(id_list, batch_size=500)
    
    # Filter and deduplicate
    print(f"\n{name}: Filtering by length [{min_len}, {max_len}]...")
    kept = 0
    fetched = 0
    seen = set()
    
    with open(out_path, "w") as out_f:
        for h, s in parse_fasta(fasta_text):
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
    
    print(f"{name}: DONE - fetched={fetched}, kept={kept}")
    return count, fetched, kept

def main():
    root = Path.cwd()
    root.mkdir(parents=True, exist_ok=True)
    
    # All validation datasets
    jobs = [
        # Flaviviridae (continued)
        ("Zika", '"Zika virus"[Organism] AND biomol_genomic[PROP] AND 9000:13000[SLEN]', 
         "zika_multi.fasta", 9000, 13000),
        ("WestNile", '"West Nile virus"[Organism] AND biomol_genomic[PROP] AND 9000:13000[SLEN]', 
         "west_nile_multi.fasta", 9000, 13000),
        ("YellowFever", '"Yellow fever virus"[Organism] AND biomol_genomic[PROP] AND 9000:13000[SLEN]', 
         "yellow_fever_multi.fasta", 9000, 13000),
        
        # Picornaviridae
        ("Poliovirus", '"Poliovirus"[Organism] AND biomol_genomic[PROP] AND 7000:8000[SLEN]', 
         "poliovirus_multi.fasta", 7000, 8000),
        ("Enterovirus", '"Enterovirus"[Organism] AND biomol_genomic[PROP] AND 7000:8000[SLEN]', 
         "enterovirus_multi.fasta", 7000, 8000),
        
        # Paramyxoviridae
        ("Measles", '"Measles morbillivirus"[Organism] AND biomol_genomic[PROP] AND 15000:17000[SLEN]', 
         "measles_multi.fasta", 15000, 17000),
        ("Mumps", '"Mumps orthorubulavirus"[Organism] AND biomol_genomic[PROP] AND 15000:17000[SLEN]', 
         "mumps_multi.fasta", 15000, 17000),
        
        # Rhabdoviridae
        ("Rabies", '"Rabies lyssavirus"[Organism] AND biomol_genomic[PROP] AND 11000:13000[SLEN]', 
         "rabies_multi.fasta", 11000, 13000),
        
        # Filoviridae
        ("Ebola", '"Ebolavirus"[Organism] AND biomol_genomic[PROP] AND 18000:20000[SLEN]', 
         "ebola_multi.fasta", 18000, 20000),
        ("Marburg", '"Marburgvirus"[Organism] AND biomol_genomic[PROP] AND 18000:20000[SLEN]', 
         "marburg_multi.fasta", 18000, 20000),
        
        # DNA virus control (Herpesviridae - expect κ ≈ 1.25)
        ("HSV1", '"Human alphaherpesvirus 1"[Organism] AND biomol_genomic[PROP] AND 100000:200000[SLEN]', 
         "hsv1_multi.fasta", 100000, 200000),
        ("CMV", '"Human betaherpesvirus 5"[Organism] AND biomol_genomic[PROP] AND 200000:250000[SLEN]', 
         "cmv_multi.fasta", 200000, 250000),
    ]
    
    summary = {}
    for name, term, fname, mn, mx in jobs:
        out_path = root / fname
        if out_path.exists():
            print(f"\n{name}: Already exists, skipping")
            # Count existing sequences
            with open(out_path) as f:
                existing_count = sum(1 for line in f if line.startswith(">"))
            summary[name] = {"status": "existing", "kept": existing_count, "path": str(out_path)}
            continue
        
        try:
            cnt, fet, kep = fetch_and_filter(name, term, str(out_path), mn, mx)
            summary[name] = {
                "esearch_count": cnt, 
                "fetched": fet, 
                "kept": kep, 
                "path": str(out_path)
            }
        except Exception as e:
            print(f"\nERROR {name}: {e}")
            import traceback
            traceback.print_exc()
            summary[name] = {"error": str(e)}
    
    # Save summary
    summary_file = root / "fetch_summary_all_validation.json"
    summary_file.write_text(json.dumps(summary, indent=2))
    
    print("\n" + "="*60)
    print("FINAL SUMMARY:")
    for name, info in summary.items():
        if "error" in info:
            print(f"  {name}: ERROR - {info['error']}")
        elif "status" in info and info["status"] == "existing":
            print(f"  {name}: ✓ Existing ({info['kept']} seqs)")
        else:
            print(f"  {name}: ✓ Fetched {info['kept']} seqs (from {info['esearch_count']} available)")
    print("="*60)
    print(f"Summary saved: {summary_file}")

if __name__ == "__main__":
    main()
