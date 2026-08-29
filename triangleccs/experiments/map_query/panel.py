"""Marker-gene panel: 16S/18S accessions spanning three domains.

This is a prototype of the shared reference set ENCODER.md requires:
sequences the map can embed and (once aligned) the sextant can place.
It is not that alignment, and it is not a freeze-gate witness panel.
"""

from __future__ import annotations

import random
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
FASTA = HERE / "sequences.fasta"

# Roles match Form anchors where possible. Prime meridian is E. coli K-12
# (GCF_000005845.2); chirality is M. jannaschii (GCF_000091665.1).
PANEL: list[dict[str, str]] = [
    {
        "id": "ecoli_k12",
        "accession": "J01859.1",
        "domain": "Bacteria",
        "clade": "Enterobacteriaceae",
        "organism": "Escherichia coli",
        "role": "prime_meridian",
        "marker": "16S",
    },
    {
        "id": "ecoli_nbrc",
        "accession": "NR_114042.1",
        "domain": "Bacteria",
        "clade": "Enterobacteriaceae",
        "organism": "Escherichia coli NBRC 102203",
        "role": "same_species",
        "marker": "16S",
    },
    {
        "id": "citrobacter",
        "accession": "NR_118106.1",
        "domain": "Bacteria",
        "clade": "Enterobacteriaceae",
        "organism": "Citrobacter amalonaticus",
        "role": "sibling",
        "marker": "16S",
    },
    {
        "id": "shigella",
        "accession": "NR_026331.1",
        "domain": "Bacteria",
        "clade": "Enterobacteriaceae",
        "organism": "Shigella flexneri",
        "role": "sibling",
        "marker": "16S",
    },
    {
        "id": "klebsiella",
        "accession": "NR_117686.1",
        "domain": "Bacteria",
        "clade": "Enterobacteriaceae",
        "organism": "Klebsiella pneumoniae",
        "role": "sibling",
        "marker": "16S",
    },
    {
        "id": "salmonella",
        "accession": "NR_074910.1",
        "domain": "Bacteria",
        "clade": "Enterobacteriaceae",
        "organism": "Salmonella enterica Typhimurium",
        "role": "sibling",
        "marker": "16S",
    },
    {
        "id": "mtb",
        "accession": "NR_044826.2",
        "domain": "Bacteria",
        "clade": "Mycobacteriaceae",
        "organism": "Mycobacterium tuberculosis H37Rv",
        "role": "distant_bacteria",
        "marker": "16S",
    },
    {
        "id": "saureus",
        "accession": "NR_118997.2",
        "domain": "Bacteria",
        "clade": "Staphylococcaceae",
        "organism": "Staphylococcus aureus",
        "role": "distant_bacteria",
        "marker": "16S",
    },
    {
        "id": "bsubtilis",
        "accession": "NR_102783.2",
        "domain": "Bacteria",
        "clade": "Bacillaceae",
        "organism": "Bacillus subtilis 168",
        "role": "distant_bacteria",
        "marker": "16S",
    },
    {
        "id": "paeruginosa",
        "accession": "OP830365.1",
        "domain": "Bacteria",
        "clade": "Pseudomonadaceae",
        "organism": "Pseudomonas aeruginosa",
        "role": "distant_bacteria",
        "marker": "16S",
    },
    {
        "id": "streptomyces",
        "accession": "NR_119342.1",
        "domain": "Bacteria",
        "clade": "Streptomycetaceae",
        "organism": "Streptomyces coelicolor",
        "role": "distant_bacteria",
        "marker": "16S",
    },
    {
        "id": "tkodak",
        "accession": "NR_028216.1",
        "domain": "Archaea",
        "clade": "Thermococcaceae",
        "organism": "Thermococcus kodakarensis KOD1",
        "role": "archaea",
        "marker": "16S",
    },
    {
        "id": "mjannaschii",
        "accession": "NR_074233.1",
        "domain": "Archaea",
        "clade": "Methanocaldococcaceae",
        "organism": "Methanocaldococcus jannaschii DSM 2661",
        "role": "chirality_anchor",
        "marker": "16S",
    },
    {
        "id": "sulfolobus",
        "accession": "NR_119198.1",
        "domain": "Archaea",
        "clade": "Sulfolobaceae",
        "organism": "Saccharolobus solfataricus",
        "role": "archaea",
        "marker": "16S",
    },
    {
        "id": "archaeoglobus",
        "accession": "NR_074334.1",
        "domain": "Archaea",
        "clade": "Archaeoglobaceae",
        "organism": "Archaeoglobus fulgidus",
        "role": "archaea",
        "marker": "16S",
    },
    {
        "id": "haloferax",
        "accession": "NR_074218.1",
        "domain": "Archaea",
        "clade": "Haloferacaceae",
        "organism": "Haloferax volcanii DS2",
        "role": "archaea",
        "marker": "16S",
    },
    {
        "id": "yeast",
        "accession": "NR_132213.1",
        "domain": "Eukaryota",
        "clade": "Saccharomycetaceae",
        "organism": "Saccharomyces cerevisiae S288C",
        "role": "eukaryote",
        "marker": "18S",
    },
    {
        "id": "arabidopsis",
        "accession": "NR_141642.1",
        "domain": "Eukaryota",
        "clade": "Brassicaceae",
        "organism": "Arabidopsis thaliana",
        "role": "eukaryote",
        "marker": "18S",
    },
]


def fetch_fasta(path: Path = FASTA) -> None:
    """Fetch the accession panel from NCBI when the local cache is absent."""
    accessions = [item["accession"] for item in PANEL]
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": ",".join(accessions),
            "rettype": "fasta",
            "retmode": "text",
        }
    )
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "triangleccs-map-query/1.0 (research@sentry.bio)"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        data = response.read()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def parse_fasta(path: Path = FASTA) -> dict[str, str]:
    if not path.exists():
        fetch_fasta(path)
    out: dict[str, str] = {}
    acc = None
    chunks: list[str] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if acc is not None:
                    out[acc] = "".join(chunks)
                acc = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip().upper().replace("U", "T"))
    if acc is not None:
        out[acc] = "".join(chunks)
    return out


def load_panel(path: Path = FASTA) -> list[dict[str, str]]:
    seqs = parse_fasta(path)
    rows = []
    for item in PANEL:
        seq = seqs.get(item["accession"])
        if not seq:
            raise KeyError(f"missing FASTA for {item['accession']}")
        rows.append({**item, "sequence": seq, "length": str(len(seq))})
    return rows


def clean_dna(seq: str) -> str:
    return "".join(ch for ch in seq.upper().replace("U", "T") if ch in "ACGT")


def canonical_kmers(seq: str, k: int = 21) -> set[str]:
    s = clean_dna(seq)
    table = str.maketrans("ACGT", "TGCA")
    rc = s.translate(table)[::-1]
    out: set[str] = set()
    for i in range(len(s) - k + 1):
        kmer = s[i : i + k]
        rkmer = rc[len(s) - i - k : len(s) - i]
        out.add(min(kmer, rkmer))
    return out


def mash_distance(a: set[str], b: set[str], k: int = 21) -> float:
    """Mash-style distance from Jaccard of canonical k-mers. INSTRUMENT."""
    if not a or not b:
        return 1.0
    inter = len(a & b)
    union = len(a | b)
    j = inter / union if union else 0.0
    if j <= 0.0:
        return 1.0
    import math

    return float(max(0.0, -math.log(2.0 * j / (1.0 + j)) / k))


def length_ladder(seq: str, lengths: tuple[int, ...] = (100, 200, 400, 800)) -> list[dict[str, str]]:
    s = clean_dna(seq)
    out = [{"kind": "prefix", "length": str(n), "sequence": s[:n]} for n in lengths if n < len(s)]
    out.append({"kind": "full", "length": str(len(s)), "sequence": s})
    if len(s) >= 400:
        out.append({"kind": "suffix400", "length": "400", "sequence": s[-400:]})
    return out


def null_sequences(ecoli: str, seed: int = 0) -> list[dict[str, str]]:
    rng = random.Random(seed)
    s = clean_dna(ecoli)
    shuffled = list(s)
    rng.shuffle(shuffled)
    alphabet = "ACGT"
    rnd = "".join(rng.choice(alphabet) for _ in range(len(s)))
    short = "".join(rng.choice(alphabet) for _ in range(80))
    return [
        {"id": "shuffled_ecoli", "kind": "shuffled", "sequence": "".join(shuffled)},
        {"id": "random_dna", "kind": "random", "sequence": rnd},
        {"id": "poly_at", "kind": "low_complexity", "sequence": ("AT" * (len(s) // 2 + 1))[: len(s)]},
        {"id": "short_random", "kind": "short", "sequence": short},
    ]
