# Data Directory

This directory contains **manifests only**, not raw genomic data.
Raw data must be downloaded separately due to size (~200GB).

## Manifests

| File | Records | Description |
|------|---------|-------------|
| `genomic_5627.tsv` | 5,627 | Main training set (RefSeq) |
| `viral_15_datasets.tsv` | 89,247 | 15 viral validation sets |
| `public_refseq.tsv` | 1,540 | Public-only subset (Tier 1) |

## Manifest Format

```tsv
accession    organism    taxonomy_id    clade    source
GCF_000001    E. coli    562    Bacteria    RefSeq
...
```

## Downloading Data

### Option 1: Public Subset (1,540 genomes)

```bash
make fetch-data
```

This downloads publicly available RefSeq genomes (~50GB).

### Option 2: Full Dataset (5,627 genomes)

Contact authors for access to the complete training set.

### Option 3: Viral Datasets

```bash
cd validation/viral/scripts/data_acquisition
python fetch_all_validation_datasets.py
```

## Data Tiers

| Tier | Genomes | Access | Purpose |
|------|---------|--------|---------|
| 1 | 1,540 | Public | Independent replication |
| 2 | 5,627 | Request | Full reproduction |
| 3 | 89,247 | Public | Viral validation |

## Checksums

All downloaded files should be validated against provided SHA256 checksums
in `manifests/checksums.txt`.
