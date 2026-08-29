# Biosphere Address Registry v0.1.0 snapshot

Commit: `038529e55c439b1dfa16fbda3e95dd251a7a05aa`  
Branch: `cursor/establish-biosphere-registry-887b`  
Tests: 13 passed locally (`python3 -m pytest -q`)

This snapshot is the BAR standard and reference code. It is **not** a certified
map, freeze-gate, or species globe.

## Files

| File | Use |
|---|---|
| `biosphere-address-registry-v0.1.0.zip` | source tree |
| `biosphere-address-registry-v0.1.0.tar.gz` | source tree |
| `biosphere-address-registry-v0.1.0.bundle` | cloneable git history |
| `biosphere-address-registry-v0.1.0.SHA256SUMS` | checksums |

## Unpack

```bash
unzip biosphere-address-registry-v0.1.0.zip
cd biosphere-address-registry
python3 -m pip install -e ".[dev]"
python3 -m pytest -q
bar validate examples/organism-record.example.json
python3 scripts/run_local_conformance.py
```

## Restore git history from the bundle

```bash
git clone biosphere-address-registry-v0.1.0.bundle biosphere-address-registry
cd biosphere-address-registry
git checkout cursor/establish-biosphere-registry-887b
```

Then push to `https://github.com/sentry-bio/biosphere-address-registry` from an
account with write access. The Cursor GitHub App could not push there (403).

## What this does and does not claim

Done: constitution, v1 schemas, content-addressed records, local registry CLI,
synthetic linked example, local BAR-Core / BAR-Frame-candidate checks.

Not claimed: BAR-Map-real against the 234k v10.9 indexes, BAR-Certified, freeze
gate, or species pins. Follow `docs/END_TO_END_CONFORMANCE.md` on the internal
`:8100` stack for the real map run.
