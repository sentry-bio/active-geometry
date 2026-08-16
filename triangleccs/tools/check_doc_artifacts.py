#!/usr/bin/env python3
"""Fail when tracked Markdown references a repo path that is not on disk."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
ENFORCED_ROOTS = ("docs", "triangleccs", "tests", "examples", "schemas", "tools")
PATH_EXTENSIONS = {
    ".py", ".md", ".json", ".sh", ".yaml", ".yml", ".toml", ".txt",
}
ALLOWLIST: set[str] = set()
BACKTICK = re.compile(r"`([^`]+)`")
EXTERNAL_PREFIXES = ("http://", "https://", "mailto:", "~", "/Users", "/home", "/tmp")


def looks_like_repo_path(token: str) -> bool:
    if not token or any(ch.isspace() for ch in token):
        return False
    if token.startswith(EXTERNAL_PREFIXES):
        return False
    if any(bad in token for bad in ("*", "…", "...", "://")):
        return False
    return Path(token).suffix in PATH_EXTENSIONS


def candidate_tokens(span: str) -> list[str]:
    tokens: list[str] = []
    for raw in re.split(r"[\s,;(){}\[\]]+", span.strip()):
        cleaned = raw.strip().strip("\"'`").rstrip(".,:;")
        if looks_like_repo_path(cleaned):
            tokens.append(cleaned)
    return tokens


def exists(token: str, doc: Path) -> bool:
    if token in ALLOWLIST:
        return True
    if (REPO_ROOT / token).exists():
        return True
    if (doc.parent / token).exists():
        return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    roots = [REPO_ROOT] if args.all else [REPO_ROOT / r for r in ENFORCED_ROOTS]
    missing: list[str] = []
    for root in roots:
        if not root.exists():
            continue
        for doc in root.rglob("*.md"):
            text = doc.read_text(encoding="utf-8")
            for span in BACKTICK.findall(text):
                for token in candidate_tokens(span):
                    if not exists(token, doc):
                        missing.append(f"{doc.relative_to(REPO_ROOT)}: `{token}`")
    if missing:
        print("Missing artifacts:")
        for m in missing:
            print(f"  {m}")
        return 1
    print("doc artifacts ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
