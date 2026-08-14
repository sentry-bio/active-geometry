#!/usr/bin/env python3
"""Fail when a tracked Markdown file references a repo path that is not on disk.

This guards against the recurring failure mode in which documentation
describes an artifact (a script, a data file, a proof) that does not exist.
It scans backticked spans in tracked ``.md`` files, extracts tokens that look
like repository-relative paths, and checks that each one exists — resolved
against the repository root, the referring document's directory, or by
basename anywhere in the repository.

By default it enforces only the theory program's directories (``theory``,
``tools``, ``tests``), where correctness is guaranteed. ``--all`` scans the
whole repository as a non-fatal advisory, which surfaces pre-existing
data/figure references elsewhere without blocking.

A token is treated as a repository path when it has no URL scheme, is not an
absolute or home path outside the repo, contains no glob or ellipsis, and ends
in a known code/data extension. External references (URLs, DOIs, absolute
machine paths) and in-prose identifiers (``Metric.packingNumber``) are ignored.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
ENFORCED_ROOTS = ("theory", "tools", "tests")

PATH_EXTENSIONS = {
    ".py", ".lean", ".md", ".json", ".sh", ".yaml", ".yml",
    ".npz", ".npy", ".csv", ".tsv", ".txt", ".wl", ".ipynb", ".toml",
}

# Prescribed output filenames and placeholders — files a protocol tells the
# operator to *create*, not artifacts a document claims already exist.
ALLOWLIST: set[str] = {"preregistration.json"}
PLACEHOLDER_SUBSTRINGS = ("path/to", "your/", "example/")

BACKTICK = re.compile(r"`([^`]+)`")
EXTERNAL_PREFIXES = ("http://", "https://", "mailto:", "~", "/Users", "/home",
                     "/zfs", "/private", "/opt", "/tmp", "/var", "/etc")


def looks_like_repo_path(token: str) -> bool:
    if not token or any(ch.isspace() for ch in token):
        return False
    if token.startswith(EXTERNAL_PREFIXES):
        return False
    if any(bad in token for bad in ("*", "…", "...", "://")):
        return False
    if any(sub in token for sub in PLACEHOLDER_SUBSTRINGS):
        return False
    return Path(token).suffix in PATH_EXTENSIONS


def candidate_tokens(span: str) -> list[str]:
    tokens: list[str] = []
    for raw in re.split(r"[\s,;(){}\[\]]+", span.strip()):
        cleaned = raw.strip().strip("\"'`")
        cleaned = cleaned.rstrip(".,:;")
        cleaned = re.sub(r":\d+(?::\d+)?$", "", cleaned)
        if looks_like_repo_path(cleaned):
            tokens.append(cleaned)
    return tokens


def repo_basenames() -> set[str]:
    return {p.name for p in REPO_ROOT.rglob("*") if p.is_file()}


def make_existence_check(basenames: set[str]):
    def exists(token: str, doc: Path) -> bool:
        if token in ALLOWLIST:
            return True
        if (REPO_ROOT / token).exists():
            return True
        if (doc.parent / token).resolve().exists():
            return True
        return Path(token).name in basenames
    return exists


def tracked_markdown() -> list[Path]:
    try:
        out = subprocess.run(
            ["git", "ls-files", "*.md"],
            cwd=REPO_ROOT, capture_output=True, text=True, check=True,
        ).stdout
        files = [REPO_ROOT / line for line in out.splitlines() if line.strip()]
        if files:
            return files
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return list(REPO_ROOT.rglob("*.md"))


def in_enforced_scope(doc: Path) -> bool:
    rel = doc.relative_to(REPO_ROOT)
    return rel.parts and rel.parts[0] in ENFORCED_ROOTS


def scan(scope_all: bool) -> list[tuple[Path, str]]:
    exists = make_existence_check(repo_basenames())
    violations: list[tuple[Path, str]] = []
    for doc in tracked_markdown():
        if not scope_all and not in_enforced_scope(doc):
            continue
        try:
            text = doc.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        seen: set[str] = set()
        for span in BACKTICK.findall(text):
            for token in candidate_tokens(span):
                if token in seen:
                    continue
                seen.add(token)
                if not exists(token, doc):
                    violations.append((doc.relative_to(REPO_ROOT), token))
    return violations


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all", action="store_true",
        help="scan the whole repository as a non-fatal advisory",
    )
    args = parser.parse_args(argv)

    violations = scan(scope_all=args.all)
    scope = "repository (advisory)" if args.all else "theory/tools/tests"
    if violations:
        print(f"Dangling artifact references in {scope}:")
        for doc, token in violations:
            print(f"  {doc}: `{token}`")
        if args.all:
            print("(advisory mode: not failing)")
            return 0
        return 1
    print(f"check_doc_artifacts: all backticked paths in {scope} exist.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
