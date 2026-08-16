"""Prefer the inner package over a vendored repo folder named triangleccs."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
