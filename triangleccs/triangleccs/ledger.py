"""Honesty tags for every public quantity."""

from __future__ import annotations

from enum import Enum
from typing import Mapping


class Tag(str, Enum):
    THEOREM = "THEOREM"
    CONVENTION = "CONVENTION"
    INSTRUMENT = "INSTRUMENT"
    EMPIRICAL = "EMPIRICAL"
    CANDIDATE = "CANDIDATE"
    ADVISORY = "ADVISORY"
    CIRCULAR = "CIRCULAR"
    OVERLAY = "OVERLAY"


ALLOWED = frozenset(t.value for t in Tag)


def validate_tags(tags: Mapping[str, str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, value in tags.items():
        if value not in ALLOWED:
            raise ValueError(f"illegal tag for {key!r}: {value!r}")
        out[str(key)] = value
    return out
