"""Public Address tuple — what consumers emit and compare."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal, Mapping

from triangleccs.datum.form import Form
from triangleccs.ledger import Tag, validate_tags


@dataclass(frozen=True)
class Address:
    theta: float
    r: float
    form_version: str
    form_hash: str
    r_proxy: str
    theta_status: Literal["candidate", "certified"] = "candidate"
    delta: float | None = None
    resolvable: float | None = None
    block_sep: float | None = None
    residual: float | None = None
    tags: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "tags", validate_tags(dict(self.tags)))
        if self.r < 0:
            raise ValueError("r must be non-negative")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def default_address_tags(
    *,
    theta_status: Literal["candidate", "certified"] = "candidate",
) -> dict[str, str]:
    theta_tag = (
        Tag.CANDIDATE.value
        if theta_status == "candidate"
        else Tag.EMPIRICAL.value
    )
    return {
        "theta": theta_tag,
        "r": Tag.ADVISORY.value,
        "kappa": Tag.CONVENTION.value,
        "dim": Tag.CONVENTION.value,
        "delta": Tag.INSTRUMENT.value,
        "resolvable": Tag.INSTRUMENT.value,
        "block_sep": Tag.INSTRUMENT.value,
        "residual": Tag.INSTRUMENT.value,
    }


def make_address(
    *,
    form: Form,
    theta: float,
    r: float,
    theta_status: Literal["candidate", "certified"] = "candidate",
    delta: float | None = None,
    resolvable: float | None = None,
    block_sep: float | None = None,
    residual: float | None = None,
    extra_tags: Mapping[str, str] | None = None,
) -> Address:
    tags = default_address_tags(theta_status=theta_status)
    if extra_tags:
        tags.update(validate_tags(extra_tags))
    return Address(
        theta=float(theta),
        theta_status=theta_status,
        r=float(r),
        r_proxy=form.radial_proxy,
        delta=None if delta is None else float(delta),
        resolvable=None if resolvable is None else float(resolvable),
        block_sep=None if block_sep is None else float(block_sep),
        residual=None if residual is None else float(residual),
        form_version=form.version,
        form_hash=form.form_hash,
        tags=tags,
    )
