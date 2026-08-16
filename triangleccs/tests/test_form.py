"""Form immutability and honesty tags."""

from __future__ import annotations

import pytest

from triangleccs.address import make_address
from triangleccs.datum.form import Form
from triangleccs.ledger import Tag


def test_form_hash_stable():
    a = Form()
    b = Form()
    assert a.form_hash == b.form_hash
    assert a.version.startswith("triangleccs-")


def test_form_rejects_non_h2():
    with pytest.raises(ValueError):
        Form(dim=3)  # type: ignore[arg-type]


def test_kappa_is_convention_not_state_equation():
    f = Form(kappa=1.25)
    s = f.summary()
    assert s["kappa_status"] == Tag.CONVENTION.value
    assert "LUCA" not in s["note"] or "not LUCA" in s["note"]
    # Changing kappa changes epoch hash
    assert Form(kappa=1.0).form_hash != f.form_hash


def test_address_tags_advisory_radius():
    f = Form()
    addr = make_address(form=f, theta=0.1, r=0.2, delta=0.01, resolvable=0.8, block_sep=1.0, residual=0.05)
    assert addr.tags["r"] == Tag.ADVISORY.value
    assert addr.tags["kappa"] == Tag.CONVENTION.value
    assert addr.theta_status == "candidate"
    assert addr.r_proxy == f.radial_proxy
    assert addr.delta == 0.01
    assert addr.resolvable == 0.8

