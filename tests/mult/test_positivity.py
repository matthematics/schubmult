"""Tests for ``schubmult.mult.positivity`` (manifestly positive double structure constants)."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.abc import x, y, z
from schubmult.mult.double import schubmult_double, schubmult_double_pair
from schubmult.mult.positivity import compute_positive_rep, dualcoeff, forwardcoeff, posify
from schubmult.symbolic import expand, sympify_sympy
from schubmult.symbolic.poly.schub_poly import schubpoly
from schubmult.utils.schub_lib import will_formula_work

S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
sp = sympify_sympy


def test_compute_positive_rep_reconstructs_known_polynomial():
    val = 2 * (y[1] - z[1]) + 3 * (y[2] - z[1]) * (y[1] - z[2])
    rep = compute_positive_rep(val, y, z)
    assert expand(sp(rep) - sp(val)) == 0


def test_compute_positive_rep_integer_passthrough():
    assert compute_positive_rep(5, y, z) == 5


def test_posify_reconstructs_raw_value_for_s3():
    for u in S3:
        for v in S3:
            coeff_dict = schubmult_double({u: 1}, v, y, z)
            for w, val in coeff_dict.items():
                positive = posify(val, u, v, w, y, z)
                assert expand(sp(positive) - sp(val)) == 0, (u, v, w)


def test_dualcoeff_identity_u_matches_schubpoly():
    u = Permutation([])
    for v in S3:
        coeff_dict = schubmult_double_pair(u, v, y, z)
        for w, val in coeff_dict.items():
            vp = v * (~w)
            if vp.inv == v.inv - w.inv:
                expected = sp(schubpoly(vp, y, z))
                assert expand(sp(dualcoeff(u, v, w, y, z)) - expected) == 0, (v, w)


def test_forwardcoeff_matches_raw_value_when_formula_applies():
    found = 0
    for u in S3:
        for v in S3:
            if not will_formula_work(u, v):
                continue
            coeff_dict = schubmult_double_pair(u, v, y, z)
            for w, val in coeff_dict.items():
                found += 1
                assert expand(sp(forwardcoeff(u, v, w, y, z)) - sp(val)) == 0, (u, v, w)
    assert found > 0
