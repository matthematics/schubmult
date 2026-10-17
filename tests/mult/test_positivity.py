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


def _manifestly_positive(expr):
    # a minus sign is only allowed directly in front of a z-variable term
    import re

    return re.search(r"-\s*(\d+\*)?[y(]", str(expr)) is None


def test_posify_one_dominates_recursion_uses_shifted_alphabet():
    # Regression: u one-dominates w, so posify recurses on (phi1(u), v3, phi1(w)) over the
    # shifted alphabet var2[1:] (a plain tuple). Previously the reduced coefficient was computed
    # in the unshifted alphabet and compute_positive_rep crashed on tuple.index.
    u, v, w = Permutation([2, 4, 3, 1]), Permutation([1, 4, 3, 2]), Permutation([3, 4, 5, 1, 2])
    val = schubmult_double({u: 1}, v, y, z)[w]
    positive = posify(val, u, v, w, y, z)
    assert expand(sp(positive) - sp(val)) == 0
    assert _manifestly_positive(positive)


def test_compute_positive_rep_accepts_plain_tuple_alphabets():
    val = (y[2] - z[1]) + 2 * (y[3] - z[2])
    rep = compute_positive_rep(val, tuple(y[1:5]), tuple(z[1:5]))
    assert expand(sp(rep) - sp(val)) == 0


def test_posify_reconstructs_raw_value_for_s4():
    S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
    for u in S4:
        for v in S4:
            for w, val in schubmult_double({u: 1}, v, y, z).items():
                positive = posify(val, u, v, w, y, z)
                assert expand(sp(positive) - sp(val)) == 0, (u, v, w)
                assert _manifestly_positive(positive), (u, v, w, positive)
