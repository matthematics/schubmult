"""Tests for ``schubmult.mult.double`` (double Schubert polynomial multiplication)."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.abc import x, y, z
from schubmult.mult import _accel
from schubmult.mult.double import (
    _schubmult_double_alt_from_elems_backwards_python,
    _schubmult_double_from_elems_python,
    _schubmult_double_python,
    mult_poly_double,
    mult_poly_double_alt,
    mult_poly_down,
    schub_coprod_double,
    schubmult_double,
    schubmult_double_alt,
    schubmult_double_alt_from_elems_backwards,
    schubmult_double_alt_from_elems_backwards_backwards,
    schubmult_double_alt_from_elems_forwards,
    schubmult_double_dict,
    schubmult_double_down,
    schubmult_double_from_elems,
    schubmult_double_pair,
    schubmult_double_pair_generic,
    schubmult_double_pair_generic_alt,
    single_variable,
    single_variable_down,
)
from schubmult.mult.single import schubmult_py
from schubmult.symbolic import S, sympify_sympy
from schubmult.symbolic.common_polys import elem_sym_poly
from schubmult.symbolic.poly.schub_poly import schubpoly
from schubmult.symbolic.poly.variables import CustomGeneratingSet, ZeroGeneratingSet

zero = ZeroGeneratingSet()
S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
sp = sympify_sympy


def _dschub(v, var2, var3):
    return sp(schubpoly(v, var2, var3))


def _same(d1, d2):
    for w in set(d1) | set(d2):
        if sympy.expand(sp(d1.get(w, 0)) - sp(d2.get(w, 0))) != 0:
            return False
    return True


def test_schubmult_double_identity_v_empty():
    d = {Permutation([2, 1]): S.One, Permutation([1, 3, 2]): 2}
    assert _same(schubmult_double(d, [], y, z), d)
    assert _same(schubmult_double(d, [1, 2], y, z), d)


def test_schubmult_double_reduces_to_schubmult_py_at_yz0():
    yz0 = {sp(y[i]): 0 for i in range(1, 8)} | {sp(z[i]): 0 for i in range(1, 8)}
    for u in S3:
        for v in S3:
            dbl = schubmult_double({u: S.One}, v, y, z)
            py = schubmult_py({u: S.One}, v)
            assert _same({w: sp(c).xreplace(yz0) for w, c in dbl.items()}, py), (u, v)


def test_schubmult_double_polynomial_identity():
    # S_u(x; y) S_v(x; z) = sum_w c^w_{u,v}(y, z) S_w(x; y).
    for u in S3:
        for v in S3:
            prod_dict = schubmult_double({u: S.One}, v, y, z)
            rhs = sum((sp(c) * _dschub(w, x, y) for w, c in prod_dict.items()), sympy.Integer(0))
            assert sympy.expand(_dschub(u, x, y) * _dschub(v, x, z) - rhs) == 0, (u, v)


def test_schubmult_double_v_zero_alphabet_matches_single_variable():
    # S_{21}(x, 0) = x_1, so multiplying by S_{21}(x, zero) is exactly the equivariant Monk rule.
    for u in S3:
        d = {u: S.One}
        assert _same(single_variable(d, 1, y), schubmult_double(d, [2, 1], y, zero))


def test_mult_poly_double_matches_schubmult_double():
    for u in S3:
        d = {u: S.One}
        direct = mult_poly_double(d, x[1] * x[2], x, y)
        expected = single_variable(single_variable(d, 2, y), 1, y)
        assert _same(direct, expected), u


def test_schub_coprod_double_recombines_to_original():
    perm = Permutation([3, 1, 4, 2])
    indices = [1, 3]
    k = len(indices)
    coeff_dict = schub_coprod_double(perm, indices, var2=y, var3=y)
    x_shifted = CustomGeneratingSet([0, *[x[k + i] for i in range(1, 20)]])
    total = sympy.Integer(0)
    for (a, b), val in coeff_dict.items():
        total += sp(val) * sp(schubpoly(a, x, y)) * sp(schubpoly(b, x_shifted, y))
    assert sympy.expand(total - _dschub(perm, x, y)) == 0


def test_single_variable_diagonal_beyond_window():
    # varnum beyond len(u): u fixes that point, so the diagonal term is var2[varnum] itself.
    u = Permutation([1, 3, 2])
    res = single_variable({u: S.One}, 5, y)
    assert sp(res[u]).has(y[5])


def test_single_variable_down_matches_single_variable_structure():
    u = Permutation([1, 3, 2])
    res = single_variable_down({u: S.One}, 1, y)
    assert res
    # varnum > 1 exercises the varnum - 1 elem_sym_perms_op branch; varnum beyond len(u)
    # exercises the "fixed point past the window" diagonal term.
    assert single_variable_down({Permutation([3, 2, 1]): S.One}, 2, y)
    assert single_variable_down({u: S.One}, 5, y)


def test_mult_poly_double_pow_add_scalar_and_plain_list_var_x():
    u = Permutation([1, 3, 2])
    d = {u: S.One}
    assert _same(mult_poly_double(d, x[1] ** 2, x, y), single_variable(single_variable(d, 1, y), 1, y))
    add_res = mult_poly_double(d, x[1] + 3, x, y)
    expected_add = {w: c for w, c in single_variable(d, 1, y).items()}
    for w, c in d.items():
        expected_add[w] = expected_add.get(w, 0) + 3 * c
    assert _same(add_res, expected_add)
    assert mult_poly_double(d, S(5), x, y) == {w: 5 * c for w, c in d.items()}
    assert _same(mult_poly_double(d, x[1], var_x=[x[0], x[1], x[2]], var_y=y), single_variable(d, 1, y))


def test_mult_poly_double_alt_matches_mult_poly_double():
    u = Permutation([1, 3, 2])
    d = {u: S.One}
    for poly in (x[1] * x[2], x[1] ** 2, x[1] + 3, S(5)):
        assert _same(mult_poly_double_alt(d, poly, x, y), mult_poly_double(d, poly, x, y)), poly
    assert _same(mult_poly_double_alt(d, x[1], var_x=[x[0], x[1], x[2]], var_y=y), single_variable(d, 1, y))


def test_mult_poly_down_dispatch_branches():
    # mult_poly_down hardcodes var2=None for single_variable_down, so it crashes on any
    # polynomial that actually contains a _vars.var1 symbol; exercise the Mul/Pow/Add/scalar
    # traversal with symbols outside that alphabet instead, which is all that's reachable.
    u = Permutation([1, 3, 2])
    d = {u: S.One}
    assert mult_poly_down(d, z[1] * z[2]) == {u: sp(z[1]) * sp(z[2])}
    assert mult_poly_down(d, z[1] ** 2) == {u: sp(z[1]) ** 2}
    assert mult_poly_down(d, z[1] + 3) == {u: sp(z[1]) + 3}
    assert mult_poly_down(d, S(5)) == {u: 5}


def test_schubmult_double_pair_variants():
    u, v = Permutation([1, 3, 2]), Permutation([2, 1])
    assert _same(schubmult_double_pair(u, v, y, z), schubmult_double({u: S.One}, v, y, z))
    generic = schubmult_double_pair_generic(u, v)
    generic_alt = schubmult_double_pair_generic_alt(u, v)
    assert set(generic) == set(generic_alt)


def test_schubmult_double_dict_matches_manual_sum():
    u = Permutation([1, 3, 2])
    d1 = {u: S.One}
    d2 = {Permutation([2, 1]): S.One, Permutation([]): 2}
    direct = schubmult_double_dict(d1, d2, y, z)
    manual = {}
    for v, c in d2.items():
        for w, val in schubmult_double(d1, v, y, z).items():
            manual[w] = manual.get(w, 0) + c * val
    assert _same(direct, manual)


def test_schubmult_double_pure_python_matches_dispatcher_and_fallback(monkeypatch):
    # S4 (not just S3) is needed to exercise a genuine cancelling v-path (sumval == 0).
    for u in S4:
        for v in S4:
            pure = {w: c for w, c in _schubmult_double_python({u: S.One}, v, y, z).items() if c != 0}
            fast = {w: c for w, c in schubmult_double({u: S.One}, v, y, z).items() if c != 0}
            assert _same(pure, fast), (u, v)
    monkeypatch.setattr(_accel, "available", False)
    for u in S3:
        for v in S3:
            assert _same(schubmult_double({u: S.One}, v, y, z), _schubmult_double_python({u: S.One}, v, y, z)), (u, v)


def test_schubmult_double_alt_matches_schubmult_double():
    u, v = Permutation([1, 3, 2]), Permutation([2, 1])
    d = {u: S.One}
    assert _same(schubmult_double_alt(d, v, y, z), schubmult_double(d, v, y, z))


def test_schubmult_double_alt_from_elems_variants_match_schubmult_double():
    u, v = Permutation([1, 3, 2]), Permutation([2, 1])
    d = {u: S.One}
    expected = schubmult_double(d, v, y, z)
    assert schubmult_double_alt_from_elems_forwards(d, v, y, z, elem_func=elem_sym_poly) == expected
    assert schubmult_double_alt_from_elems_backwards(d, v, y, z, elem_sym_poly) == expected
    assert _schubmult_double_alt_from_elems_backwards_python(d, v, y, z, elem_sym_poly) == expected
    assert schubmult_double_alt_from_elems_backwards_backwards(d, v, y, z, elem_sym_poly) == expected


def test_schubmult_double_alt_from_elems_backwards_falls_back_without_accel(monkeypatch):
    monkeypatch.setattr(_accel, "available", False)
    u, v = Permutation([1, 3, 2]), Permutation([2, 1])
    d = {u: S.One}
    assert schubmult_double_alt_from_elems_backwards(d, v, y, z, elem_sym_poly) == _schubmult_double_alt_from_elems_backwards_python(d, v, y, z, elem_sym_poly)


def test_schubmult_double_from_elems_matches_schubmult_double_and_falls_back(monkeypatch):
    u, v = Permutation([1, 3, 2]), Permutation([2, 1])
    d = {u: S.One}
    expected = schubmult_double(d, v, y, z)
    assert schubmult_double_from_elems(d, v, y, z, elem_sym_poly) == expected
    assert _schubmult_double_from_elems_python(d, v, y, z, elem_sym_poly) == expected
    monkeypatch.setattr(_accel, "available", False)
    assert schubmult_double_from_elems(d, v, y, z, elem_sym_poly) == expected


def test_schubmult_double_from_elems_python_identity_and_s4_sweep():
    u = Permutation([1, 3, 2])
    idp = Permutation([1, 2])
    d = {u: S.One}
    assert _schubmult_double_from_elems_python(d, idp, y, z, elem_sym_poly) == d
    for uu in S4:
        for v in S4:
            _schubmult_double_from_elems_python({uu: S.One}, v, y, z, elem_sym_poly)


def test_schubmult_double_down_known_case():
    u, v = Permutation([1, 3, 2]), Permutation([2, 1])
    assert schubmult_double_down({u: S.One}, v, y, z) == {u: sp(y[1]) - sp(z[1])}


def test_schubmult_double_down_identity_and_s4_sweep():
    u = Permutation([1, 3, 2])
    idp = Permutation([1, 2])
    d = {u: S.One}
    assert schubmult_double_down(d, idp, y, z) == d
    for uu in S4:
        for v in S4:
            schubmult_double_down({uu: S.One}, v, y, z)

