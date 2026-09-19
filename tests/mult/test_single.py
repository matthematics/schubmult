"""Tests for ``schubmult.mult.single`` (ordinary Schubert polynomial multiplication)."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.mult import _accel
from schubmult.mult.single import _schubmult_py_python, mult_poly_py, schub_coprod_py, schubmult_py, schubmult_py_down, single_variable
from schubmult.symbolic import S, sympify_sympy
from schubmult.symbolic.poly.schub_poly import schubpoly
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet, ZeroGeneratingSet

x = GeneratingSet("x")
zero = ZeroGeneratingSet()
S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
sp = sympify_sympy


def _schub(v):
    return sp(schubpoly(v, x, zero))


def _same(d1, d2):
    for w in set(d1) | set(d2):
        if sympy.expand(sp(d1.get(w, 0)) - sp(d2.get(w, 0))) != 0:
            return False
    return True


def test_schubmult_py_identity_v_empty():
    d = {Permutation([2, 1]): 3, Permutation([1, 3, 2]): 1}
    assert schubmult_py(d, []) == d
    assert schubmult_py(d, [1, 2]) == d


def test_schubmult_py_known_monk_value():
    # S_{21} * S_{21} = S_{312}: the only length-2 permutation covering [2,1] reachable
    # by crossing position 1 twice via the Monk/v-path recursion.
    res = schubmult_py({Permutation([2, 1]): 1}, [2, 1])
    assert res == {Permutation([3, 1, 2]): 1}


def test_schubmult_py_polynomial_identity():
    for u in S3:
        for v in S3:
            prod_dict = schubmult_py({u: S.One}, v)
            rhs = sum((c * _schub(w) for w, c in prod_dict.items()), sympy.Integer(0))
            assert sympy.expand(_schub(u) * _schub(v) - rhs) == 0, (u, v)


def test_schubmult_py_matches_polynomial_identity_s4():
    u, v = Permutation([2, 1, 4, 3]), Permutation([3, 4, 1, 2])
    prod_dict = schubmult_py({u: S.One}, v)
    rhs = sum((c * _schub(w) for w, c in prod_dict.items()), sympy.Integer(0))
    assert sympy.expand(_schub(u) * _schub(v) - rhs) == 0


def test_single_variable_matches_schubmult_py_with_s1():
    # S_{x_k} * S_u = single_variable(coeff_dict, k); at k = 1 this is S_{21} * S_u
    for u in S3:
        d = {u: S.One}
        assert _same(single_variable(d, 1), schubmult_py(d, [2, 1]))


def test_mult_poly_py_matches_single_variable_chain():
    for u in S3:
        d = {u: S.One}
        direct = mult_poly_py(d, x[1] * x[2])
        expected = single_variable(single_variable(d, 2), 1)
        assert _same(direct, expected), u


def test_schubmult_py_down_identity_v_empty():
    d = {Permutation([2, 1]): 3, Permutation([1, 3, 2]): 1}
    assert schubmult_py_down(d, [1, 2]) == d


def test_schubmult_py_down_nontrivial_v():
    # no dedicated ground truth for the down variant; just exercise the full v-path loop
    # (including cancellation to a zero partial sum, only reachable for S4-sized inputs)
    # and check it stays a valid (nonempty) coefficient dict.
    res = schubmult_py_down({Permutation([3, 2, 1]): S.One}, [2, 3, 1])
    assert res
    for u in S4:
        for v in S4:
            schubmult_py_down({u: S.One}, v)


def test_schubmult_py_pure_python_matches_dispatcher():
    # the pure-Python kernel can leave explicit zero entries (from cancelling v-paths)
    # that the C++ kernel omits, so compare with zeros stripped.
    for u in S4:
        for v in S4:
            pure = {w: c for w, c in _schubmult_py_python({u: S.One}, v).items() if c != 0}
            fast = {w: c for w, c in schubmult_py({u: S.One}, v).items() if c != 0}
            assert pure == fast, (u, v)


def test_schubmult_py_falls_back_to_pure_python_when_accel_unavailable(monkeypatch):
    monkeypatch.setattr(_accel, "available", False)
    for u in S3:
        for v in S3:
            assert schubmult_py({u: S.One}, v) == _schubmult_py_python({u: S.One}, v), (u, v)


def test_mult_poly_py_pow_add_and_scalar_branches():
    d = {Permutation([2, 1]): S.One}
    pow_res = mult_poly_py(d, x[1] ** 2)
    expected_pow = single_variable(single_variable(d, 1), 1)
    assert pow_res == expected_pow

    add_res = mult_poly_py(d, x[1] + 3)
    expected_add = {w: c for w, c in single_variable(d, 1).items()}
    for w, c in d.items():
        expected_add[w] = expected_add.get(w, 0) + 3 * c
    assert add_res == expected_add

    scalar_res = mult_poly_py(d, S(5))
    assert scalar_res == {w: 5 * c for w, c in d.items()}


def test_mult_poly_py_accepts_plain_list_var_x():
    d = {Permutation([2, 1]): S.One}
    assert mult_poly_py(d, x[1], var_x=[x[0], x[1], x[2]]) == single_variable(d, 1)


def test_schub_coprod_py_recombines_to_original():
    # sum_{(a,b)} coeff * S_a(x_1..x_k) S_b(x_{k+1}..) = S_perm(x), after shifting b's
    # variables up by k (the split introduced by schub_coprod_py).
    perm = Permutation([3, 1, 4, 2])
    indices = [1, 3]
    k = len(indices)
    coeff_dict = schub_coprod_py(perm, indices)
    total = sympy.Integer(0)
    x_shifted = CustomGeneratingSet([0, *[x[k + i] for i in range(1, 20)]])
    for (a, b), val in coeff_dict.items():
        b_shifted_poly = sp(schubpoly(b, x_shifted, zero))
        total += val * sp(schubpoly(a, x, zero)) * b_shifted_poly
    assert sympy.expand(total - _schub(perm)) == 0
