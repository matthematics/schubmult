"""Tests for ``schubmult.mult._accel`` (the compiled ``schubmult_cpp`` dispatch layer)."""

import sys
import warnings

import sympy

from schubmult import Permutation
from schubmult.mult import _accel
from schubmult.mult.double import schubmult_double, schubmult_double_alt_from_elems, schubmult_double_from_elems
from schubmult.symbolic.common_polys import elem_sym_poly
from schubmult.abc import y, z
from schubmult.symbolic.poly.schub_poly import _vars

q = _vars.q_var


def test_available_in_normal_environment():
    assert _accel.available is True


def test_is_int():
    assert _accel._is_int(3) is True
    assert _accel._is_int(sympy.Integer(3)) is True
    assert _accel._is_int(sympy.Rational(1, 2)) is False
    assert _accel._is_int(sympy.Symbol("a")) is False


def test_call_falls_back_to_none_when_extension_raises_runtime_error():
    from schubmult import schubmult_cpp as cpp

    n = cpp.MAXN + 2
    big_v = list(range(n, 0, -1))
    assert _accel._call(cpp.schubmult_py, {tuple(range(1, n + 1)): 1}, big_v) is None


def test_schubmult_py_non_int_coefficient_returns_none():
    d = {Permutation([2, 1]): sympy.Symbol("a")}
    assert _accel.schubmult_py(d, [2, 1]) is None


def test_schubmult_py_int_coefficient_direct_call():
    u, v = Permutation([2, 1]), Permutation([2, 1])
    assert _accel.schubmult_py({u: 1}, v) == {Permutation([3, 1, 2]): 1}


def test_schubmult_double_direct_call():
    u, v = Permutation([2, 1, 3]), Permutation([1, 3, 2])
    assert _accel.schubmult_double({u: 1}, v, y, z) == schubmult_double({u: 1}, v, y, z)


def test_schubmult_q_fast_int_and_fractional_coefficients():
    u = v = Permutation([2, 1])
    int_res = _accel.schubmult_q_fast({u: 1}, v, q)
    frac_res = _accel.schubmult_q_fast({u: sympy.Rational(1, 2)}, v, q)
    assert frac_res == {w: c / 2 for w, c in int_res.items()}


def test_schubmult_q_fast_fractional_coefficient_none_on_oversized_permutation():
    from schubmult import schubmult_cpp as cpp

    n = cpp.MAXN + 2
    big_v = tuple(range(n, 0, -1))
    u = Permutation([2, 1])
    assert _accel.schubmult_q_fast({u: sympy.Rational(1, 2)}, big_v, q) is None


def test_schubmult_q_double_fast_direct_call():
    from schubmult.mult.quantum_double import schubmult_q_double_fast

    u, v = Permutation([2, 1]), Permutation([2, 1])
    assert _accel.schubmult_q_double_fast({u: 1}, v, y, z, q) == schubmult_q_double_fast({u: 1}, v, y, z, q)


def test_schubmult_double_from_elems_matches_schubmult_double():
    u, v = Permutation([2, 1, 3]), Permutation([1, 3, 2])
    d = {u: 1}
    assert schubmult_double_from_elems(d, v, y, z, elem_sym_poly) == schubmult_double(d, v, y, z)
    assert schubmult_double_alt_from_elems(d, v, y, z, elem_sym_poly) == schubmult_double(d, v, y, z)


def test_no_cpp_env_var_warns_and_disables_extension(monkeypatch):
    # reimport in-process (rather than via subprocess) so coverage instrumentation sees it
    import schubmult.mult as mult_pkg

    monkeypatch.setenv("SCHUBMULT_NO_CPP", "1")
    sys.modules.pop("schubmult.mult._accel", None)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            import schubmult.mult._accel as reimported
            assert reimported.available is False
            assert any("SCHUBMULT_NO_CPP" in str(w.message) for w in caught)
    finally:
        # undo the reimport's side effect of overwriting schubmult.mult._accel with the stub
        sys.modules["schubmult.mult._accel"] = _accel
        mult_pkg._accel = _accel
