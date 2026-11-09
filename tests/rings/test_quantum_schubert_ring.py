import pytest

from schubmult import QDSx, QSx, Sx, GeneratingSet, Permutation
from schubmult.rings.schubert.quantum_double_schubert_ring import QuantumDoubleSchubertRing, q_var
from schubmult.rings.schubert.quantum_schubert_ring import QuantumSingleSchubertRing
from schubmult.symbolic import expand, expand_func, sympify, sympy_Add
from schubmult.symbolic.poly.variables import CustomGeneratingSet, poly_genset

x, y, z, t = (GeneratingSet(name) for name in ("x", "y", "z", "t"))


@pytest.fixture(params=["double", "single"])
def ring_factory(request):
    if request.param == "double":
        return lambda **kwargs: QuantumDoubleSchubertRing(x, y, **kwargs)
    return lambda **kwargs: QuantumSingleSchubertRing(x, **kwargs)


@pytest.fixture(params=[t, CustomGeneratingSet([t[i] ** 2 for i in range(10)])])
def quantum_alphabet(request):
    return request.param


def replace_q(expr, alphabet):
    expr = sympify(expr)
    return expr.xreplace({s: alphabet[q_var.index(s)] for s in expr.free_symbols if s in q_var})


def assert_coefficients(actual, expected, alphabet):
    for w in actual.keys() | expected.keys():
        assert expand(actual.get(w, 0) - replace_q(expected.get(w, 0), alphabet)) == 0


def test_default_alphabet(ring_factory):
    ring = ring_factory()
    assert ring.q_var is q_var
    assert ring == ring_factory(q_var=q_var)
    assert hash(ring) == hash(ring_factory(q_var=q_var))
    assert QSx.q_var is q_var
    assert QDSx([]).ring.q_var is q_var


def test_custom_products_and_cache_isolation(ring_factory, quantum_alphabet):
    default, custom = ring_factory(), ring_factory(q_var=quantum_alphabet)
    assert custom.q_var is quantum_alphabet
    assert default != custom
    assert len({default, custom, ring_factory(q_var=quantum_alphabet)}) == 2
    for ring in (default, custom, default):
        actual = ring([2, 1]) * ring([2, 1])
        assert_coefficients(actual, default([2, 1]) * default([2, 1]), ring.q_var)
        assert expand(actual[Permutation([])] - ring.q_var[1]) == 0


def test_mixed_coefficient_alphabets(ring_factory, quantum_alphabet):
    default, custom = ring_factory(), ring_factory(q_var=quantum_alphabet)
    right = QuantumDoubleSchubertRing(x, z, q_var=quantum_alphabet)
    default_right = QuantumDoubleSchubertRing(x, z)
    left, other = custom([1, 3, 2]), right([2, 1])
    assert_coefficients(left * other, default([1, 3, 2]) * default_right([2, 1]), quantum_alphabet)
    assert_coefficients(other * left, default_right([2, 1]) * default([1, 3, 2]), quantum_alphabet)
    assert expand((left * other).expand() - left.expand() * other.expand()) == 0


def test_polynomial_multiplication(ring_factory, quantum_alphabet):
    default, custom = ring_factory(), ring_factory(q_var=quantum_alphabet)
    elem, reference = custom([2, 1]), default([2, 1])
    poly = (x[1] + x[2]) * x[1] ** 2
    for actual, expected in (
        (elem * poly, reference * poly),
        (poly * elem, poly * reference),
        (elem.mult_poly(poly), reference.mult_poly(poly)),
    ):
        assert_coefficients(actual, expected, quantum_alphabet)
        assert expand(actual.expand() - elem.expand() * poly) == 0


def test_expansion_conversion_and_substitution(ring_factory, quantum_alphabet):
    default, custom = ring_factory(), ring_factory(q_var=quantum_alphabet)
    elem, reference = custom([3, 2, 1]), default([3, 2, 1])
    assert expand(elem.expand() - replace_q(reference.expand(), quantum_alphabet)) == 0
    assert expand(elem.as_classical().expand() - elem.expand()) == 0
    assert expand(custom.from_expr(elem.expand()).expand() - elem.expand()) == 0
    substituted = elem.subs(x[1], x[2])
    assert substituted.ring == custom
    assert expand(substituted.expand() - elem.expand().subs(x[1], x[2])) == 0
    assert expand(custom.elem_sym_subs(2)[sympify("e_2_2")] - x[1] * x[2] - quantum_alphabet[1]) == 0


def test_symmetric_function_multiplication(ring_factory, quantum_alphabet):
    default, custom = ring_factory(), ring_factory(q_var=quantum_alphabet)
    elem = custom([2, 1])
    symmetric = custom.elem_sym(2, 2, [x[1], x[2]], poly_genset(0))
    reference = default.elem_sym(2, 2, [x[1], x[2]], poly_genset(0))
    assert symmetric != reference
    assert expand(expand_func(symmetric) - x[1] * x[2] - quantum_alphabet[1]) == 0
    assert expand(expand_func(reference) - x[1] * x[2] - q_var[1]) == 0
    assert expand(expand_func(sympy_Add(symmetric, reference)) - 2 * x[1] * x[2] - quantum_alphabet[1] - q_var[1]) == 0
    assert symmetric.func(*symmetric.args) is symmetric
    assert expand(expand_func(symmetric.xreplace({x[1]: x[3]})) - x[2] * x[3] - quantum_alphabet[2]) == 0
    assert_coefficients(custom.elem_mul(elem, symmetric), default.elem_mul(default([2, 1]), reference), quantum_alphabet)
    assert expand(expand_func(elem.in_SEM_basis()) - elem.expand()) == 0


def test_positive_products(quantum_alphabet):
    default = QuantumDoubleSchubertRing(x, y)
    custom = QuantumDoubleSchubertRing(x, y, q_var=quantum_alphabet)
    u = Permutation([2, 1])
    assert_coefficients(custom.cached_positive_product(u, u, custom), default.cached_positive_product(u, u, default), quantum_alphabet)
    assert expand(custom.cached_positive_product(u, u, custom)[Permutation([])] - quantum_alphabet[1]) == 0


def test_public_constructor_and_quantization(quantum_alphabet):
    assert QDSx([2, 1], "z", q_var=quantum_alphabet).ring == QuantumDoubleSchubertRing(x, z, q_var=quantum_alphabet)
    single = QuantumSingleSchubertRing(x, q_var=quantum_alphabet)
    poly = x[1] ** 2 * x[2]
    assert expand(single.quantize(poly) - replace_q(QSx.quantize(poly), quantum_alphabet)) == 0
    converted = single(Sx([3, 2, 1]))
    assert converted.ring == single
    assert expand(converted.expand() - poly) == 0
