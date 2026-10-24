import pytest

from schubmult import DGx, QDGx, GeneratingSet, Permutation, QuantumDoubleGrothendieckElement, QuantumDoubleGrothendieckRing
from schubmult.abc import beta, x, y, z
from schubmult.mult.groth_quantum import grothmult_q
from schubmult.mult.groth_quantum_double import grothmult_q_double
from schubmult.symbolic import S, Symbol, expand, sympify, sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars

q = _vars.q_var


def assert_zero(expr):
    assert sympify_sympy(expr).cancel() == 0


def assert_same_coefficients(actual, expected):
    assert all(sympify_sympy(actual.get(w, 0) - expected.get(w, 0)).cancel() == 0 for w in actual.keys() | expected.keys())


def test_public_constructors():
    from schubmult.rings.schubert import QDGx as public_constructor
    from schubmult.rings.schubert import QuantumDoubleGrothendieckRing as public_ring

    assert public_constructor is QDGx
    assert public_ring is QuantumDoubleGrothendieckRing
    elem = QDGx([2, 1])
    assert isinstance(elem, QuantumDoubleGrothendieckElement)
    assert elem.ring.beta == beta
    assert elem == QDGx(Permutation([2, 1]))
    assert elem.ring(elem) is elem
    assert elem.ring(0) == elem.ring.zero
    assert elem.ring(1) == elem.ring.one
    assert elem.ring.from_dict({(2, 1): 1, (): 0}) == elem
    assert elem.as_quantum() is elem
    assert "QuantumDoubleGrothendieckRing" in repr(elem.ring)
    assert "Quantum double Grothendieck" in str(elem.ring)


@pytest.mark.parametrize("value", [None, -1, 0, sympify("1/2"), Symbol("t") + 1])
def test_multiplication_uses_kernel(value, monkeypatch):
    import schubmult.rings.schubert.quantum_double_grothendieck_ring as module

    left, right = QDGx([1, 3, 2], beta=value), QDGx([2, 1], beta=value)
    expected = grothmult_q_double(left, Permutation([2, 1]), y, y, beta=left.ring.beta)

    def forbid_expansion(*args, **kwargs):
        pytest.fail("Multiplication must not expand polynomials or solve for coefficients")

    monkeypatch.setattr(module, "qgroth_poly", forbid_expansion)
    assert_same_coefficients(left * right, expected)


@pytest.mark.parametrize("value", [-1, 0, None])
def test_explicit_polynomial_product_identity(value):
    left, right = QDGx([1, 3, 2], beta=value), QDGx([2, 1], beta=value)
    product = left * right
    assert_zero(left.expand() * right.expand() - product.expand())
    assert_same_coefficients(product.simplify(), product)
    assert_same_coefficients(product.simplify(factor=False), product)


def test_s1_expansion():
    expected = x[1] + y[1] + beta * x[1] * y[1] - beta * q[1] * (1 + beta * x[1]) * (1 + beta * y[1])
    assert_zero(QDGx([2, 1]).expand() - expected)


def test_mixed_alphabets():
    left, right = QDGx([1, 3, 2], beta=-1), QDGx([2, 1], "z", beta=-1)
    product = left * right
    assert product.ring == left.ring
    assert_zero(left.expand() * right.expand() - product.expand())
    assert_same_coefficients(product, grothmult_q_double(left, Permutation([2, 1]), y, z, beta=-1))
    assert_same_coefficients(right * left, grothmult_q_double(right, Permutation([1, 3, 2]), z, y, beta=-1))


def test_zero_alphabet():
    left, right = QDGx([1, 3, 2], "0", beta=-1), QDGx([2, 1], "0", beta=-1)
    assert_same_coefficients(left * right, grothmult_q(left, Permutation([2, 1]), beta=-1))
    assert_zero(left.expand() * right.expand() - (left * right).expand())
    assert not any(symbol in y for symbol in left.expand().free_symbols)


def test_zero_alphabet_rings_are_compatible():
    left, right = QDGx([2, 1], "0"), QDGx([1, 3, 2], "0")
    assert left.ring == right.ring
    assert left + right == left.ring.from_dict({Permutation([2, 1]): 1, Permutation([1, 3, 2]): 1})
    assert left - right == left.ring.from_dict({Permutation([2, 1]): 1, Permutation([1, 3, 2]): -1})
    assert left.ring(right) == right


def test_classical_specialization():
    left, right = QDGx([1, 3, 2], beta=-1), QDGx([2, 1], beta=-1)
    product = (left * right).simplify()
    classical = DGx([1, 3, 2]).ring
    from schubmult.rings.schubert.double_grothendieck_ring import DoubleGrothendieckRing

    classical = DoubleGrothendieckRing(classical.genset, classical.coeff_genset, beta=-1)
    expected = classical([1, 3, 2]) * classical([2, 1])
    actual = {w: c.subs({q[i]: 0 for i in range(1, 5)}) for w, c in product.items()}
    assert_same_coefficients(actual, expected)


def test_linear_combinations_and_scalars():
    a, b = QDGx([2, 1], beta=-1), QDGx([1, 3, 2], beta=-1)
    scalar = 1 + y[1] + q[1]
    assert a * a.ring.one == a
    assert a.ring.one * a == a
    assert a * a.ring.zero == a.ring.zero
    assert a ** 0 == a.ring.one
    assert_same_coefficients((a + b) * a, a * a + b * a)
    assert a * scalar == scalar * a
    assert a.mult_poly(scalar) == scalar * a
    assert_zero((scalar * a + 2 - b).expand() - (scalar * a.expand() + 2 - b.expand()))
    assert (a - a).simplify() == a.ring.zero
    assert a.expand(deep=False).ring == a.ring


@pytest.mark.parametrize("mixed", [False, True])
def test_associativity(mixed):
    a, b = QDGx([2, 1], beta=-1), QDGx([1, 3, 2], beta=-1)
    c = QDGx([2, 1], "z" if mixed else "y", beta=-1)
    assert_same_coefficients((a * b) * c, a * (b * c))


def test_stable_product_retains_terms_outside_s4():
    product = QDGx([1, 2, 4, 3], beta=-1) * QDGx([3, 2, 4, 1], beta=-1)
    assert Permutation([3, 2, 5, 1, 4]) in product
    assert Permutation([3, 2, 1, 5, 4]) in product
    assert_zero(product[Permutation([3, 2, 1])] - q[3] * (1-y[1]) / (1-y[4]))


def test_simplify_cancels_hidden_zeros():
    ring = QDGx([], beta=-1).ring
    hidden = (y[1] - y[2]) * (1 - y[3]) - (y[1] - y[2] - y[1] * y[3] + y[2] * y[3])
    elem = ring.from_dict({Permutation([2, 1]): hidden, Permutation([]): S.One})
    assert elem.simplify() == ring.one


@pytest.mark.parametrize("operation", ["add", "sub", "mul", "convert"])
@pytest.mark.parametrize("other", [QDGx([2, 1], beta=0), DGx([2, 1])])
def test_incompatible_rings_are_rejected(operation, other):
    elem = QDGx([2, 1], beta=-1)
    with pytest.raises(ValueError):
        if operation == "add":
            elem + other
        elif operation == "sub":
            elem - other
        elif operation == "mul":
            elem * other
        else:
            elem.ring(other)


def test_mixed_alphabet_addition_is_rejected():
    with pytest.raises(ValueError):
        QDGx([2, 1]) + QDGx([2, 1], "z")


def test_different_x_alphabets_are_rejected():
    other = QuantumDoubleGrothendieckRing(GeneratingSet("a"), y)([2, 1])
    with pytest.raises(ValueError):
        QDGx([2, 1]) * other


@pytest.mark.parametrize("operation", ["construct", "right_multiply", "left_multiply", "mult_poly"])
def test_unsupported_x_polynomial_conversion_is_explicit(operation):
    elem = QDGx([2, 1])
    with pytest.raises(NotImplementedError, match="x-polynomials"):
        if operation == "construct":
            elem.ring(x[1])
        elif operation == "right_multiply":
            elem * x[1]
        elif operation == "left_multiply":
            x[1] * elem
        else:
            elem.mult_poly(x[1])


def test_ring_equality_and_hash_include_beta():
    a, b = QDGx([], beta=-1).ring, QDGx([], beta=-1).ring
    assert a == b
    assert hash(a) == hash(b)
    assert a != QDGx([], beta=0).ring
    assert a != QDGx([], "z", beta=-1).ring
    assert a != DGx([]).ring


def test_printing():
    from schubmult.symbolic import latex, sstr

    elem = QDGx([2, 1])
    assert "QDGx" in str(elem)
    assert "QDGx" in sstr(elem.as_expr())
    assert "Q" in latex(elem)
    assert str(elem.ring.zero) == "0"
    assert elem.as_expr() != sympify_sympy(DGx([2, 1]).as_terms()[0])
    assert elem.as_expr() != QDGx([1, 3, 2]).as_expr()
    assert len((elem + QDGx([1, 3, 2])).as_expr().args) == 2


def test_classical_basis_conversion():
    elem = QDGx([2, 1], beta=-1)
    classical = elem.as_classical()
    assert classical.ring.beta == -1
    assert_zero(classical.expand() - elem.expand())


def test_fourth_appendix_product_is_an_explicit_polynomial_identity():
    left, right = QDGx([2, 3, 4, 1], beta=-1), QDGx([3, 2, 1], beta=-1)
    actual = (left * right).simplify(factor=False)
    expected = {
        Permutation([4, 3, 2, 1]): (1-y[1])**2 / ((1-y[2])*(1-y[3])),
        Permutation([3, 4, 2, 1]): (y[1]-y[3])*(1-y[1]) / ((1-y[2])*(1-y[3])),
        Permutation([4, 2, 3, 1]): (y[1]-y[2])*(1-y[1]) / ((1-y[2])*(1-y[3])),
        Permutation([3, 2, 4, 1]): (y[2]-y[1])*(y[3]-y[1]) / ((1-y[2])*(1-y[3])),
    }
    assert_same_coefficients(actual, expected)
    assert_zero(expand(left.expand() * right.expand() - actual.expand()))
