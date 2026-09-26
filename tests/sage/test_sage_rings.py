"""Property tests for the SageMath integration (``schubmult.sage``).

Skipped entirely when Sage is not importable; CI runs them in the ``sage-doctests`` job. Every product
is checked against polynomial multiplication after ``expand()``, which is independent of the kernels
under test (Sage does the polynomial arithmetic).
"""

import itertools

import pytest

sage_all = pytest.importorskip("sage.all")

from sage.all import QQ, ZZ, Permutation, Permutations, PolynomialRing, SchubertPolynomialRing, SymmetricFunctions  # noqa: E402
from sage.misc.sage_unittest import TestSuite as SageTestSuite  # noqa: E402  (aliased so pytest does not try to collect it)

from schubmult.sage import DoubleSchubertPolynomialRing, QuantumDoubleSchubertPolynomialRing, QuantumSchubertPolynomialRing  # noqa: E402

S3 = list(Permutations(3))
S4_SAMPLE = [Permutation(p) for p in ([2, 1, 4, 3], [3, 1, 4, 2], [1, 4, 2, 3], [4, 2, 3, 1], [2, 4, 1, 3])]


def common(*polys):
    """Coerce polynomials from different finite rings into one ring on the union of their variables."""
    names = sorted({str(g) for p in polys for g in p.parent().gens()})
    T = PolynomialRing(QQ, names)
    return [T(p) for p in polys]


def product_matches(elem_product, *factors):
    """``elem_product.expand()`` equals the product of the factors' expansions, in a common ring."""
    lhs, *fs = common(elem_product.expand(), *(f.expand() for f in factors))
    rhs = fs[0]
    for f in fs[1:]:
        rhs *= f
    return lhs == rhs


def zero_out(poly, prefix):
    """Set every variable whose name starts with ``prefix`` to zero."""
    return poly.subs({g: 0 for g in poly.parent().gens() if str(g).startswith(prefix)})


@pytest.fixture(scope="module")
def X():
    return DoubleSchubertPolynomialRing(QQ)


@pytest.fixture(scope="module")
def Z():
    return DoubleSchubertPolynomialRing(QQ, "z")


@pytest.fixture(scope="module")
def S():
    return SchubertPolynomialRing(QQ)


@pytest.fixture(scope="module")
def Q():
    return QuantumSchubertPolynomialRing(QQ)


@pytest.fixture(scope="module")
def QD():
    return QuantumDoubleSchubertPolynomialRing(QQ)


# --- Sage category framework -----------------------------------------------------------------------


@pytest.mark.parametrize(
    "make",
    [
        lambda: DoubleSchubertPolynomialRing(ZZ),
        lambda: DoubleSchubertPolynomialRing(QQ, "z"),
        lambda: QuantumSchubertPolynomialRing(QQ),
        lambda: QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3)),
        lambda: QuantumDoubleSchubertPolynomialRing(QQ),
        lambda: QuantumDoubleSchubertPolynomialRing(QQ, parabolic=(1, 2)),
    ],
    ids=["double_ZZ", "double_z", "quantum", "quantum_parabolic", "quantum_double", "quantum_double_parabolic"],
)
def test_sage_testsuite(make):
    """Sage's generic parent/element/algebra axioms (associativity, distributivity, pickling, coercion)."""
    SageTestSuite(make()).run(raise_on_failure=True)


def test_unique_representation():
    assert DoubleSchubertPolynomialRing(QQ) is DoubleSchubertPolynomialRing(QQ, "y", ("z", "y"))
    assert QuantumSchubertPolynomialRing(QQ, parabolic=[2, 3]) is QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3))


# --- double Schubert -------------------------------------------------------------------------------


def test_double_specializes_to_single(X, S):
    """S_w(x; 0) = S_w(x)."""
    for w in S3 + S4_SAMPLE:
        lhs, rhs = common(zero_out(X(w).expand(), "y"), S(w).expand())
        assert lhs == rhs, w


def test_double_product_matches_polynomial_product(X):
    for u, v in itertools.product(S3, S3):
        assert product_matches(X(u) * X(v), X(u), X(v)), (u, v)


def test_double_product_matches_polynomial_product_S4_sample(X):
    for u, v in itertools.product(S4_SAMPLE, S3):
        assert product_matches(X(u) * X(v), X(u), X(v)), (u, v)


def test_mixed_alphabet_product(X, Z):
    """S_u(x; y) S_v(x; z) via coercion of the z-ring into the y-ring."""
    for u, v in itertools.product(S3, S3):
        assert product_matches(X(u) * Z(v), X(u), Z(v)), (u, v)


def test_polynomial_round_trip(X, Z):
    for u, v in itertools.product(S3, S3):
        p = X(u) * X(v)
        assert X(p.expand()) == p, (u, v)
        m = X(u) * Z(v)
        assert X(m.expand()) == m, (u, v)


def test_other_alphabet_round_trip(X, Z):
    for w in S3 + S4_SAMPLE:
        assert Z(X(Z(w))) == Z(w), w
        assert X(Z(X(w))) == X(w), w


def test_single_schubert_coerces(X, S):
    for w in S3 + S4_SAMPLE:
        lhs, rhs = common(X(S(w)).expand(), S(w).expand())
        assert lhs == rhs, w
    assert (X([2, 1]) + S([2, 1])).parent() is X


def test_divided_difference_is_the_operator(X):
    """partial_i f = (f - s_i f) / (x_{i-1} - x_i) on the x variables, coefficients untouched."""
    for w in S4_SAMPLE:
        f = X(w) * X([2, 1, 3])
        g = f.expand()
        xs = [v for v in g.parent().gens() if str(v).startswith("x")]
        for i in range(1, len(xs)):
            swapped = g.subs({xs[i - 1]: xs[i], xs[i]: xs[i - 1]})
            expected = (g - swapped) // (xs[i - 1] - xs[i])
            lhs, rhs = common(f.divided_difference(i).expand(), expected)
            assert lhs == rhs, (w, i)


def test_invalid_permutation_rejected(X):
    with pytest.raises(ValueError):
        X([1, 2, 1])


# --- quantum ---------------------------------------------------------------------------------------


def test_quantum_product_matches_polynomial_product(Q):
    for u, v in itertools.product(S3, S3):
        assert product_matches(Q(u) * Q(v), Q(u), Q(v)), (u, v)


def test_quantum_reduces_to_classical_at_q_zero(Q, S):
    for u, v in itertools.product(S3, S3):
        lhs, a, b = common(zero_out((Q(u) * Q(v)).expand(), "q"), S(u).expand(), S(v).expand())
        assert lhs == a * b, (u, v)


def test_quantum_monk(Q):
    """S^q_{s_1} S^q_{s_1} = S^q_{312} + q_1 (Fomin-Gelfand-Postnikov)."""
    q0 = Q.base_ring().gen(0)[0]
    assert Q([2, 1]) * Q([2, 1]) == Q([3, 1, 2]) + q0 * Q.one()


def test_quantum_double_product_matches_polynomial_product(QD):
    for u, v in itertools.product(S3, S3):
        assert product_matches(QD(u) * QD(v), QD(u), QD(v)), (u, v)


def test_quantum_double_specializations(QD, Q, X):
    """y = 0 gives quantum Schubert; q = 0 gives double Schubert."""
    for w in S3 + S4_SAMPLE:
        p = QD(w).expand()
        lhs, rhs = common(zero_out(p, "y"), Q(w).expand())
        assert lhs == rhs, w
        lhs, rhs = common(zero_out(p, "q"), X(w).expand())
        assert lhs == rhs, w


# --- parabolic -------------------------------------------------------------------------------------


def test_parabolic_rejects_non_parabolic_permutations():
    P = QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3))
    with pytest.raises(ValueError):
        P([2, 1, 3])  # descent at 1 is inside the first block
    assert P([1, 3, 2, 4, 6, 5]) is not None  # descent at N_2 = 5, increasing beyond


def test_grassmannian_basis_is_schur_inside_the_box():
    """Blocks (2, 3) = Gr(2, 5): sigma_lambda for lambda in the 2 x 3 box is s_lambda(x0, x1), no q."""
    G = QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3))
    s = SymmetricFunctions(QQ).schur()
    for lam, w in [((1,), [1, 3, 2]), ((2, 1), [2, 4, 1, 3]), ((3, 2), [3, 5, 1, 2, 4]), ((3, 3), [4, 5, 1, 2, 3])]:
        schur = s[lam].expand(2)  # in x0, x1
        lhs, rhs = common(G(w).expand(), schur)
        assert lhs == rhs, lam
        assert not any(str(g).startswith("q") for g in G(w).expand().variables()), lam


def test_grassmannian_quantum_pieri_gr_2_5():
    """sigma_32 sigma_1 = sigma_33 + q sigma_1 in QH^*(Gr(2,5)), plus sigma_42 which lives past the box."""
    G = QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3))
    q0 = G.base_ring().gen(0)[0]
    s32, s1, s33, s42 = G([3, 5, 1, 2, 4]), G([1, 3, 2]), G([4, 5, 1, 2, 3]), G([3, 6, 1, 2, 4, 5])
    assert s32 * s1 == s33 + q0 * s1 + s42
    # the q's in sigma_42 cancel against q sigma_1 on expansion
    assert product_matches(s32 * s1, s32, s1)


@pytest.mark.parametrize("blocks", [(2, 3), (1, 2), (2, 2)])
def test_parabolic_product_matches_polynomial_product(blocks):
    P = QuantumSchubertPolynomialRing(QQ, parabolic=blocks)
    bounds = set(itertools.accumulate(blocks))
    perms = [w for n in (3, 4, 5) for w in Permutations(n) if set(w.descents()) <= bounds and w.length() <= 3]
    for u, v in itertools.product(perms, perms):
        assert product_matches(P(u) * P(v), P(u), P(v)), (blocks, u, v)


def test_parabolic_result_independent_of_implicit_block_size():
    """Extending (2, 3) by an explicit extra block of any size gives the same structure constants."""
    base = QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3))
    u, v = [3, 5, 1, 2, 4], [2, 4, 1, 3]
    reference = {tuple(w): str(c) for w, c in base(u) * base(v)}
    for extra in (1, 2, 3):
        P = QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3, extra))
        assert {tuple(w): str(c) for w, c in P(u) * P(v)} == reference, extra


def test_parabolic_quantum_double_product_matches_polynomial_product():
    PD = QuantumDoubleSchubertPolynomialRing(QQ, parabolic=(1, 2))
    perms = [w for w in Permutations(3) if set(w.descents()) <= {1, 3}]
    for u, v in itertools.product(perms, perms):
        assert product_matches(PD(u) * PD(v), PD(u), PD(v)), (u, v)
