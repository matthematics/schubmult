"""Property tests for ``schubmult.sage.PolynomialAlgebra`` (the polynomial ring with its combinatorial bases).

Skipped when Sage is not importable. Every basis is checked against polynomial arithmetic after
``expand()`` (done by Sage), and the key, Schubert and Grothendieck bases against Sage's own
implementations resp. the Grothendieck ring.
"""

import itertools
import re

import pytest

sage_all = pytest.importorskip("sage.all")

from sage.all import QQ, ZZ, KeyPolynomials, Permutation, Permutations, PolynomialRing, SchubertPolynomialRing  # noqa: E402
from sage.misc.sage_unittest import TestSuite as SageTestSuite  # noqa: E402

from schubmult.sage import DoubleSchubertPolynomialRing, GrothendieckPolynomialRing, PolynomialAlgebra  # noqa: E402

COMPOSITIONS = [(2, 0, 1), (1, 2), (0, 1, 1), (0, 0, 2), (1, 0, 0, 1)]
PERMS = [Permutation(p) for p in ([2, 1], [1, 3, 2], [3, 1, 2], [2, 3, 1], [1, 2, 4, 3], [2, 1, 4, 3])]
COMPOSITION_BASES = ["monomial", "key", "fundamental_slide", "monomial_slide", "forest", "glide", "lascoux", "grove"]
PERMUTATION_BASES = ["schubert", "grothendieck"]
T = PolynomialRing(QQ, ["x0", "x1", "x2", "x3", "x4"])


@pytest.fixture(scope="module")
def A():
    return PolynomialAlgebra(QQ)


@pytest.fixture(scope="module")
def B3():
    return PolynomialAlgebra(QQ, 3)


def as_T(p):
    return T(str(p).replace("z_", "x"))  # Sage's key polynomials live in z_i


def keys_for(name):
    return PERMS if name in PERMUTATION_BASES else COMPOSITIONS


def fits(B, key):
    try:
        B._key(key)
    except ValueError:
        return False
    return True


# --- Sage category framework -----------------------------------------------------------------------


@pytest.mark.parametrize("name", COMPOSITION_BASES + PERMUTATION_BASES)
def test_sage_testsuite(A, name):
    SageTestSuite(A.basis(name)).run(raise_on_failure=True)


@pytest.mark.parametrize("name", [*COMPOSITION_BASES, *PERMUTATION_BASES, "elementary"])
def test_sage_testsuite_finite(B3, name):
    SageTestSuite(B3.basis(name)).run(raise_on_failure=True)


def test_algebra_testsuites():
    SageTestSuite(PolynomialAlgebra(ZZ)).run(raise_on_failure=True)
    SageTestSuite(PolynomialAlgebra(ZZ, 2)).run(raise_on_failure=True)
    assert PolynomialAlgebra(QQ, 3) is PolynomialAlgebra(QQ, ZZ(3))
    with pytest.raises(ValueError):
        PolynomialAlgebra(QQ).elementary()
    with pytest.raises(ValueError):
        PolynomialAlgebra(QQ, 0)


def test_coefficients_in_a_polynomial_base_ring():
    """Coefficients that are not rational numbers survive the round trip through schubmult (regression:
    the conversion back assumed every coefficient was an integer or a rational)."""
    R = PolynomialRing(QQ, "a")
    a = R.gen()
    A = PolynomialAlgebra(R)
    x, k, S, G = A.monomial(), A.key(), A.schubert(), A.grothendieck()
    c = a**2 - QQ((1, 2))
    f = a * x[1] + c * x[0, 2]
    assert k(f) == a * k[1] + c * (k[0, 2] - k[2] - k[1, 1])  # x1^2 = kappa_02 - kappa_20 - kappa_11
    assert x(k(f)) == f
    assert S(G(f)) == S(f)
    assert (a * S[2, 1]) * (a * S[1, 3, 2]) == a**2 * (S[2, 3, 1] + S[3, 1, 2])
    p = f.expand()
    assert p.parent().base_ring() is R
    assert p == a * p.parent()("x0") + c * p.parent()("x1^2")
    assert x(p) == f
    SageTestSuite(k).run(raise_on_failure=True)
    with pytest.raises(ValueError):
        PolynomialAlgebra(QQ).key()._scalar("a")  # a symbol the base ring does not have


# --- products and basis changes --------------------------------------------------------------------


@pytest.mark.parametrize("name", COMPOSITION_BASES + PERMUTATION_BASES)
def test_products_match_polynomial_products(A, name):
    B = A.basis(name)
    for a, b in itertools.product(keys_for(name), keys_for(name)):
        assert as_T((B(a) * B(b)).expand()) == as_T(B(a).expand()) * as_T(B(b).expand()), (name, a, b)


@pytest.mark.parametrize("target", COMPOSITION_BASES + PERMUTATION_BASES)
def test_basis_changes_preserve_the_polynomial(A, target):
    """Coercion source -> target (through the monomial basis) keeps the polynomial and round-trips."""
    Bt = A.basis(target)
    for source in COMPOSITION_BASES + PERMUTATION_BASES:
        Bs = A.basis(source)
        for a in keys_for(source):
            f = Bs(a)
            g = Bt(f)
            assert as_T(g.expand()) == as_T(f.expand()), (source, target, a)
            assert Bs(g) == f, (source, target, a)
        assert (Bs(keys_for(source)[0]) + Bt(keys_for(target)[1])).parent() in (Bs, Bt)


@pytest.mark.parametrize("name", COMPOSITION_BASES + PERMUTATION_BASES)
def test_finite_algebra_agrees_with_the_infinite_one(A, B3, name):
    """In 3 variables every basis element is the same polynomial as in the infinite ring, products agree,
    and everything round-trips through the elementary symmetric basis."""
    B, E = B3.basis(name), B3.elementary()
    keys = [k for k in keys_for(name) if fits(B, k)]
    for a in keys:
        f = B(a)
        assert as_T(f.expand()) == as_T(A.basis(name)(a).expand()), (name, a)
        assert as_T(E(f).expand()) == as_T(f.expand()), (name, a)
        assert B(E(f)) == f, (name, a)
    for a, b in itertools.product(keys, keys):
        assert as_T((B(a) * B(b)).expand()) == as_T(B(a).expand()) * as_T(B(b).expand()), (name, a, b)


def test_elementary_basis(B3):
    E = B3.elementary()
    keys = [(1, 0, 2), (0, 1, 1), (0, 0, 1, 1), (1, 2, 0), (0, 0, 3), (1, 0, 0)]
    for a, b in itertools.product(keys, keys):
        assert as_T((E(a) * E(b)).expand()) == as_T(E(a).expand()) * as_T(E(b).expand()), (a, b)
    assert E([1, 0, 2, 0, 1]) == E[1, 0, 1, 2]  # e_0 factors drop out, full-alphabet degrees sorted
    assert E.one() == E[0, 0, 0]
    assert E[0, 0, 1].expand() == T("x0 + x1 + x2")
    assert E[1, 1, 0].expand() == T("x0^2 + x0*x1")
    assert E[0, 0, 1, 1].expand() == T("(x0 + x1 + x2)^2")
    with pytest.raises(ValueError):
        E([2, 0, 0])  # e_2 of a single variable
    with pytest.raises(ValueError):
        E([0, 0, 4])  # e_4 of three variables


def test_finite_algebra_rejects_too_many_variables(B3):
    with pytest.raises(ValueError):
        B3.key()([1, 0, 0, 2])
    with pytest.raises(ValueError):
        B3.schubert()([1, 2, 3, 5, 4])
    with pytest.raises(ValueError):
        B3.monomial()(T("x3"))
    assert B3.key()([2, 0, 1, 0, 0]) == B3.key()[2, 0, 1]
    assert (B3.key()[1, 0, 0] + PolynomialAlgebra(QQ, 2).key()[0, 1]).parent() is B3.key()  # fewer variables embed


# --- indexing --------------------------------------------------------------------------------------


def test_composition_keys_are_trimmed_and_graded(A):
    k = A.key()
    assert k([2, 0, 1, 0, 0]) == k[2, 0, 1]
    assert k([]) == k.one() == k(1)
    assert k[2, 0, 1].degree() == 3
    assert A.lascoux()[2, 0, 1].expand().degree() > 3  # inhomogeneous (beta = -1)
    with pytest.raises(TypeError):
        k(Permutation([3, 1, 2]))


def test_permutation_keys(A):
    S = A.schubert()
    assert S([3, 1, 2, 4, 5]) == S[3, 1, 2] == S(Permutation([3, 1, 2]))
    assert S.one() == S(Permutation([1]))
    assert S[3, 1, 2].degree() == 2
    assert S[2, 1, 4, 3].expand() == T("x0^2 + x0*x1 + x0*x2")


# --- against Sage and the other schubmult rings ----------------------------------------------------


def test_key_basis_matches_sage(A):
    k, K = A.key(), KeyPolynomials(QQ)
    for alpha in COMPOSITIONS + [(3, 0, 1, 1), (0, 2, 1)]:
        assert as_T(k(alpha).expand()) == as_T(K(list(alpha)).expand()), alpha
        assert k(K(list(alpha))) == k(alpha), alpha
    assert (k[1] + K([0, 1])).parent() is k


def test_schubert_basis_matches_sage(A):
    S, X = A.schubert(), SchubertPolynomialRing(QQ)
    for w in PERMS + list(Permutations(4)):
        assert as_T(S(w).expand()) == as_T(X(w).expand()), w
        assert S(X(w)) == S(w), w
    assert X(S[3, 1, 2].expand()) == X([3, 1, 2])
    assert (S[2, 1] + X([1, 3, 2])).parent() is S


def test_grothendieck_basis_matches_grothendieck_ring(A):
    """The Grothendieck basis is the Grothendieck ring at beta = -1 (the classical convention)."""
    G, S = A.grothendieck(), A.schubert()
    GR = GrothendieckPolynomialRing(QQ)
    for w in PERMS:
        p = GR(w).expand()
        assert as_T(G(w).expand()) == as_T(p.subs({p.parent()("beta"): -1})), w
    for u, v in itertools.product(PERMS[:4], PERMS[:4]):
        f = GR(u) * GR(v)
        expected = {tuple(w): c(-1) for w, c in f}
        assert {tuple(w): c for w, c in G(u) * G(v)} == expected, (u, v)
    assert G[1, 3, 2].expand() == T("x0 + x1 - x0*x1")
    assert S(G[1, 3, 2]) == S[1, 3, 2] - S[2, 3, 1]


def test_k_theoretic_bases_are_at_beta_minus_one(A):
    """Independent of the sign twist: the grove basis agrees with schubmult's ``GrovePolyBasis(beta=-1)``,
    and every K-theoretic basis element is its beta = 1 version with x -> -x and the sign (-1)^|alpha|."""
    from schubmult.rings.polynomial_algebra import GlidePolyBasis, GrovePolyBasis, LascouxPolyBasis
    from schubmult.rings.polynomial_algebra import PolynomialAlgebra as SPA
    from schubmult.symbolic.poly.variables import GeneratingSet

    x = GeneratingSet("x")
    grove_minus = SPA(GrovePolyBasis(x, beta=-1))
    negate = {T.gen(i): -T.gen(i) for i in range(5)}

    def sage_poly(expr):  # schubmult's 1-indexed x_i -> Sage's 0-indexed xi
        return T(re.sub(r"x_(\d+)", lambda m: f"x{int(m.group(1)) - 1}", str(expr)).replace("**", "^"))

    for a in COMPOSITIONS:
        assert as_T(A.grove()(a).expand()) == sage_poly(grove_minus.from_dict({a: 1}).expand()), a
        for name, basis in (("glide", GlidePolyBasis), ("lascoux", LascouxPolyBasis), ("grove", GrovePolyBasis)):
            plus = sage_poly(SPA(basis(x)).from_dict({a: 1}).expand())
            assert as_T(A.basis(name)(a).expand()) == (-1) ** sum(a) * plus.subs(negate), (name, a)


def test_k_theoretic_bases_specialize_to_their_classical_ones(A):
    """The lowest-degree part of a glide / Lascoux / grove polynomial is the slide / key / forest one."""
    pairs = [("glide", "fundamental_slide"), ("lascoux", "key"), ("grove", "forest")]
    for kname, cname in pairs:
        for a in COMPOSITIONS:
            p = as_T(A.basis(kname)(a).expand())
            low = sum((c * m for c, m in zip(p.coefficients(), p.monomials()) if m.degree() == sum(a)), T.zero())
            assert low == as_T(A.basis(cname)(a).expand()), (kname, a)


def test_coerces_into_schubert_family_rings(A):
    X = DoubleSchubertPolynomialRing(QQ)
    k = A.key()
    f = X(k[0, 2])
    assert as_T(f.expand().subs({g: 0 for g in f.expand().parent().gens() if str(g).startswith("y")})) == as_T(k[0, 2].expand())
    assert (X([2, 1]) + k[1]).parent() is X
    assert X(A.schubert()[3, 1, 2]) == X(SchubertPolynomialRing(QQ)([3, 1, 2]))
