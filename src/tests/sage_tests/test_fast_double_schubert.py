from sage.all import ZZ
from schubmult.sage_integration import (
    FastDoubleSchubertPolynomialRing,
    FastSchubertPolynomialRing,
    FastQuantumSchubertPolynomialRing,
    FastQuantumDoubleSchubertPolynomialRing,
)


def test_schub_expand():
    """
    Test expand
    """
    X = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", "y")
    assert X._base_polynomial_ring(X([3, 1, 2]).expand()) == X._base_polynomial_ring(
        "(x1 - y1)*(x1 - y2)"
    )
    assert X(X([3, 4, 1, 2]).expand() * X([4, 1, 2, 3]).expand()) == X(
        [3, 4, 1, 2]
    ) * X([4, 1, 2, 3])
    assert (
        X(X._base_polynomial_ring("x1")) * X([3, 4, 1, 2])
    ).expand() == X._base_polynomial_ring("y3") * (X([3, 4, 1, 2]).expand()) + X(
        [4, 3, 1, 2]
    ).expand()


def test_coproduct():
    """
    Test coproduct
    """
    X = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", "y")
    R = X._base_polynomial_ring
    indices = [0, 1, 3]
    indices2 = [0, 2, 4]
    X.set_coproduct_indices(indices[1:])
    subs_dict1 = {R.gens()[i]: R.gens()[indices[i]] for i in range(len(indices))}
    subs_dict2 = {R.gens()[i]: R.gens()[indices2[i]] for i in range(len(indices))}

    perm = [3, 5, 1, 4, 2]
    assert X(perm) == X(
        sum(
            [
                X._base_polynomial_ring(v)
                * (X(*k[0]).expand().subs(subs_dict1))
                * (X(*k[1]).expand().subs(subs_dict2))
                for k, v in X(perm).coproduct().monomial_coefficients().items()
            ]
        )
    )

    # double coproduct
    X = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
    R = X._base_polynomial_ring
    X.set_coproduct_indices(indices[1:])
    subs_dict1 = {R.gens()[i]: R.gens()[indices[i]] for i in range(len(indices))}
    subs_dict2 = {R.gens()[i]: R.gens()[indices2[i]] for i in range(len(indices))}
    assert X(perm) == X(
        sum(
            [
                X._base_polynomial_ring(v)
                * (X(*k[0]).expand().subs(subs_dict1))
                * (X(*k[1]).expand().subs(subs_dict2))
                for k, v in X(perm).coproduct().monomial_coefficients().items()
            ]
        )
    )


def test_associative():
    """
    Test associative on some large perms
    """
    X = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", "y")
    perm1 = [1, 5, 3, 4, 2]
    perm2 = [1, 3, 2, 5, 4]
    perm3 = [3, 1, 4, 2, 5]
    assert (X(perm1) * X(perm2)) * X(perm3) == X(perm1) * (X(perm2) * X(perm3))

    X = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
    assert (X(perm1) * X(perm2)) * X(perm3, "z") == X(perm1) * (
        X(perm2) * X(perm3, "z")
    )


def test_coerce():
    XD = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
    XDQ = FastQuantumDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
    XS = FastSchubertPolynomialRing(ZZ, 100, "x")
    XSQ = FastQuantumSchubertPolynomialRing(ZZ, 100, "x")
    R = XD._base_polynomial_ring

    assert (
        XD([2, 3, 5, 4, 1], "z") * R("x2^2") - R("x2^2") * XD([2, 3, 5, 4, 1], "z") == 0
    )

    assert XS([3,1,4,2]) * XD([4,1,3,2], "z") == XD([4,1,3,2], "z") * XS([3,1,4,2])
    assert XS([3,1,4,2]) * XDQ([1], "z") == XSQ([1]) * XS([3,1,4,2])

def test_mixed_equal():
    """
    Test mixed equality
    """
    X = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
    R = X._coeff_polynomial_ring
    perms = [
        [1],
        [1, 2, 4, 3],
        [2, 1],
        [1, 2, 5, 3, 4],
        [2, 1, 4, 3],
        [2, 1, 5, 3, 4],
        [3, 1, 4, 2],
        [3, 1, 5, 2, 4],
    ]
    coeffs = [
        "(y1 - z3)*((y1 - z1)*(y1 - z2)*(y3 - z2) + (y2 - z1)*(y1 - z1)*(y1 - z2)) + (y3 - z4)*((y1 - z1)*(y1 - z2)*(y3 - z2) + (y2 - z1)*(y1 - z1)*(y1 - z2)) + (y2 - z1)*(y1 - z1)*(y1 - z2)*(y2 - z2)",
        "(y1 - z1)*(y1 - z2)*(y3 - z2) + (y1 - z1)*(y1 - z2)*(y4 - z4) + (y1 - z1)*(y1 - z3)*(y1 - z2) + (y2 - z1)*(y1 - z1)*(y1 - z2)",
        "(y2 - z3)*((y1 - z1)*(y3 - z2) + (y2 - z1)*(y1 - z1) + (y3 - z2)*(y2 - z2)) + (y3 - z4)*((y1 - z1)*(y3 - z2) + (y2 - z1)*(y1 - z1) + (y3 - z2)*(y2 - z2)) + (y1 - z1)*(y1 - z2)*(y3 - z2) + (y2 - z1)*(y1 - z1)*(y1 - z2)",
        "(y1 - z1)*(y1 - z2)",
        "(y1 - z1)*(y1 - z2) + (y1 - z1)*(y3 - z2) + (y2 - z1)*(y1 - z1) + (y3 - z2)*(y2 - z2) + (y1 + y2 - z1 - z2)*(y2 - z3) + (y1 + y2 - z1 - z2)*(y4 - z4)",
        "y1 + y2 - z1 - z2",
        "y1 + y2 + y3 + y4 - z1 - z2 - z3 - z4",
        1,
    ]

    other_perm = [3, 1, 5, 2, 4]

    assert len(perms) == len(coeffs)

    assert sum([R(coeffs[i]) * X(perms[i]) for i in range(len(coeffs))]) == X(
        other_perm, "z"
    )

    assert X([3, 1, 5, 2, 4]) * X([5, 3, 1, 2, 4], "z") == X([5, 3, 1, 2, 4], "z") * X(
        [3, 1, 5, 2, 4]
    )

    assert X([3, 1, 5, 2, 4], "z") * X([5, 3, 1, 2, 4], "z") != X(
        [5, 3, 1, 2, 4], "z"
    ) * X([3, 1, 5, 2, 4])

    assert X([3, 1, 5, 2, 4], "z") * X([5, 3, 1, 2, 4]) == X([5, 3, 1, 2, 4]) * X(
        [3, 1, 5, 2, 4], "z"
    )
