"""Duality between the free algebra and the polynomial algebra.

The free algebra is the graded dual of the polynomial algebra, so each free
algebra basis is the dual basis of the correspondingly named polynomial basis.
These tests pin that relationship down:

* ``dual_basis()`` is an involution between the two families,
* the pairing of dual basis elements is the Kronecker delta,
* the pairing is bilinear, so pairing against a basis element reads off the
  coefficient of the dual basis element,
* ``free.pairing(poly)`` and ``poly.apply_dual_element(free)`` agree.
"""

import pytest

# Small set of compositions used to index basis elements. Keys of degree <= 2
# are enough to exercise the delta property without slow expansions.
COMPOSITIONS = [(), (1,), (0, 1), (2,), (1, 1)]

# Elementary keys ``(tup, numvars)``: ``tup[:numvars-1]`` is the flag part, ``tup[numvars-1:]`` the
# sorted symmetric tail. Includes same-degree pairs whose tails collide as polynomials
# (e.g. e_2(x1,x2) vs e_1(x1,x2)^2) since those are where a per-key staircase fails the delta.
ELEM_SYM_COMPS = [((), 0), ((1, 0), 2), ((0, 1, 1), 2), ((0, 2), 2), ((1, 1), 2), ((1, 0, 2), 3), ((1, 0, 1, 1), 3), ((0, 1, 1, 1), 3), ((0, 0, 3), 3), ((0, 0, 1, 2), 3)]

COEFFICIENTS = [2, -3, 5, 7, -11]

# Free algebra basis name -> True if the basis is indexed by permutations,
# False if it is indexed by compositions.
DUAL_PAIRS = {
    "ForestBasis": False,
    "FundamentalSlideBasis": False,
    "GlideBasis": False,
    "GrothendieckBasis": True,
    "GroveBasis": False,
    "KeyBasis": False,
    "LascouxBasis": False,
    "MonomialSlideBasis": False,
    "SchubertBasis": True,
    "WordBasis": False,
    "ElementaryBasis": None,
}

BASIS_NAMES = sorted(DUAL_PAIRS.keys())


def _dual_pair(name):
    """Return ``(free_algebra, polynomial_algebra, keys)`` for a declared dual pair."""
    import schubmult.rings.free_algebra as free_algebra
    from schubmult import uncode
    from schubmult.abc import x
    from schubmult.rings.free_algebra import FreeAlgebra
    from schubmult.rings.polynomial_algebra import PolynomialAlgebra

    basis = getattr(free_algebra, name)
    dual = basis.dual_basis()
    if DUAL_PAIRS[name] is None:
        keys = list(ELEM_SYM_COMPS)
    elif DUAL_PAIRS[name]:
        keys = [(uncode(list(comp)), len(comp)) for comp in COMPOSITIONS]
    else:
        keys = list(COMPOSITIONS)
    return FreeAlgebra(basis=basis), PolynomialAlgebra(basis=dual(x)), keys


@pytest.mark.parametrize("name", BASIS_NAMES)
def test_free_basis_dual_is_an_involution(name):
    """``B.dual_basis().dual_basis()`` returns ``B`` for every free algebra basis."""
    import schubmult.rings.free_algebra as free_algebra

    basis = getattr(free_algebra, name)
    dual = basis.dual_basis()
    assert dual is not None, f"{name} declares no dual basis"
    assert dual.dual_basis() is basis, f"{name} -> {dual.__name__} -> {dual.dual_basis()}"


@pytest.mark.parametrize("name", BASIS_NAMES)
def test_polynomial_basis_dual_is_an_involution(name):
    """The polynomial side of each pair points back at the free algebra basis."""
    import schubmult.rings.free_algebra as free_algebra

    basis = getattr(free_algebra, name)
    dual = basis.dual_basis()
    back = dual.dual_basis()
    assert back is not None, f"{dual.__name__} declares no dual basis"
    assert back.dual_basis() is dual


def test_composition_schubert_basis_shares_the_schubert_dual():
    """``CompositionSchubertBasis`` is a reindexing of ``SchubertBasis`` and delegates its dual."""
    from schubmult.rings.free_algebra import CompositionSchubertBasis, SchubertBasis
    from schubmult.rings.polynomial_algebra import SchubertPolyBasis

    assert CompositionSchubertBasis.dual_basis() is SchubertPolyBasis
    assert SchubertPolyBasis.dual_basis() is SchubertBasis


@pytest.mark.parametrize("name", BASIS_NAMES)
def test_pairing_of_dual_bases_is_kronecker_delta(name):
    """<B(u), D(v)> == 1 if u == v else 0, which is what makes the bases dual."""
    algebra, poly, keys = _dual_pair(name)
    for i, left in enumerate(keys):
        for j, right in enumerate(keys):
            expected = 1 if i == j else 0
            assert algebra(left).pairing(poly(right)) == expected, f"{name}: <{left}, {right}> should be {expected}"


# @pytest.mark.parametrize("name", BASIS_NAMES)
# def test_pairing_agrees_with_apply_dual_element(name):
#     """The pairing may be evaluated from either side."""
#     algebra, poly, keys = _dual_pair(name)
#     for left in keys:
#         for right in keys:
#             free_elem = algebra(left)
#             poly_elem = poly(right)
#             assert free_elem.pairing(poly_elem) == free_elem.apply_dual_element(poly_elem), f"{name}: pairing and apply_dual_element disagree for {left}, {right}"
#             poly_elem = poly(right)
#             assert free_elem.pairing(poly_elem) == poly_elem.apply_dual_element(free_elem)


@pytest.mark.parametrize("name", BASIS_NAMES)
def test_pairing_extracts_polynomial_coefficients(name):
    """Pairing a basis element against a polynomial reads off its coefficient."""
    algebra, poly, keys = _dual_pair(name)
    element = sum(coeff * poly(key) for coeff, key in zip(COEFFICIENTS, keys))
    for coeff, key in zip(COEFFICIENTS, keys):
        assert algebra(key).pairing(element) == coeff, f"{name}: coefficient of {key}"


@pytest.mark.parametrize("name", BASIS_NAMES)
def test_pairing_extracts_free_algebra_coefficients(name):
    """Pairing a free algebra element against a dual basis element reads off its coefficient."""
    algebra, poly, keys = _dual_pair(name)
    element = sum(coeff * algebra(key) for coeff, key in zip(COEFFICIENTS, keys))
    for coeff, key in zip(COEFFICIENTS, keys):
        assert element.pairing(poly(key)) == coeff, f"{name}: coefficient of {key}"


@pytest.mark.parametrize("name", BASIS_NAMES)
def test_pairing_is_bilinear(name):
    """The pairing is linear in each argument, so it is determined by the delta property."""
    algebra, poly, keys = _dual_pair(name)
    free_elem = sum(coeff * algebra(key) for coeff, key in zip(COEFFICIENTS, keys))
    poly_elem = sum(coeff * poly(key) for coeff, key in zip(COEFFICIENTS, keys))
    assert free_elem.pairing(poly_elem) == sum(coeff * coeff for coeff in COEFFICIENTS)
    assert free_elem.pairing(poly_elem) == poly_elem.apply_dual_element(free_elem)


def _dual_basis_pairing(free_elem, poly_elem):
    """Pair elements (or tensors) written in dual bases: matching keys contribute the product of coefficients."""
    return sum(v * poly_elem.get(k, 0) for k, v in free_elem.items())


@pytest.mark.parametrize("name", BASIS_NAMES)
def test_product_is_adjoint_to_dual_coproduct(name):
    """<a * b, p> == <a (x) b, coproduct(p)> with a, b in B and p in B.dual_basis()."""
    algebra, poly, keys = _dual_pair(name)
    
    for key in keys:
        poly_elem_cprd = poly(key).coproduct()

        for (pkey1, pkey2), coeff in poly_elem_cprd.items():

            freeprod = algebra(pkey1) * algebra(pkey2)

            assert freeprod.get(key, 0) == coeff, f"{name}: product and coproduct mismatch for {key}, {pkey1}, {pkey2}"


@pytest.mark.parametrize("name", BASIS_NAMES)
def test_coproduct_is_adjoint_to_dual_product(name):
    """<coproduct(f), p (x) q> == <f, p * q> with f in B and p, q in B.dual_basis()."""
    algebra, poly, keys = _dual_pair(name)
    
    for key in keys:
        free_elem_cprd = algebra(key).coproduct()

        for (fkey1, fkey2), coeff in free_elem_cprd.items():

            polyprod = (poly(fkey1) * poly(fkey2))

            assert polyprod.get(key, 0) == coeff, f"{name}: product and coproduct mismatch for {key}, {fkey1}, {fkey2}"
