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
    "SchubertBasis": True,
    "WordBasis": False,
}

BASIS_NAMES = sorted(DUAL_PAIRS)


def _dual_pair(name):
    """Return ``(free_algebra, polynomial_algebra, keys)`` for a declared dual pair."""
    import schubmult.rings.free_algebra as free_algebra
    from schubmult import uncode
    from schubmult.abc import x
    from schubmult.rings.free_algebra import FreeAlgebra
    from schubmult.rings.polynomial_algebra import PolynomialAlgebra

    basis = getattr(free_algebra, name)
    dual = basis.dual_basis()
    if DUAL_PAIRS[name]:
        keys = [uncode(list(comp)) for comp in COMPOSITIONS]
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


@pytest.mark.parametrize("name", BASIS_NAMES)
def test_pairing_agrees_with_apply_dual_element(name):
    """The pairing may be evaluated from either side."""
    algebra, poly, keys = _dual_pair(name)
    for left in keys:
        for right in keys:
            free_elem = algebra(left)
            poly_elem = poly(right)
            assert free_elem.pairing(poly_elem) == poly_elem.apply_dual_element(free_elem)


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
