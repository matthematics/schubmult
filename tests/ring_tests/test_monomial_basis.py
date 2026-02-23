def test_from_expr():
    from schubmult.rings.polynomial_algebra import PolynomialAlgebra, MonomialBasis
    from schubmult.symbolic import expand
    from schubmult.abc import x
    monom = PolynomialAlgebra(MonomialBasis(x))
    test_poly = (x[1]**2 + 2*x[1]*x[2] + x[2]**2)**2
    assert expand(monom.from_expr(test_poly).expand() - expand(test_poly)) == 0

def test_product():
    from schubmult.rings.polynomial_algebra import PolynomialAlgebra, MonomialBasis
    from schubmult.symbolic import expand
    from schubmult.abc import x
    monom = PolynomialAlgebra(MonomialBasis(x))
    poly1 = x[1]**2 + 2*x[1]*x[2] + x[2]**2
    poly2 = x[1] + x[2]
    result = poly1 * poly2
    monom1 = monom.from_expr(poly1)
    monom2 = monom.from_expr(poly2)
    result_monom = monom.from_expr(result)
    assert (monom1 * monom2).almosteq(result_monom)
    assert ((monom1 * monom2).expand() - (result_monom).expand()).expand() == 0
