


def test_schub_expand():
    """
    Test expand
    """
    from schubmult import Sx
    from schubmult.symbolic import symbols, expand
    from schubmult.symbolic import sympify
    x_1, x_2 = symbols("x_1 x_2")
    assert Sx([3,1,2]).expand() == x_1**2
    assert expand((Sx([5,3,4,1,2]).expand() * Sx([4,1,5,2,3]).expand() - (Sx([5,3,4,1,2]) * Sx([4,1,5,2,3])).expand()).expand()) == 0
    assert (x_1*Sx([3,4,1,2])).expand() == x_1**3*x_2**2

def test_coproduct():
    from schubmult import Sx
    from schubmult.symbolic import expand
    perm = [3, 1, 5, 2, 6, 7, 4]
    indices = [2, 4, 6]
    assert expand(Sx(perm).coproduct(*indices).expand() - Sx(perm).expand()) == 0

def test_expr_trans():
    from schubmult import Sx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1]*x[2])**5
    schub = Sx([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero


def test_associative():
    """
    Test associative on some large perms
    """
    from schubmult import Sx
    perm1, perm2, perm3 = [6,1,5,3,4,2], [1,3,2,6,5,4], [3,1,6,4,2,5]
    assert (((Sx(perm1)*Sx(perm2))*Sx(perm3)) - ((Sx(perm1)*(Sx(perm2)*Sx(perm3))))).almosteq(0)

if __name__ == "__main__":
    test_schub_expand()
