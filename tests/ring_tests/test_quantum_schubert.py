def test_schub_expand():
    """
    Test expand
    """
    from schubmult import QSx
    from schubmult.symbolic import symbols, expand
    from schubmult.symbolic import sympify
    x_1, x_2, q_1, q_2 = symbols("x_1 x_2 q_1 q_2")
    assert QSx([3,1,2]).expand() == x_1**2 - q_1
    assert expand((QSx([5,3,4,1,2]).expand() * QSx([4,1,5,2,3]).expand() - (QSx([5,3,4,1,2]) * QSx([4,1,5,2,3])).expand()).expand()) == 0
    assert (x_1*QSx([3,4,1,2])).expand() == q_1**2*x_1 + q_1*q_2*x_1 + 2*q_1*x_1**2*x_2 - q_2*x_1**3 + x_1**3*x_2**2

def test_expr_trans():
    from schubmult import QSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1]*x[2])**5
    schub = QSx([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero

def test_expr_trans_parabolic():
    from schubmult import QPSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1]*(x[2]+x[3])**2)**2
    schub = QPSx(1,2)([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero

def test_parabolic():
    from schubmult import QPSx
    from schubmult import uncode
    from schubmult.symbolic import S, expand
    # QPDSx = make_parabolic_quantum_basis([2, 3, 4])
    A = QPSx(2,3,4)(uncode([1,2,0,2,3]))
    B = QPSx(2,3,4)(uncode([1,3,0,1,2]))
    C = A*B
    assert expand(C.as_polynomial()-A.as_polynomial()*B.as_polynomial()) == S.Zero
    # assert expand(C.as_polynomial()-QPDSx(C.as_polynomial(), 0).as_polynomial()) == S.Zero


def test_associative():
    """
    Test associative on some large perms
    """
    from schubmult import QSx
    from schubmult.symbolic import expand
    perm1, perm2, perm3 = [1,5,3,4,2], [1,3,2,5,4], [3,1, 4,2,5]
    assert (((QSx(perm1)*QSx(perm2))*QSx(perm3)) - ((QSx(perm1)*(QSx(perm2)*QSx(perm3))))).almosteq(0)

if __name__ == "__main__":
    test_schub_expand()
