
def test_schub_expand():
    """
    Test expand
    """
    from schubmult import QDSx
    from schubmult.symbolic import symbols, expand
    x_1, y_1, y_2, y_3, q_1 = symbols("x_1 y_1 y_2 y_3 q_1")
    # QDSx = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", "y")
    assert expand(QDSx([3, 1, 2]).expand()) == -q_1 + x_1**2 - x_1*y_1 - x_1*y_2 + y_1*y_2
    print(f"{expand(QDSx([3, 4, 1, 2]).expand() * QDSx([4, 1, 2, 3]).expand() - ((QDSx([3, 4, 1, 2]) * QDSx([4, 1, 2, 3])).expand()))=}")
    assert expand(QDSx([3, 4, 1, 2]).expand() * QDSx([4, 1, 2, 3]).expand() - ((QDSx([3, 4, 1, 2]) * QDSx([4, 1, 2, 3])).expand())) == 0
    assert expand((x_1* QDSx([3, 4, 1, 2])).expand()) == expand(y_3 * QDSx([3, 4, 1, 2]).expand() + QDSx([4, 3, 1, 2]).expand())

def test_parabolic():
    from schubmult import QPDSx, uncode
    from schubmult.symbolic import S, expand
    #QPDSx = make_parabolic_quantum_basis([2, 3])
    A = QPDSx(2,3)(uncode([1,2]))
    B = QPDSx(2,3)(uncode([1,3]), "z")
    C = A*B
    assert expand(C.as_polynomial()-A.as_polynomial()*B.as_polynomial()) == S.Zero
    assert expand(C.as_polynomial()-QPDSx(2,3)(C.as_polynomial()).as_polynomial()) == S.Zero

def test_expr_trans():
    from schubmult import QDSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1]*x[2])**5
    schub = QDSx([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero

def test_expr_trans_parabolic():
    from schubmult import QPDSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1]*(x[2]+x[3])**2)**2
    schub = QPDSx(1,2)([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero


def test_associative():
    """
    Test associative on some large perms
    """
    from schubmult import QDSx
    from schubmult.symbolic import expand, S
    perm1, perm2, perm3 = [1, 5, 3, 4, 2], [1, 3, 2, 5, 4], [3, 1, 4, 2, 5]
    assert (((QDSx(perm1) * QDSx(perm2)) * QDSx(perm3)) - (QDSx(perm1) * (QDSx(perm2) * QDSx(perm3)))).almosteq(S.Zero)

    assert (((QDSx(perm1) * QDSx(perm2)) * QDSx(perm3, "z")) - (QDSx(perm1) * (QDSx(perm2) * QDSx(perm3, "z")))).almosteq(S.Zero)

#

if __name__ == "__main__":
    test_subs()