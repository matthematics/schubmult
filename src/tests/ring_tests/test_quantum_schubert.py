def test_schub_expand():
    """
    Test expand
    """
    from schubmult.rings import QSx
    from sympy import symbols, expand
    from symengine import sympify
    x_1, x_2, q_1, q_2 = symbols("x_1 x_2 q_1 q_2")
    assert QSx([3,1,2]).expand() == x_1**2 - q_1
    assert expand((QSx([5,3,4,1,2]).expand() * QSx([4,1,5,2,3]).expand() - QSx([5,3,4,1,2]) * QSx([4,1,5,2,3])).expand()) == 0
    assert (x_1*QSx([3,4,1,2])).expand() == q_1**2*x_1 + q_1*q_2*x_1 + 2*q_1*x_1**2*x_2 - q_2*x_1**3 + x_1**3*x_2**2

# def test_coproduct():
#     from schubmult.rings import QSx
#     from sympy import expand
#     perm = [3, 1, 5, 2, 6, 7, 4]
#     indices = [2, 4, 6]
#     assert expand(QSx(perm).coproduct(indices).expand() - QSx(perm).expand()) == 0

def test_parabolic():
    from schubmult.rings import make_parabolic_quantum_basis
    from schubmult.perm_lib.perm_lib import uncode
    from symengine import S, expand
    QPDSx = make_parabolic_quantum_basis([2, 3])
    A = QPDSx(uncode([1,2]), 0)
    B = QPDSx(uncode([1,3]), 0)
    C = A*B
    assert expand(C.as_polynomial()-A.as_polynomial()*B.as_polynomial()) == S.Zero


def test_associative():
    """
    Test associative on some large perms
    """
    from schubmult.rings import QSx
    from symengine import expand
    perm1, perm2, perm3 = [1,5,3,4,2], [1,3,2,5,4], [3,1, 4,2,5]
    assert expand(((QSx(perm1)*QSx(perm2))*QSx(perm3)).as_polynomial() - ((QSx(perm1)*(QSx(perm2)*QSx(perm3)))).as_polynomial()) == 0

if __name__ == "__main__":
    test_schub_expand()
