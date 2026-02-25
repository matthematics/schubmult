def test_schubert_schub_expand():
    from schubmult import Sx
    from schubmult.symbolic import symbols, expand

    x_1, x_2 = symbols("x_1 x_2")
    assert Sx([3, 1, 2]).expand() == x_1 ** 2
    assert expand((Sx([5, 3, 4, 1, 2]).expand() * Sx([4, 1, 5, 2, 3]).expand() - (Sx([5, 3, 4, 1, 2]) * Sx([4, 1, 5, 2, 3])).expand()).expand()) == 0
    assert (x_1 * Sx([3, 4, 1, 2])).expand() == x_1 ** 3 * x_2 ** 2


def test_schubert_coproduct():
    from schubmult import Sx
    from schubmult.symbolic import expand

    perm = [3, 1, 5, 2, 6, 7, 4]
    indices = [2, 4, 6]
    assert expand(Sx(perm).coproduct(*indices).expand() - Sx(perm).expand()) == 0


def test_schubert_expr_trans():
    from schubmult import Sx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1] * x[2]) ** 5
    schub = Sx([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero


def test_schubert_associative():
    from schubmult import Sx

    perm1, perm2, perm3 = [6, 1, 5, 3, 4, 2], [1, 3, 2, 6, 5, 4], [3, 1, 6, 4, 2, 5]
    assert (((Sx(perm1) * Sx(perm2)) * Sx(perm3)) - (Sx(perm1) * (Sx(perm2) * Sx(perm3)))).almosteq(0)


def test_double_schubert_expand():
    from schubmult import DSx
    from schubmult.symbolic import symbols, expand

    x_1, y_1, y_2, y_3 = symbols("x_1 y_1 y_2 y_3")
    assert expand(DSx([3, 1, 2]).expand()) == expand((x_1 - y_1) * (x_1 - y_2))
    assert expand(DSx([3, 4, 1, 2]).expand() * DSx([4, 1, 2, 3]).expand() - ((DSx([3, 4, 1, 2]) * DSx([4, 1, 2, 3])).expand())) == 0
    assert expand((x_1 * DSx([3, 4, 1, 2])).expand()) == expand(y_3 * DSx([3, 4, 1, 2]).expand() + DSx([4, 3, 1, 2]).expand())


def test_double_schubert_associative():
    from schubmult import DSx
    from schubmult.symbolic import S

    perm1, perm2, perm3 = [1, 5, 3, 4, 2], [1, 3, 2, 5, 4], [3, 1, 4, 2, 5]
    assert ((DSx(perm1) * DSx(perm2)) * DSx(perm3) - DSx(perm1) * (DSx(perm2) * DSx(perm3))).almosteq(S.Zero)
    assert (((DSx(perm1) * DSx(perm2)) * DSx(perm3, "z")) - (DSx(perm1) * (DSx(perm2) * DSx(perm3, "z")))).almosteq(S.Zero)


def test_double_schubert_expr_trans():
    from schubmult import DSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1] * x[2]) ** 5
    schub = DSx([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero


def test_double_schubert_elem_sym():
    from schubmult import DSx, Permutation
    from schubmult.symbolic import expand, S

    perm1 = Permutation([3, 4, 1, 5, 6, 2])
    perm2 = Permutation([4, 3, 6, 5, 2, 1])
    dct1 = DSx(perm1) * DSx(perm2, "z")
    dct2 = DSx(perm1, elem_sym=True) * DSx(perm2, "z")
    assert all(expand(dct1[k] - dct2.get(k, S.Zero), func=True) == S.Zero for k in dct1.keys())


def test_double_schubert_coproduct():
    from schubmult import DSx
    from schubmult.symbolic import expand

    perm = [3, 1, 5, 2, 6, 7, 4]
    indices = [2, 4, 6]
    assert expand(DSx(perm).coproduct(*indices).expand() - DSx(perm).expand()) == 0


def test_double_schubert_subs():
    from schubmult.symbolic import expand, S
    from schubmult import DSx, GeneratingSet

    x = GeneratingSet("x")
    z = GeneratingSet("z")
    perm = [4, 6, 1, 2, 3, 5]
    old = x[2]
    new = x[1]
    assert expand(DSx(perm, "z").subs(old, new) - DSx(perm, "z").as_polynomial().subs(old, new)) == S.Zero
    old = z[1]
    new = 3
    assert expand(DSx(perm, "z").subs(old, new) - DSx(perm, "z").as_polynomial().subs(old, new)) == S.Zero


def test_quantum_schub_expand():
    from schubmult import QSx
    from schubmult.symbolic import symbols, expand

    x_1, x_2, q_1, q_2 = symbols("x_1 x_2 q_1 q_2")
    assert QSx([3, 1, 2]).expand() == x_1 ** 2 - q_1
    assert expand((QSx([5, 3, 4, 1, 2]).expand() * QSx([4, 1, 5, 2, 3]).expand() - (QSx([5, 3, 4, 1, 2]) * QSx([4, 1, 5, 2, 3])).expand()).expand()) == 0
    assert (x_1 * QSx([3, 4, 1, 2])).expand() == q_1 ** 2 * x_1 + q_1 * q_2 * x_1 + 2 * q_1 * x_1 ** 2 * x_2 - q_2 * x_1 ** 3 + x_1 ** 3 * x_2 ** 2


def test_quantum_schubert_expr_trans():
    from schubmult import QSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1] * x[2]) ** 5
    schub = QSx([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero


def test_quantum_schubert_expr_trans_parabolic():
    from schubmult import QPSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1] * (x[2] + x[3]) ** 2) ** 2
    schub = QPSx(1, 2)([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero


def test_quantum_schubert_parabolic():
    from schubmult import QPSx
    from schubmult import uncode
    from schubmult.symbolic import S, expand

    a = QPSx(2, 3, 4)(uncode([1, 2, 0, 2, 3]))
    b = QPSx(2, 3, 4)(uncode([1, 3, 0, 1, 2]))
    c = a * b
    assert expand(c.as_polynomial() - a.as_polynomial() * b.as_polynomial()) == S.Zero


def test_quantum_schubert_associative():
    from schubmult import QSx

    perm1, perm2, perm3 = [1, 5, 3, 4, 2], [1, 3, 2, 5, 4], [3, 1, 4, 2, 5]
    assert (((QSx(perm1) * QSx(perm2)) * QSx(perm3)) - ((QSx(perm1) * (QSx(perm2) * QSx(perm3))))).almosteq(0)


def test_quantum_double_schub_expand():
    from schubmult import QDSx
    from schubmult.symbolic import symbols, expand

    x_1, y_1, y_2, y_3, q_1 = symbols("x_1 y_1 y_2 y_3 q_1")
    assert expand(QDSx([3, 1, 2]).expand()) == -q_1 + x_1 ** 2 - x_1 * y_1 - x_1 * y_2 + y_1 * y_2
    assert expand(QDSx([3, 4, 1, 2]).expand() * QDSx([4, 1, 2, 3]).expand() - ((QDSx([3, 4, 1, 2]) * QDSx([4, 1, 2, 3])).expand())) == 0
    assert expand((x_1 * QDSx([3, 4, 1, 2])).expand()) == expand(y_3 * QDSx([3, 4, 1, 2]).expand() + QDSx([4, 3, 1, 2]).expand())


def test_quantum_double_schubert_parabolic():
    from schubmult import QPDSx, uncode
    from schubmult.symbolic import S, expand

    a = QPDSx(2, 3)(uncode([1, 2]))
    b = QPDSx(2, 3)(uncode([1, 3]), "z")
    c = a * b
    assert expand(c.as_polynomial() - a.as_polynomial() * b.as_polynomial()) == S.Zero
    assert expand(c.as_polynomial() - QPDSx(2, 3)(c.as_polynomial()).as_polynomial()) == S.Zero


def test_quantum_double_schubert_expr_trans():
    from schubmult import QDSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1] * x[2]) ** 5
    schub = QDSx([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero


def test_quantum_double_schubert_expr_trans_parabolic():
    from schubmult import QPDSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1] * (x[2] + x[3]) ** 2) ** 2
    schub = QPDSx(1, 2)([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero


def test_quantum_double_schubert_associative():
    from schubmult import QDSx
    from schubmult.symbolic import S

    perm1, perm2, perm3 = [1, 5, 3, 4, 2], [1, 3, 2, 5, 4], [3, 1, 4, 2, 5]
    assert (((QDSx(perm1) * QDSx(perm2)) * QDSx(perm3)) - (QDSx(perm1) * (QDSx(perm2) * QDSx(perm3)))).almosteq(S.Zero)
    assert (((QDSx(perm1) * QDSx(perm2)) * QDSx(perm3, "z")) - (QDSx(perm1) * (QDSx(perm2) * QDSx(perm3, "z")))).almosteq(S.Zero)


def test_elem_sym():
    from schubmult.abc import x, y, E
    from schubmult.symbolic import expand, S

    elem = E(2, 5, [x[1], x[4], x[3], x[5], x[7]], [y[3], y[1], y[4], y[2]])
    assert expand(elem - elem.split_out_vars([x[3], x[5]]), func=True) == S.Zero
    assert elem.xreplace({y[2]: x[4]}) == E(2, 4, [x[1], x[3], x[5], x[7]], [y[3], y[1], y[4]])


def test_complete_sym():
    from schubmult.abc import x, y, H
    from schubmult.symbolic import expand, S

    elem = H(5, 4, x[1], x[4], x[3], x[5], y[3], y[1], y[4], y[2], y[7], y[9], y[11], y[15])
    assert expand(elem - elem.split_out_vars([x[3], x[1]]), func=True) == S.Zero
    assert elem.xreplace({y[2]: x[4]}) == H(5, 3, [x[1], x[3], x[5]], [y[3], y[1], y[4], y[7], y[9], y[11], y[15]])
