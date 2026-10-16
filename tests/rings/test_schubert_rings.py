def test_schubert_schub_expand():
    from schubmult import Sx
    from schubmult.symbolic import symbols, expand

    x_1, x_2 = symbols("x_1 x_2")
    assert Sx([3, 1, 2]).expand() == x_1 ** 2
    assert expand((Sx([5, 3, 4, 1, 2]).expand() * Sx([4, 1, 5, 2, 3]).expand() - (Sx([5, 3, 4, 1, 2]) * Sx([4, 1, 5, 2, 3])).expand()).expand()) == 0
    assert (x_1 * Sx([3, 4, 1, 2])).expand() == x_1 ** 3 * x_2 ** 2


def test_strip_zeros_exact_finds_hidden_zeros():
    """Unexpanded coefficients that cancel are invisible to the structural strip_zeros(); exact=True
    removes them by random-point evaluation, agreeing with full coefficient expansion."""
    from schubmult import DSx, QDSx
    from schubmult.abc import y, z
    from schubmult.symbolic import S, sympify

    ring = DSx([3, 1, 2]).ring
    w = next(iter(DSx([3, 1, 2])))
    hidden = sympify((y[1] - z[2]) * (y[3] - z[1]) - (y[1] * y[3] - y[1] * z[1] - z[2] * y[3] + z[1] * z[2]))
    elem = ring.from_dict({w: hidden, next(iter(DSx([2, 1]))): y[1] - z[1]})
    assert len(elem.strip_zeros()) == 2
    assert dict(elem.strip_zeros(exact=True)) == {next(iter(DSx([2, 1]))): y[1] - z[1]}

    # a product with hidden zeros: structural stripping sees none, exact stripping matches expansion
    # (from_dict drops literal zeros) and never drops a nonzero term. The number of hidden zeros
    # depends on how the kernel groups its unexpanded terms, so only the nonzero count is pinned.
    prod = QDSx([2, 4, 1, 3]) * QDSx([3, 1, 4, 2])
    assert len(prod.strip_zeros()) == len(prod) > 10
    stripped = prod.strip_zeros(exact=True)
    assert len(stripped) == 10
    assert set(stripped) == set(prod.expand(deep=False))
    assert all(k in prod for k in stripped)


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


def test_grothendieck_expand():
    from schubmult import Gx
    from schubmult.symbolic import symbols, expand

    x_1, x_2 = symbols("x_1 x_2")
    beta = Gx._beta
    assert Gx([3, 1, 2]).expand() == x_1 ** 2
    assert expand(Gx([1, 3, 2]).expand()) == expand(x_1 + x_2 + beta * x_1 * x_2)
    assert expand((Gx([5, 3, 4, 1, 2]).expand() * Gx([4, 1, 5, 2, 3]).expand() - (Gx([5, 3, 4, 1, 2]) * Gx([4, 1, 5, 2, 3])).expand()).expand()) == 0
    assert (x_1 * Gx([3, 4, 1, 2])).expand() == x_1 ** 3 * x_2 ** 2


def test_grothendieck_expr_trans():
    from schubmult import Gx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1] * x[2]) ** 5
    groth = Gx([]).ring.from_expr(expr)
    assert expand(groth.as_polynomial() - expr) == S.Zero


def test_grothendieck_associative():
    from schubmult import Gx
    from schubmult.symbolic import expand

    perm1, perm2, perm3 = [1, 5, 3, 4, 2], [1, 3, 2, 5, 4], [3, 1, 4, 2, 5]
    diff = ((Gx(perm1) * Gx(perm2)) * Gx(perm3)) - (Gx(perm1) * (Gx(perm2) * Gx(perm3)))
    assert all(expand(v) == 0 for v in diff.values())


def test_grothendieck_mult_poly():
    from schubmult import Gx
    from schubmult.abc import x
    from schubmult.symbolic import expand, S

    elem = Gx([1, 3, 2])
    assert expand((x[1] * elem).expand() - expand(x[1] * elem.expand())) == S.Zero


def test_grothendieck_beta_zero_is_schubert():
    from schubmult import Gx, Sx, uncode
    from schubmult.symbolic import expand, sympify_sympy

    perm = uncode([1, 1])
    at_zero = sympify_sympy(Gx(perm).expand()).xreplace({sympify_sympy(Gx._beta): 0})
    assert expand(at_zero - sympify_sympy(Sx(perm).expand())) == 0


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


def test_double_grothendieck_expand():
    from schubmult import DGx
    from schubmult.symbolic import symbols, expand, S

    x_1, y_1, y_2 = symbols("x_1 y_1 y_2")
    beta = DGx([]).ring.beta
    # [3, 1, 2] is dominant with code (2, 0), so G is the product over its two cells.
    assert expand(DGx([3, 1, 2]).expand() - expand((x_1 + y_1 + beta * x_1 * y_1) * (x_1 + y_2 + beta * x_1 * y_2))) == S.Zero


def test_double_grothendieck_product_polynomial():
    from sympy import cancel
    from schubmult import DGx
    from schubmult.symbolic import expand, sympify_sympy

    for perm1, perm2 in (([2, 1], [1, 3, 2]), ([1, 3, 2], [1, 3, 2]), ([3, 1, 2], [2, 1])):
        resid = expand((DGx(perm1) * DGx(perm2)).expand() - expand(DGx(perm1).expand() * DGx(perm2).expand()))
        assert cancel(sympify_sympy(resid)) == 0


def test_double_grothendieck_mult_poly():
    from sympy import cancel
    from schubmult import DGx
    from schubmult.abc import x
    from schubmult.symbolic import expand, sympify_sympy

    elem = DGx([1, 3, 2])
    resid = expand((x[1] * elem).expand() - expand(x[1] * elem.expand()))
    assert cancel(sympify_sympy(resid)) == 0


def test_double_grothendieck_associative():
    from sympy import cancel
    from schubmult import DGx
    from schubmult.symbolic import sympify_sympy

    # The right bracketing multiplies by a multi-term non-Monk element, which is the
    # path that bypasses grothmult_double.
    perm1, perm2, perm3 = [2, 1], [1, 3, 2], [2, 1]
    diff = ((DGx(perm1) * DGx(perm2)) * DGx(perm3)) - (DGx(perm1) * (DGx(perm2) * DGx(perm3)))
    assert all(cancel(sympify_sympy(v)) == 0 for v in diff.values())


def test_double_grothendieck_second_genset():
    from sympy import cancel
    from schubmult import DGx
    from schubmult.symbolic import expand, sympify_sympy

    a, b = DGx([2, 1]), DGx([1, 3, 2], "z")
    resid = expand((a * b).as_polynomial() - a.as_polynomial() * b.as_polynomial())
    assert cancel(sympify_sympy(resid)) == 0


def test_double_grothendieck_expr_trans():
    from sympy import cancel
    from schubmult import DGx
    from schubmult.abc import x, y
    from schubmult.symbolic import expand, sympify_sympy

    expr = (x[1] + y[1] * x[2]) ** 2
    groth = DGx([]).ring.from_expr(expr)
    assert cancel(sympify_sympy(expand(groth.as_polynomial() - expr))) == 0


def test_double_grothendieck_poly_round_trip():
    from schubmult import DGx
    from schubmult.utils.test_utils import vanishes

    elem = DGx([2, 1]) * DGx([1, 3, 2])
    diff = elem.ring.from_expr(elem.as_polynomial()) - elem
    # coefficients are large unsimplified rational functions; sympy.cancel took ~90s here
    assert vanishes(diff.values())


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

    # keep the second factor small: [1, 3, 0, 1, 2] made this an 82-term, 12s product
    a = QPSx(2, 3, 4)(uncode([1, 2, 0, 2, 3]))
    b = QPSx(2, 3, 4)(uncode([1, 1, 0, 1, 1]))
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


# --- Skew divided differences (Macdonald; B.-H. Hwang, arXiv:2610.21053) --------------------------------------


def _partial(poly, perm, genset):
    """partial_w on a raw polynomial, peeling the last letter of a reduced word (rightmost acts first)."""
    from schubmult import Permutation
    from schubmult.symbolic.poly.schub_poly import div_diff

    perm = Permutation(perm)
    while perm.inv:
        d = max(perm.descents(zero_indexed=False))
        perm = perm.swap(d - 1, d)
        poly = div_diff(poly, genset[d], genset[d + 1])
    return poly


def test_skew_div_diff_generalized_leibniz_rule():
    """partial_w(fg) = sum_v (partial_{w/v} f)(partial_v g) (Macdonald II.2.17 in the Gao--Xiong normalization)."""
    from schubmult import Permutation
    from schubmult.abc import x
    from schubmult.symbolic import expand
    from schubmult.symbolic.poly.schub_poly import skew_div_diff

    f = x[1] ** 2 * x[2] + 3 * x[1] * x[3] - x[2] * x[3] ** 2
    g = x[2] ** 2 + x[1] * x[2] * x[3] + 2 * x[1]
    S3 = Permutation.all_permutations(3)
    for w in S3:
        lhs = _partial(expand(f * g), w, x)
        rhs = sum(skew_div_diff(f, w, v, x) * _partial(g, v, x) for v in S3)
        assert expand(lhs - rhs) == 0, w


def test_skew_div_diff_endpoints():
    """partial_{w/id} = partial_w, partial_{w/w} = w (variables x_i -> x_{w(i)}), and zero unless v <= w."""
    from schubmult import Permutation
    from schubmult.abc import x
    from schubmult.symbolic import expand
    from schubmult.symbolic.poly.schub_poly import skew_div_diff

    f = x[1] ** 2 * x[2] + 3 * x[1] * x[3]
    w = Permutation([3, 1, 2])
    assert expand(skew_div_diff(f, w, [], x) - _partial(f, w, x)) == 0
    assert expand(skew_div_diff(f, w, w, x) - f.xreplace({x[1]: x[3], x[2]: x[1], x[3]: x[2]})) == 0
    assert skew_div_diff(f, [2, 1, 3], [1, 3, 2], x) == 0
    assert skew_div_diff(f, [], [2, 1], x) == 0
    assert skew_div_diff(f, [], [], x) == f


def test_skew_divdiff_matches_polynomial_operator_and_paper_example():
    """DSx.skew_divdiff (via c^w_{u,v}(y, x), Prop. 2.7) agrees with the polynomial recurrence; running example of the paper."""
    from schubmult import DSx, Permutation, Sx
    from schubmult.abc import x, y
    from schubmult.symbolic import expand
    from schubmult.symbolic.poly.schub_poly import skew_div_diff

    S3 = Permutation.all_permutations(3)
    for u in S3:
        poly = DSx(u).expand()
        for w in S3:
            for v in S3:
                assert expand(DSx(u).skew_divdiff(w, v).expand() - skew_div_diff(poly, w, v, x)) == 0, (u, w, v)

    u, v, w = [4, 3, 5, 1, 2, 6], [1, 4, 2, 6, 3, 5], [1, 6, 4, 3, 2, 5]
    e1, e2, e3 = y[1] + y[2] + y[3], y[1] * y[2] + y[1] * y[3] + y[2] * y[3], y[1] * y[2] * y[3]
    expected = (x[1] ** 3 - e1 * x[1] ** 2 + e2 * x[1] - e3) * (x[4] + x[6] - y[1] - y[2])
    assert expand(DSx(u).skew_divdiff(w, v).expand() - expected) == 0
    assert expand(Sx(u).skew_divdiff(w, v).expand() - x[1] ** 3 * (x[4] + x[6])) == 0
    # d = 0: the skew divided difference is the structure constant itself
    assert DSx([2, 1, 3]).skew_divdiff([3, 1, 2], [1, 3, 2]).expand() == 1
    assert DSx([2, 1, 3]).skew_divdiff([1, 3, 2], [2, 1, 3]).expand() == 0
    # linearity over the coefficient ring
    elem = 2 * DSx([2, 3, 1]) + y[1] * DSx([3, 1, 2])
    assert expand(elem.skew_divdiff([3, 1, 2], [2, 1, 3]).expand() - 2 * DSx([2, 3, 1]).skew_divdiff([3, 1, 2], [2, 1, 3]).expand() - y[1] * DSx([3, 1, 2]).skew_divdiff([3, 1, 2], [2, 1, 3]).expand()) == 0


def test_skew_divdiff_single_coefficients_are_block_lift_structure_constants():
    """Corollary 5.1: partial_{w/v} S_u(x) = sum_{|alpha| = d} c^{w~}_{u~,v~} x^alpha, so it lies in N[x]."""
    import itertools

    import sympy

    from schubmult import Permutation, Sx
    from schubmult.abc import x
    from schubmult.mult.single import schubmult_py
    from schubmult.symbolic import sympify_sympy
    from schubmult.utils.perm_utils import block_lift

    n = 4
    gens = [sympify_sympy(x[i]) for i in range(1, n + 1)]
    S4 = Permutation.all_permutations(n)
    for u, v, w in [(S4[k], S4[(7 * k + 3) % 24], S4[(5 * k + 11) % 24]) for k in range(24)]:
        d = u.inv + v.inv - w.inv
        skew = Sx(u).skew_divdiff(w, v).expand()
        if d < 0 or not v.bruhat_leq(w):
            assert skew == 0
            continue
        poly = sympy.Poly(sympify_sympy(skew), *gens)
        for alpha in itertools.product(range(d + 1), repeat=n):
            if sum(alpha) != d:
                continue
            ul, vl, wl = block_lift(u, v, w, alpha)
            coeff = poly.coeff_monomial(sympy.Mul(*[g**a for g, a in zip(gens, alpha)]))
            assert coeff >= 0
            assert coeff == schubmult_py({ul: 1}, vl).get(wl, 0), (u, v, w, alpha)
