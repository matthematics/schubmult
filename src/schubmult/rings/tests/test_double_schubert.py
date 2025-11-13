


def test_schub_expand():
    """
    Test expand
    """
    from schubmult import DSx
    from schubmult.symbolic import symbols, expand
    x_1, y_1, y_2, y_3 = symbols("x_1 y_1 y_2 y_3")
    # DSx = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", "y")
    assert expand(DSx([3, 1, 2]).expand()) == expand((x_1 - y_1)*(x_1 - y_2))
    print(f"{expand(DSx([3, 4, 1, 2]).expand() * DSx([4, 1, 2, 3]).expand() - ((DSx([3, 4, 1, 2]) * DSx([4, 1, 2, 3])).expand()))=}")
    assert expand(DSx([3, 4, 1, 2]).expand() * DSx([4, 1, 2, 3]).expand() - ((DSx([3, 4, 1, 2]) * DSx([4, 1, 2, 3])).expand())) == 0
    assert expand((x_1* DSx([3, 4, 1, 2])).expand()) == expand(y_3 * DSx([3, 4, 1, 2]).expand() + DSx([4, 3, 1, 2]).expand())


def test_associative():
    """
    Test associative on some large perms
    """
    from schubmult import DSx
    from schubmult.symbolic import expand, S
    perm1, perm2, perm3 = [1, 5, 3, 4, 2], [1, 3, 2, 5, 4], [3, 1, 4, 2, 5]
    assert ((DSx(perm1) * DSx(perm2)) * DSx(perm3) - DSx(perm1) * (DSx(perm2) * DSx(perm3))).almosteq(S.Zero)

    assert (((DSx(perm1) * DSx(perm2)) * DSx(perm3, "z")) - (DSx(perm1) * (DSx(perm2) * DSx(perm3, "z")))).almosteq(S.Zero)

# def test_change_vars():
#     """
#     Test associative on some large perms
#     """
#     from schubmult import DSx
#     from schubmult.symbolic import expand
#     perm = [1, 5, 3, 4, 2]
    
#     assert expand(DSx([], "theta").basis.express(DSx(perm)) - DSx(perm).change_vars("gamama")) == 0
#     #assert (DSx(perm)*DSx(perm,"z")).change_vars("theta") == ((DSx(perm)*DSx(perm,"z")).change_vars("yourmom"))

def test_expr_trans():
    from schubmult import DSx
    from schubmult.symbolic import expand, S
    from schubmult.abc import x, y

    expr = (x[1] + y[1]*x[2])**5
    schub = DSx([]).ring.from_expr(expr)
    assert expand(schub.as_polynomial() - expr) == S.Zero

def test_elem_sym():
    from schubmult import DSx, Permutation
    from schubmult.symbolic import expand, S
    perm1 = Permutation([3,4,1,5,6,2])
    perm2 = Permutation([4,3,6,5,2,1])
    dct1 = DSx(perm1) * DSx(perm2, "z")
    dct2 = DSx(perm1,elem_sym=True) * DSx(perm2, "z")
    assert all(expand(dct1[k] - dct2.get(k, S.Zero), func=True) == S.Zero for k in dct1.keys())

def test_coproduct():
    from schubmult import DSx
    from schubmult.symbolic import expand
    perm = [3, 1, 5, 2, 6, 7, 4]
    indices = [2, 4, 6]
    assert expand(DSx(perm).coproduct(*indices).expand() - DSx(perm).expand()) == 0

def test_subs():
    from schubmult.symbolic import expand, S
    from schubmult import DSx, GeneratingSet
    x = GeneratingSet("x")
    z = GeneratingSet("z")
    perm = [4, 6, 1, 2, 3, 5]
    old = x[2]
    new = x[1]
    assert expand(DSx(perm,"z").subs(old, new) - DSx(perm,"z").as_polynomial().subs(old, new)) == S.Zero
    old = z[1]
    new = 3
    assert expand(DSx(perm,"z").subs(old, new) - DSx(perm,"z").as_polynomial().subs(old, new)) == S.Zero

# def test_elem_sym_schub():
#     from schubmult import DSx
#     from schubmult.abc import x, y, E, z
#     from schubmult.symbolic import expand, S
#     ring = DSx([]).ring
#     em = ring(E(2,3,[x[1],x[5],x[7]],[y[1],y[5],y[7],y[9]])*E(1,3,[x[1],x[4],x[5]],[y[8],y[14],y[29]]) - E(3,3,[x[1],x[7],x[3]],[z[4]]))
#     #em = E(2,3,x[1:],y[1:]) - E(1,3,x[1:],y[1:])
#     assert expand(em - ring.from_sympy(em),func=True) == S.Zero

# def test_complete_sym_schub():
#     from schubmult import DSx
#     from schubmult.abc import x, y, H, z
#     from schubmult.symbolic import expand, S
#     ring = DSx([]).ring
#     em = ring(H(2,3,[x[1],x[5],x[7]],[y[1],y[5],y[7],y[9],y[13],y[15]])*H(1,3,[x[1],x[4],x[5]],[y[8],y[14],y[29]]) - H(3,3,[x[1],x[7],x[3]],[z[4],z[15],y[1],y[12],z[2]]))
#     assert expand(em - ring.from_sympy(em),func=True) == S.Zero