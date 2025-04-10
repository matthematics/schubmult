


def test_schub_expand():
    """
    Test expand
    """
    from schubmult.rings import QDSx
    from symengine import symbols, expand
    x_1, y_1, y_2, y_3, q_1 = symbols("x_1 y_1 y_2 y_3 q_1")
    # QDSx = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", "y")
    assert expand(QDSx([3, 1, 2]).expand()) == -q_1 + x_1**2 - x_1*y_1 - x_1*y_2 + y_1*y_2
    print(f"{expand(QDSx([3, 4, 1, 2]).expand() * QDSx([4, 1, 2, 3]).expand() - ((QDSx([3, 4, 1, 2]) * QDSx([4, 1, 2, 3])).expand()))=}")
    assert expand(QDSx([3, 4, 1, 2]).expand() * QDSx([4, 1, 2, 3]).expand() - ((QDSx([3, 4, 1, 2]) * QDSx([4, 1, 2, 3])).expand())) == 0
    assert expand((x_1* QDSx([3, 4, 1, 2])).expand()) == expand(y_3 * QDSx([3, 4, 1, 2]).expand() + QDSx([4, 3, 1, 2]).expand())


def test_subs():
    from sympy import expand
    from schubmult.rings import QDSx
    from schubmult import GeneratingSet
    x = GeneratingSet("x")
    z = GeneratingSet("z")
    perm = [4, 6, 1, 2, 3, 5]
    old = x[5]
    new = x[1]
    assert expand(QDSx(perm,"z").subs(old, new) - QDSx(perm,"z").as_polynomial().subs(old, new)) == 0
    old = z[1]
    new = 3
    assert expand(QDSx(perm,"z").subs(old, new) - QDSx(perm,"z").as_polynomial().subs(old, new)) == 0

# def test_coproduct():
#     """
#     Test coproduct
#     """
#     from sympy import expand

#     from schubmult.rings import QDSx
#     indices = [0, 1, 3]
#     indices2 = [0, 2, 4]
    
#     perm = [3, 5, 1, 4, 2]
#     assert expand(QDSx(perm).expand() - QDSx(
#         sum(
#             [
#                 QDSx._base_polynomial_ring(v)
#                 * (QDSx(*k[0]).expand().subs(subs_dict1))
#                 * (QDSx(*k[1]).expand().subs(subs_dict2))
#                 for k, v in QDSx(perm).coproduct().monomial_coefficients().items()
#             ],
#         ),
#     )

#     # double coproduct
#     QDSx = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
#     R = QDSx._base_polynomial_ring
#     QDSx.set_coproduct_indices(indices[1:])
#     subs_dict1 = {R.gens()[i]: R.gens()[indices[i]] for i in range(len(indices))}
#     subs_dict2 = {R.gens()[i]: R.gens()[indices2[i]] for i in range(len(indices))}
#     assert QDSx(perm) == QDSx(
#         sum(
#             [
#                 QDSx._base_polynomial_ring(v)
#                 * (QDSx(*k[0]).expand().subs(subs_dict1))
#                 * (QDSx(*k[1]).expand().subs(subs_dict2))
#                 for k, v in QDSx(perm).coproduct().monomial_coefficients().items()
#             ],
#         ),
#     )


def test_associative():
    """
    Test associative on some large perms
    """
    from schubmult.rings import QDSx
    from sympy import expand
    perm1, perm2, perm3 = [1, 5, 3, 4, 2], [1, 3, 2, 5, 4], [3, 1, 4, 2, 5]
    assert expand(((QDSx(perm1) * QDSx(perm2)) * QDSx(perm3)) - (QDSx(perm1) * (QDSx(perm2) * QDSx(perm3)))) == 0

    assert expand(((QDSx(perm1) * QDSx(perm2)) * QDSx(perm3, "z")) - (QDSx(perm1) * (QDSx(perm2) * QDSx(perm3, "z")))) == 0

# def test_change_vars():
#     """
#     Test associative on some large perms
#     """
#     from schubmult.rings import QDSx
#     from sympy import expand
#     perm = [1, 5, 3, 4, 2]
    
#     assert expand(QDSx(perm).change_vars("theta") - QDSx(perm).change_vars("gamama")) == 0
#     #assert (QDSx(perm)*QDSx(perm,"z")).change_vars("theta") == ((QDSx(perm)*QDSx(perm,"z")).change_vars("yourmom"))


# def test_coproduct():
#     from schubmult.rings import QDSx
#     from sympy import expand
#     perm = [3, 1, 5, 2, 6, 7, 4]
#     indices = [2, 4, 6]
#     assert expand(QDSx(perm).coproduct(indices).expand() - QDSx(perm).expand()) == 0

# def test_coerce():
#     from sage.all import ZZ

#     from schubmult.sage_integration import (
#         FastDoubleSchubertPolynomialRing,
#         FastQuantumDoubleSchubertPolynomialRing,
#         FastQuantumSchubertPolynomialRing,
#         FastSchubertPolynomialRing,
#     )
#     QDSxD = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
#     QDSxDQ = FastQuantumDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
#     QDSxS = FastSchubertPolynomialRing(ZZ, 100, "x")
#     QDSxSQ = FastQuantumSchubertPolynomialRing(ZZ, 100, "x")
#     R = QDSxD._base_polynomial_ring

#     assert (
#         QDSxD([2, 3, 5, 4, 1], "z") * R("x2^2") - R("x2^2") * QDSxD([2, 3, 5, 4, 1], "z") == 0
#     )

#     assert QDSxS([3,1,4,2]) * QDSxD([4,1,3,2], "z") == QDSxD([4,1,3,2], "z") * QDSxS([3,1,4,2])
#     assert QDSxS([3,1,4,2]) * QDSxDQ([1], "z") == QDSxSQ([1]) * QDSxS([3,1,4,2])

# def test_mixed_equal():
#     """
#     Test mixed equality
#     """
#     from sage.all import ZZ

#     from schubmult.sage_integration import (
#         FastDoubleSchubertPolynomialRing,
#     )
#     QDSx = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
#     R = QDSx._coeff_polynomial_ring
#     perms = [
#         [1],
#         [1, 2, 4, 3],
#         [2, 1],
#         [1, 2, 5, 3, 4],
#         [2, 1, 4, 3],
#         [2, 1, 5, 3, 4],
#         [3, 1, 4, 2],
#         [3, 1, 5, 2, 4],
#     ]
#     coeffs = [
#         "(y1 - z3)*((y1 - z1)*(y1 - z2)*(y_3 - z2) + (y_2 - z1)*(y1 - z1)*(y1 - z2)) + (y_3 - z4)*((y1 - z1)*(y1 - z2)*(y_3 - z2) + (y_2 - z1)*(y1 - z1)*(y1 - z2)) + (y_2 - z1)*(y1 - z1)*(y1 - z2)*(y_2 - z2)",
#         "(y1 - z1)*(y1 - z2)*(y_3 - z2) + (y1 - z1)*(y1 - z2)*(y4 - z4) + (y1 - z1)*(y1 - z3)*(y1 - z2) + (y_2 - z1)*(y1 - z1)*(y1 - z2)",
#         "(y_2 - z3)*((y1 - z1)*(y_3 - z2) + (y_2 - z1)*(y1 - z1) + (y_3 - z2)*(y_2 - z2)) + (y_3 - z4)*((y1 - z1)*(y_3 - z2) + (y_2 - z1)*(y1 - z1) + (y_3 - z2)*(y_2 - z2)) + (y1 - z1)*(y1 - z2)*(y_3 - z2) + (y_2 - z1)*(y1 - z1)*(y1 - z2)",
#         "(y1 - z1)*(y1 - z2)",
#         "(y1 - z1)*(y1 - z2) + (y1 - z1)*(y_3 - z2) + (y_2 - z1)*(y1 - z1) + (y_3 - z2)*(y_2 - z2) + (y1 + y_2 - z1 - z2)*(y_2 - z3) + (y1 + y_2 - z1 - z2)*(y4 - z4)",
#         "y1 + y_2 - z1 - z2",
#         "y1 + y_2 + y_3 + y4 - z1 - z2 - z3 - z4",
#         1,
#     ]

#     other_perm = [3, 1, 5, 2, 4]

#     assert len(perms) == len(coeffs)

#     assert sum([R(coeffs[i]) * QDSx(perms[i]) for i in range(len(coeffs))]) == QDSx(
#         other_perm, "z",
#     )

#     assert QDSx([3, 1, 5, 2, 4]) * QDSx([5, 3, 1, 2, 4], "z") == QDSx([5, 3, 1, 2, 4], "z") * QDSx(
#         [3, 1, 5, 2, 4],
#     )

#     assert QDSx([3, 1, 5, 2, 4], "z") * QDSx([5, 3, 1, 2, 4], "z") != QDSx(
#         [5, 3, 1, 2, 4], "z",
#     ) * QDSx([3, 1, 5, 2, 4])

#     assert QDSx([3, 1, 5, 2, 4], "z") * QDSx([5, 3, 1, 2, 4]) == QDSx([5, 3, 1, 2, 4]) * QDSx(
#         [3, 1, 5, 2, 4], "z",
#     )
