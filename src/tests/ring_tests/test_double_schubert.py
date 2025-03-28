


def test_schub_expand():
    """
    Test expand
    """
    from schubmult.rings import DSx
    from sympy import symbols
    x_1, y_1, y_2, y_3 = symbols("x_1 y_1 y_2 y_3")
    # DSx = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", "y")
    assert DSx([3, 1, 2]).expand().equals((x_1 - y_1)*(x_1 - y_2))
    assert (DSx([3, 4, 1, 2]).expand() * DSx([4, 1, 2, 3]).expand()).equals((DSx([3, 4, 1, 2]) * DSx([4, 1, 2, 3])).expand())
    assert (x_1* DSx([3, 4, 1, 2])).expand().equals(y_3 * DSx([3, 4, 1, 2]).expand() + DSx([4, 3, 1, 2]).expand())


# def test_coproduct():
#     """
#     Test coproduct
#     """
#     from sage.all import ZZ

#     from schubmult.sage_integration import (
#         FastDoubleSchubertPolynomialRing,
#     )
#     DSx = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", "y")
#     R = DSx._base_polynomial_ring
#     indices = [0, 1, 3]
#     indices2 = [0, 2, 4]
#     DSx.set_coproduct_indices(indices[1:])
#     subs_dict1 = {R.gens()[i]: R.gens()[indices[i]] for i in range(len(indices))}
#     subs_dict2 = {R.gens()[i]: R.gens()[indices2[i]] for i in range(len(indices))}

#     perm = [3, 5, 1, 4, 2]
#     assert DSx(perm) == DSx(
#         sum(
#             [
#                 DSx._base_polynomial_ring(v)
#                 * (DSx(*k[0]).expand().subs(subs_dict1))
#                 * (DSx(*k[1]).expand().subs(subs_dict2))
#                 for k, v in DSx(perm).coproduct().monomial_coefficients().items()
#             ],
#         ),
#     )

#     # double coproduct
#     DSx = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
#     R = DSx._base_polynomial_ring
#     DSx.set_coproduct_indices(indices[1:])
#     subs_dict1 = {R.gens()[i]: R.gens()[indices[i]] for i in range(len(indices))}
#     subs_dict2 = {R.gens()[i]: R.gens()[indices2[i]] for i in range(len(indices))}
#     assert DSx(perm) == DSx(
#         sum(
#             [
#                 DSx._base_polynomial_ring(v)
#                 * (DSx(*k[0]).expand().subs(subs_dict1))
#                 * (DSx(*k[1]).expand().subs(subs_dict2))
#                 for k, v in DSx(perm).coproduct().monomial_coefficients().items()
#             ],
#         ),
#     )


def test_associative():
    """
    Test associative on some large perms
    """
    from schubmult.rings import DSx
    perm1 = [1, 5, 3, 4, 2]
    perm2 = [1, 3, 2, 5, 4]
    perm3 = [3, 1, 4, 2, 5]
    assert ((DSx(perm1) * DSx(perm2)) * DSx(perm3)).equals(DSx(perm1) * (DSx(perm2) * DSx(perm3)))

    assert ((DSx(perm1) * DSx(perm2)) * DSx(perm3, "z")).equals(DSx(perm1) * (DSx(perm2) * DSx(perm3, "z")))

def test_change_vars():
    """
    Test associative on some large perms
    """
    from schubmult.rings import DSx
    perm = [1, 5, 3, 4, 2]
    
    assert DSx(perm).change_vars("theta").equals(DSx(perm))


# def test_coerce():
#     from sage.all import ZZ

#     from schubmult.sage_integration import (
#         FastDoubleSchubertPolynomialRing,
#         FastQuantumDoubleSchubertPolynomialRing,
#         FastQuantumSchubertPolynomialRing,
#         FastSchubertPolynomialRing,
#     )
#     DSxD = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
#     DSxDQ = FastQuantumDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
#     DSxS = FastSchubertPolynomialRing(ZZ, 100, "x")
#     DSxSQ = FastQuantumSchubertPolynomialRing(ZZ, 100, "x")
#     R = DSxD._base_polynomial_ring

#     assert (
#         DSxD([2, 3, 5, 4, 1], "z") * R("x2^2") - R("x2^2") * DSxD([2, 3, 5, 4, 1], "z") == 0
#     )

#     assert DSxS([3,1,4,2]) * DSxD([4,1,3,2], "z") == DSxD([4,1,3,2], "z") * DSxS([3,1,4,2])
#     assert DSxS([3,1,4,2]) * DSxDQ([1], "z") == DSxSQ([1]) * DSxS([3,1,4,2])

# def test_mixed_equal():
#     """
#     Test mixed equality
#     """
#     from sage.all import ZZ

#     from schubmult.sage_integration import (
#         FastDoubleSchubertPolynomialRing,
#     )
#     DSx = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
#     R = DSx._coeff_polynomial_ring
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

#     assert sum([R(coeffs[i]) * DSx(perms[i]) for i in range(len(coeffs))]) == DSx(
#         other_perm, "z",
#     )

#     assert DSx([3, 1, 5, 2, 4]) * DSx([5, 3, 1, 2, 4], "z") == DSx([5, 3, 1, 2, 4], "z") * DSx(
#         [3, 1, 5, 2, 4],
#     )

#     assert DSx([3, 1, 5, 2, 4], "z") * DSx([5, 3, 1, 2, 4], "z") != DSx(
#         [5, 3, 1, 2, 4], "z",
#     ) * DSx([3, 1, 5, 2, 4])

#     assert DSx([3, 1, 5, 2, 4], "z") * DSx([5, 3, 1, 2, 4]) == DSx([5, 3, 1, 2, 4]) * DSx(
#         [3, 1, 5, 2, 4], "z",
#     )
