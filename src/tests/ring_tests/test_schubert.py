


def test_schub_expand():
    """
    Test expand
    """
    from schubmult.rings import Sx
    from sympy import symbols, sympify
    x_1, x_2 = symbols("x_1 x_2")
    assert Sx([3,1,2]).expand() == x_1**2
    assert (Sx([5,3,4,1,2]).expand() * Sx([4,1,5,2,3]).expand() - Sx([5,3,4,1,2]) * Sx([4,1,5,2,3])).expand() == sympify(0)
    assert (x_1*Sx([3,4,1,2])).expand() == x_1**3*x_2**2

# def test_coproduct():
#     """
#     Test coproduct
#     """
#     from sage.all import ZZ

#     from schubmult.sage_integration import FastSchubertPolynomialRing
#     Sx = FastSchubertPolynomialRing(ZZ, 100, "x")
#     R = Sx._polynomial_ring
#     indices = [0, 1, 3, 5]
#     indices2 = [0, 2, 4, 6]
#     Sx.set_coproduct_indices(indices[1:])
#     subs_dict1 = {R.gens()[i]: R.gens()[indices[i]] for i in range(len(indices))}
#     subs_dict2 = {R.gens()[i]: R.gens()[indices2[i]] for i in range(len(indices))}

#     perm = [3,6,5,1,4,7,2]
#     assert Sx(perm) == Sx(sum([v*Sx(list(k[0])).expand().subs(subs_dict1)*(Sx(list(k[1])).expand().subs(subs_dict2)) for k,v in Sx(perm).coproduct().monomial_coefficients().items()]))

def test_associative():
    """
    Test associative on some large perms
    """
    from schubmult.rings import Sx
    perm1, perm2, perm3 = [6,7,1,5,3,4,2], [1,3,2,7,6,5,4], [3,7,1,6,4,2,5]
    assert ((Sx(perm1)*Sx(perm2))*Sx(perm3)) == ((Sx(perm1)*(Sx(perm2)*Sx(perm3))))
