from sage.all import ZZ
from schubmult.sage_integration import FastSchubertPolynomialRing


def test_schub_expand():
    """
        Test expand
    """
    X = FastSchubertPolynomialRing(ZZ, 100, "x")
    R = X._polynomial_ring
    assert X([3,1,2]).expand() == R("x1^2")
    assert X([5,3,4,1,2]).expand() * X([4,1,5,2,3]).expand() - X([5,3,4,1,2]) * X([4,1,5,2,3]) == 0
    assert R("x1")*X([3,4,1,2]) == R("x1^3*x2^2")

def test_coproduct():
    """
        Test coproduct
    """
    X = FastSchubertPolynomialRing(ZZ, 100, "x")
    R = X._polynomial_ring
    indices = [0, 1, 3, 5]
    indices2 = [0, 2, 4, 6]
    X.set_coproduct_indices(indices[1:])
    subs_dict1 = {R.gens()[i]: R.gens()[indices[i]] for i in range(len(indices))}
    subs_dict2 = {R.gens()[i]: R.gens()[indices2[i]] for i in range(len(indices))}

    perm = [3,6,5,1,4,7,2]
    assert X(perm) == X(sum([v*X(list(k[0])).expand().subs(subs_dict1)*(X(list(k[1])).expand().subs(subs_dict2)) for k,v in X(perm).coproduct().monomial_coefficients().items()]))

def test_associative():
    """
    Test associative on some large perms
    """
    X = FastSchubertPolynomialRing(ZZ, 100, "x")
    perm1 = [6,7,1,5,3,4,2]
    perm2 = [1,3,2,7,6,5,4]
    perm3 = [3,7,1,6,4,2,5]
    assert (X(perm1)*X(perm2))*X(perm3) == X(perm1)*(X(perm2)*X(perm3))