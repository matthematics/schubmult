from schubmult import *
from schubmult.rings.polynomial_algebra import *
from sympy import pretty_print
def trimit(elem):
    #return sum([coeff * r(rc) @ Sx(perm) for (rc, perm), coeff in elem.items() if rc.perm == perm])
    return elem

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    r = DualRCGraphRing()
    perms = Permutation.all_permutations(n)
    k = n 
    bprod = 1
    for k in range(1, n):
        new_bprod = 0
        for p in range(k + 1):
            for rc in r.schub((uncode([0] * (k - p) + [1] * p)), n - 1):
                new_bprod += trimit(bprod * (r(rc) @ PA.from_expr(Sx.genset[n - k] ** (k - p), length=n-1)))
        bprod = new_bprod
    cprod = bprod
    # for (rc, tup), coeff in bprod.items():
    #     cprod += coeff * r(rc) @ FA(*tup).change_basis(SchubertBasis)
    # if any(v>1 for v in cprod.values()):
    #     raise ValueError("Coefficients should be 0 or 1")
    pretty_print(cprod)
                