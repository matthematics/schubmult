# from schubmult import *
# from schubmult.rings.polynomial_algebra import *
# from schubmult.rings.free_algebra import *
# from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
# from sympy import pretty_print



# if __name__ == "__main__":
#     import sys
#     import itertools

#     n = int(sys.argv[1])
#     perms = Permutation.all_permutations(n)
    
#     ForestPoly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
#     r = ForestRCGraphRing()
#     #CompSchub = PolynomialAlgebra(CompositionSchubertBasis)
#     length = n - 1
#     for perm1, perm2 in itertools.product(perms, repeat=2):
#         comp1 = perm1.pad_code(length)
#         comp2 = perm2.pad_code(length)
#         fschub1 = ForestPoly(*comp1)
#         fschub2 = ForestPoly(*comp2)
#         fprd = fschub1 * fschub2

#         cschub1 = Schub(perm1, length)
#         cschub2 = Schub(perm2, length)
#         cprd = cschub1 * cschub2
#         rc_reduce = r.zero
#         for (perm, lenn), v in cprd.items():
#             rc_reduce = 
#         nonsense = ForestPoly.from_dict({p.pad_code(lenn): v for (p, lenn), v in cprd.items() })
#         assert cprd.almosteq(nonsense), f"Failed for {perm1} and {perm2}, got {nonsense} but expected {fprd}"
        