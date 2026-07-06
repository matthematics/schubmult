# from schubmult import *
# from schubmult.symbolic.poly.variables import MaskedGeneratingSet
# from schubmult.abc import x
# from schubmult.symbolic import Symbol
# from schubmult.rings.free_algebra import *

# if __name__ == "__main__":
#     import sys
#     import itertools

#     n = int(sys.argv[1])
#     perms = Permutation.all_permutations(n)
#     #DualGroth = FreeAlgebra(GrothendieckBasis)

#     for perm1, perm2 in itertools.product(perms, repeat=2):
#         for p in range(1, n - 1):
#             if perm1.max_descent > p or perm2.max_descent > n - 1 - p:
#                 continue
#             perm1_min, perm1_J = perm1.coset_decomp(*list(range(1, p)))
#             perm2_min, perm2_J = perm2.coset_decomp(*list(range(1, n - 1 - p)))
#             result = AGx(perm1, p) * AGx(perm2_min, n - 1 - p)
#             real_result = AGx(perm1, p) * AGx(perm2, n - 1 - p)
#             for (perm, _), coeff in real_result.items():
#                 if pe
#                 assert (result.get((reduced_perm, n - 1), 0) - coeff).expand() == 0, f"Failed for {perm1}, {perm2}, {p}: got {result.get((reduced_perm, n - 1), 0)}, expected {coeff}\n{perm=}\n{reduced_perm=}\n{J=}\n{perm1_J=}\n{perm2_J=}\n{J=}"