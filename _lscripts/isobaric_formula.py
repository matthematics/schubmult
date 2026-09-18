# from schubmult import *
# from schubmult.symbolic import S, expand, Symbol
# from schubmult.abc import H
# import itertools
# from functools import cache

# def isobaric_as_schub(start, end, x, beta, schub_elem):
#     """Return the isobaric divided difference operator pi_i = partial_i(1 + beta x_{i+1}) as a Schubert element."""
#     from schubmult.symbolic.poly.variables import CustomGeneratingSet
#     from schubmult.rings.schubert.double_schubert_ring import DoubleSchubertRing
#     betan1 = Symbol("betbet")
#     bgn = CustomGeneratingSet([betan1] * len(x))
#     bg = CustomGeneratingSet([beta] * len(x))
#     bring =  DoubleSchubertRing(x,bg)
#     lst = [_ for _ in range(start, end + 1)]
#     for deg in range(end + 1 - start):
#     for clist in itertools.combinations(lst, end - start + 1):
#     this_elem = schub_elem.ring.divdiff(uncode([0] * (start - 1) + [end + 1 - start]), schub_elem)
#     if start  != 1:
#         raise ValueError("start must be 1 for isobaric_as_schub")
#     result = (bring([]) * this_elem) * H(end + 1 - start, start, x[1:], bgn)
#     return schub_elem.ring.one * result - beta * schub_elem
            
# if __name__ == "__main__":
    
#     import sys
#     from schubmult.rings.schubert.nil_hecke import NilHeckeRing
#     nh = NilHeckeRing(Sx.genset)
#     n = int(sys.argv[1])
#     x = Gx.genset
#     beta = Gx.beta
#     for perm in Permutation.all_permutations(n):
#         for i in range(1,n):
#             # When this file is executed as a module (-m), the local Permutation class
#             # can differ from the one imported inside GrothendieckRing.new; pass one-line
#             # notation explicitly to avoid isinstance mismatches.
#             if perm.inv == 0:
#                 continue
#             val = nh.g_isobaric(uncode([i])).apply(Sx(perm))
#             test_val = isobaric_as_schub(1, i, Sx.genset, beta, Sx(perm))
#             assert val.almosteq(test_val), f"Failed for {perm} with i={i}"
