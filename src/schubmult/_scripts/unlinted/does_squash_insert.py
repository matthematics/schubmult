from schubmult import *
from schubmult.schub_lib.crystal_graph import CrystalGraphTensor
from schubmult.utils.schub_lib import elem_sym_perms
from sympy import pretty_print, S

P = PlacticAlgebra()

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    # cauchy
    # cauchy_prod = (Sx@r).one
    
    # for new_n in range(n - 1, 0, -1):
    #     new_cauchy_prod = (Sx@r).zero
    #     for deg in range(new_n + 1):
    #         if deg < new_n:
    #             elem_perm = uncode([0] * (new_n - deg) + [1] * deg)
    #             for rc in RCGraph.all_rc_graphs(elem_perm, new_n):
    #             cauchy_prod += sum(r(rc)for rc in RCGraph.all_rc_graphs(, new_n))
    dom_rcs = [RCGraph.principal_rc(perm, n - 1) for perm in perms if perm.is_dominant]
    for perm1, dom_rc in itertools.product(perms, dom_rcs):
        pass
        # for k in range(1, n):
        #     sputnik = elem_sym_perms(perm, k, k)
        #     for p in range(1, k + 1):
        #         for rc1 in RCGraph.all_rc_graphs(perm, n):
        #             rc1_t = rc1.transpose()
        #             if len(perm.trimcode) > k:
                        
        #                 rc1_tt = rc1_t.rowrange(len(rc1_t) - k, len(rc1_t)).transpose(k)
        #             else:
        #                 rc1_tt = rc1.resize(k)
        #             for rc2 in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1] * p), k):
        #                 new_rc = rc1_tt.squash_product(rc2)
        #                 if len(perm.trimcode) > k:
        #                     fatbank = RCGraph([*RCGraph([*rc1_t[:len(rc1_t) - k]]).transpose(n), *new_rc.shiftup(len(rc1_t) - k)])
        #                 else:
        #                     fatbank = new_rc
        #                 pretty_print(fatbank)
        #                 print(perm)
        #                 print(sputnik)
        #                 assert (fatbank.perm, p) in sputnik
                        