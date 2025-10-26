
from schubmult.rings.ck_ring import CoxeterKnuthRing
import itertools

if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.crystal_graph import CrystalGraph, CrystalGraphTensor
    from schubmult.schub_lib.schub_lib import elem_sym_perms
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms}

    
    for perm in perms:
        if perm.inv == 0:
            continue
        result = {}
        for perm1 in perms:
            if not perm1.bruhat_leq(perm):
                continue
            dom = RCGraph.principal_rc(Permutation.w0(n-1), n-1)
            rc_check = perm1 * (~dom.perm)
            perm_top = (~rc_check) * perm
            if perm_top.inv != perm.inv + rc_check.inv:
                continue
            for rc1 in RCGraph.all_rc_graphs(perm1, n-1):
                for perm2 in perms:
                    if not perm2.bruhat_leq(perm):
                        continue
                    if perm2.inv == 0:
                        continue
                    if max(len(perm0) for perm0 in (Sx(perm1)*Sx(perm2)).keys()) > n:
                        continue
                    
                    
                    for rc2 in RCGraph.all_rc_graphs(perm2, n-1):
                        if tuple(a + b for a, b in zip(rc1.length_vector, rc2.length_vector)) != RCGraph.principal_rc(perm,n-1).length_vector:
                                continue
                        if rc1.is_dom_perm_yamanouchi(dom.perm, perm_top) and rc2.is_dom_perm_yamanouchi(dom.perm, perm_top) and tuple([a + b for a, b in zip(rc1.length_vector, rc2.length_vector)]) == RCGraph.principal_rc(perm_top,n-1).length_vector:
                            result[(perm1,perm2)] = result.get((perm1, perm2), 0) + 1
            
        matches = {}

        for (perm1, perm2) in itertools.product(perms, perms):
            if perm1.inv + perm2.inv != perm.inv:
                continue
            product = (Sx(perm1) * Sx(perm2))
            if max(len(perm0) for perm0 in product.keys()) > n:
                continue
            assert product.get((perm1, perm2), 0) == result.get((perm1, perm2), 0)
            print(f"Success {perm1=}, {perm2=} {perm=}")

                
            #     matches[(perm1, perm2)] = True
            # else:
            #     matches[(perm1, perm2)] = False
            #     print(f"Warning: mismatch! {(perm1, perm2)}: expected {product.get((perm1, perm2), 0)}, got {result.get((perm1, perm2), 0)}")
                
            #     #input()
            # #print(f"Matches: {matches}")
            
            