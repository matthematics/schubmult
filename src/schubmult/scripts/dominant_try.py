
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
    # cd = []
    # for i in range(2 * (n - 1), 0, -2):
    #     cd += [i]
    # TRANSPORT CRYSTAL
    dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms}
    #perms2n = {perm for perm in Permutation.all_permutations(2 * n - 1) if perm.bruhat_leq(uncode(cd))}
    for dom in dominant_graphs:
        if dom.perm.inv == 0:
            continue
        for perm in perms:
            if perm.inv == 0:
                continue
            try:
                cprd = ASx(perm, n-1).coproduct()
            except Exception as e:
                print(f"Error computing coproduct for {perm=}, {n=}: {e}")
                continue
            prin_rc = RCGraph.principal_rc(perm, n-1)
            weight = tuple(a - b for a, b in zip(prin_rc.length_vector, dom.length_vector))
            prin_hw = prin_rc.to_highest_weight()[0]
            hw_weight = tuple(a - b for a, b in zip(prin_hw.length_vector, dom.length_vector))
            
            for ((dm, _ ),(perm2, _)), v in cprd.items():
                if dm != dom.perm:
                    continue
                result = 0
                hw_set = set()
                max_hw = None
                for rc in RCGraph.all_rc_graphs(perm2, n-1):
                    
                    total_hw = rc.to_highest_weight()[0]
                    if not max_hw:
                        max_hw = total_hw
                    if total_hw != max_hw:
                        continue
                    rc_hw = (CrystalGraphTensor(dom,rc).to_highest_weight()[0].factors[1])
                    if (rc_hw, total_hw) in hw_set:
                        continue
                    rc_lw = (CrystalGraphTensor(dom,rc).to_lowest_weight()[0].factors[1])
                    # if rc_hw.is_dom_perm_yamanouchi(dom.perm, perm):
                    #     rc_lw = (CrystalGraphTensor(dom,rc).to_lowest_weight()[0].factors[1])
                    #     if tuple([a + b for a, b in zip(rc_lw.length_vector,dom.length_vector)]) != RCGraph.principal_rc(perm,n-1).length_vector:
                    #         continue
                    #     if (rc_hw, rc_lw) not in hw_set:
                    #         hw_set.add((rc_hw, rc_lw))
                    #         result += 1
                    if rc_hw.length_vector == hw_weight and rc_lw.length_vector == weight:
                        hw_set.add((rc_hw, total_hw))
                        result += 1
                if v == result:
                    print(f"Success {dom.perm=}, {perm2=} {perm=} {result=}")
                else:
                    print(f"Warning: mismatch! {(dom.perm, perm2, perm)}: expected {v}, got {result=}")
                    raise Exception("Mismatch found")
    # for perm1 in perms:
    #     dom = RCGraph.principal_rc(Permutation.w0(n), n-1)
    #     #rc_check = perm1 * (~dom.perm)
    #     rc_check = Permutation([])
    #     #for perm_top in perms:
    #     perm_top = dom.perm
    #     perm = (~rc_check) * perm_top
    #     if perm_top.inv != perm.inv + rc_check.inv:
    #         continue
    #     if not perm1.bruhat_leq(perm):
    #         continue
    #     result[perm] = result.get(perm, {})
    #     if perm_top.inv != perm.inv + rc_check.inv:
    #         continue
    #     for rc1 in RCGraph.all_rc_graphs(perm1, n-1):
    #         for perm2 in perms:
    #             if not perm2.bruhat_leq(perm):
    #                 continue
    #             if perm2.inv == 0:
    #                 continue
    #             if max(len(perm0) for perm0 in (Sx(perm1)*Sx(perm2)).keys()) > n:
    #                 continue
                
                
    #             for rc2 in RCGraph.all_rc_graphs(perm2, n-1):
    #                 # if tuple(a + b for a, b in zip(rc1.length_vector, rc2.length_vector)) != RCGraph.principal_rc(perm,n-1).length_vector:
    #                 #         continue
    #                 if rc2.is_dom_perm_yamanouchi(dom.perm, perm_top):
    #                     result[perm][(perm1,perm2)] = result[perm].get((perm1, perm2), 0) + 1
        
    #     matches = {}

    #     for perm in perms:
    #         for (perm1, perm2) in itertools.product(perms, perms):
    #             if perm1.inv + perm2.inv != perm.inv:
    #                 continue
    #             product = (Sx(perm1) * Sx(perm2))
    #             if max(len(perm0) for perm0 in product.keys()) > n:
    #                 continue
    #             if result.get(perm) is None:
    #                 if product.get(perm, 0) == 0:
    #                     print(f"Success {perm1=}, {perm2=} {perm=} {product.get(perm, 0)=}")
    #                 else:
    #                     print(f"Warning: mismatch! {(perm1, perm2)}: expected {product.get(perm, 0)}, got 0")
    #                 continue
    #             if product.get(perm, 0) == result[perm].get((perm1, perm2), 0):
    #                 print(f"Success {perm1=}, {perm2=} {perm=} {result[perm].get((perm1, perm2), 0)=}")
    #             else:
    #                 print(f"Warning: mismatch! {(perm1, perm2)}: expected {product.get(perm, 0)}, got {result[perm].get((perm1, perm2), 0)=}")

                
            #     matches[(perm1, perm2)] = True
            # else:
            #     matches[(perm1, perm2)] = False
            #     print(f"Warning: mismatch! {(perm1, perm2)}: expected {product.get((perm1, perm2), 0)}, got {result.get((perm1, perm2), 0)}")
                
            #     #input()
            # #print(f"Matches: {matches}")
            
            