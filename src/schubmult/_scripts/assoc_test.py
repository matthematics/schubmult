from sympy import init_printing, pretty_print





if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    from schubmult.abc import x
    from symengine import S
    from sympy import pretty_print

    from schubmult import Permutation, RCGraph, RCGraphRing
    # from schubmultutils.perm_utils import artin_sequences

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    # sKEW DIV DIFF WEIGHT
    def elem_rc(k):
        return RCGraph([(i,) for i in range(1, k + 1)])
    from schubmult import CrystalGraphTensor
    for perm1 in perms:
        start_rcs = [elem_rc(k).resize(n-1) for k in range(1, n)]
        perm = perm1
        rc_stack = [start_rcs]
        while perm.inv > 0:
            d = max(perm.descents()) + 1
            perm = perm.swap(d - 1, d)
            upd_rc_stack = []
            for rc_list in rc_stack:
                
                for i, rc in enumerate(rc_list):
                    if rc_list[i].perm[d - 1] < rc_list[i].perm[d]:
                        continue
                    reflected = [rc_list[j].weight_reflection(d).resize(n-1) for j in range(i + 1, len(rc_list))]
                    if i + 1 < len(rc_list) - 2:
                        tens = CrystalGraphTensor(rc_list[i + 1], rc_list[i+2])
                        for j in range(i + 3, len(rc_list)):
                            tens = CrystalGraphTensor(tens, rc_list[j])
                        tens = tens.weight_reflection(d)
                        reflected = tens.factors
                    else:
                        reflected = rc_list[-1]
                    if reflected.raising_operator(d) is not None:
                        continue    
                    spigot = rc_list[i].divdiff_desc(d)
                    for rc00 in spigot:
                        tenstens = CrystalGraphTensor(rc00, reflected)
                        new_rc_list = [*[rc0 for rc0 in rc_list[:i]], rc00.resize(n-1), *(reflected.factors if isinstance(reflected, CrystalGraphTensor) else [reflected])]
                        upd_rc_stack.append(new_rc_list)
                        reflected0 = tenstens.lowering_operator(d)
                        while reflected0 is not None:
                            new_rc_list = [*[rc0 for rc0 in rc_list[:i]], *tenstens.factors]
                            upd_rc_stack.append(new_rc_list)
                            reflected0 = reflected0.lowering_operator(d)
                        
                    
            rc_stack = upd_rc_stack
        print(f"{( Permutation.w0(n)*(~perm1)).trimcode}")
        pretty_print([[bacon for i, bacon in enumerate(rc_list)] for rc_list in rc_stack])
        try:         
            assert len(rc_stack) == len(RCGraph.all_rc_graphs(Permutation.w0(n)*(~perm1), n - 1))
        except Exception as e:
            print(e)
            print(f"{len(rc_stack)} vs {len(RCGraph.all_rc_graphs(Permutation.w0(n)*(~perm1), n - 1))}")
            raise e
        
