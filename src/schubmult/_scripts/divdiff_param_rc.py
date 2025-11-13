from sympy import init_printing, pretty_print





if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    from schubmult import x
    from symengine import S
    from sympy import pretty_print

    from schubmult import Permutation, RCGraph, RCGraphRing
    # from schubmult.utils.perm_utils import artin_sequences

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    # sKEW DIV DIFF WEIGHT
    def elem_rc(k):
        return RCGraph([(i,) for i in range(1, k + 1)])
    
    for perm1 in perms:
        start_rcs = [elem_rc(k) for k in range(1, n)]
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
                    spigot = rc_list[i].divdiff_desc(d)
                    for rc00 in spigot:
                        new_rc_list = [*[rc0 for rc0 in rc_list[:i]], rc00, *[rc_list[j].weight_reflection(d) for j in range(i + 1, len(rc_list))]]
                        upd_rc_stack.append(new_rc_list)
                    
            rc_stack = upd_rc_stack
        print(f"{(perm1 * Permutation.w0(n)).trimcode}")
        pretty_print([[bacon for i, bacon in enumerate(rc_list)] for rc_list in rc_stack])            
        
