from schubmult import *
from sympy import pretty_print, S

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    elems = r.one @ Sx.one
    for k in range(1, n):
        # print("Fat")
        # the_prod = Sx(perm1) * Sx(perm2)
        # plactic_elem = P.zero
        # for w, coeff in the_prod.items():
        #     plactic_elem += sum([coeff * P(rc.hw_tab_rep()[1]) for rc in RCGraph.all_rc_graphs(w, len(w.trimcode))])
        # plactic_elem2 = P.zero
        # for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n), RCGraph.all_rc_graphs(perm2, n)):
        #     # _, raise_seq1 = rc1.to_highest_weight()
        #     # _, raise_seq2 = rc2.to_highest_weight()
        new_elems = r.zero @ Sx.zero
        for p in range(0, k + 1):
            for (rc0, perm), coeff in elems.items():
                #r_elem = r_elem.resize(k)
                rc0 = rc0.resize(k)
                for elem_rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1] * p), k):
                    new_rc = rc0.squash_product(elem_rc)
                    new_elems += coeff * r(new_rc) @ (Sx(perm) * Sx(Sx.genset[n - k]**(k - p)))
                # for index_rc, coeff in index_elem.items():
                #     assert len(index_rc) == k
                #     if len(index_rc.perm) <= k + 1:
                #         new_elems[index_rc] = new_elems.get(index_rc, r.zero) + coeff * new_r_elem
        elems = new_elems

    #assert len(elems) == len(perms), f"Failed: {len(elems)} vs {len(perms)}"    
    w0 = Permutation.w0(n)
    for (rc, perm), fatness in elems.items():
        print(perm)
        pretty_print(rc)
        print(fatness)
        if fatness != 0:
            #assert fatness == 1
            assert perm == rc.perm * w0, f"{rc.perm * w0=}"
        # patsy_dict[rc1.perm] == patsy_dict.get(rc1.perm, r.zero)
        # try:
        #     assert coeff == 1
        #     #assert elem_set.almosteq(r.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(rc.perm*w0, n - 1),1))), f"Failed on {rc}"
            
        # except AssertionError as e:
        #     pretty_print(rc)
        #     pretty_print(elem_set)
        #     raise