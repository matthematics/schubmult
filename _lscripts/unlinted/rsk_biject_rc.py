from schubmult import *

if __name__ == "__main__":
    import sys
    from sympy import pretty_print

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    e = EGPlacticRing()
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
            seq = []
            for i, val in enumerate(rc.length_vector):
                seq.extend([i + 1] * val)
            nil, plac = NilPlactic.ed_insert_rsk(rc.perm_word, seq)
            plac = plac.transpose()
            elem = e((nil, len(rc)), plac)
            rc_elem = e.from_rc_graph(rc)
            assert elem.almosteq(rc_elem), f"Failed on {rc}\nExpected:\n{elem}\nGot:\n{rc_elem}"