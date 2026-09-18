if __name__ == "__main__":
    from schubmult import *
    import sys
    import itertools
    import numpy as np
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()

    for perm1, perm2 in itertools.product(perms, repeat=2):
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n-1), RCGraph.all_rc_graphs(perm2, n-1)):
            fung = r.zero
            prod_perms = Sx(perm1) * Sx(perm2)
            weight = (np.array(rc1.length_vector,dtype=int)+np.array(rc2.length_vector,dtype=int)).tolist()
            for perm, coeff in prod_perms.items():
                rc_graphs = RCGraph.all_rc_graphs(perm, n-1, weight=weight)
                fung += r.from_dict(dict.fromkeys(rc_graphs, coeff))
            print(fung)