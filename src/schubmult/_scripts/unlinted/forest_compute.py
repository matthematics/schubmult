from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic import *

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    for perm in perms:
        if perm.inv == 0:
            continue
        polystink = {}
        for rc in RCGraph.all_rc_graphs(perm):
            inv = rc.forest_invariant
            polystink[inv] = polystink.get(inv, S.Zero) + rc.polyvalue(Sx.genset)
        for inv in polystink:
            # print(expand(polystink[inv] - Forest(*inv.forest.code).expand()) )
            # print(inv)
            # print(inv.forest.code)
            # print(Forest(*inv.forest.code).expand())
            assert expand(polystink[inv] - Forest(*inv.forest.code, 0).expand()) == 0