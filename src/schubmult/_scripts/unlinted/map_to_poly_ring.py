from schubmult import *
from schubmult.symbolic import *

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    r = RCGraphRing()
    perms = Permutation.all_permutations(n)
    for perm in perms:
        
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            the_val = S.Zero
            for perm2 in perms:
                table = r.from_free_algebra_element(ASx(perm2, n - 1))
                the_val += table.get(rc, S.Zero) * Sx(perm2).expand()
            print(f"RC graph {rc} contributes {expand(the_val)} vs {rc.polyvalue(Sx.genset)}")