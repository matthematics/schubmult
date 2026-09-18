from schubmult import *
from schubmult.symbolic import S, expand
import itertools
from functools import cache

@cache
def pivot_transition(perm2, target_d=None):
    x = Gx.genset
    beta = Gx.beta
    build_groth = S.Zero
    if target_d is None:
        d = perm2.max_descent
    else:
        d = target_d
    if perm2.max_descent < d:
        return Gx(perm2).as_polynomial()
    pivots = perm2.pivots()
    
    for r in range(1, len(pivots) + 1):
        for pivot_set in itertools.combinations(pivots, r):
            pivot_set = set(pivot_set)
            ptrans = perm2.pivot_transition(pivot_set)
            print("Debug: ", perm2, pivot_set, ptrans)
            build_groth += (beta**(len(pivot_set) - 1)) * pivot_transition(ptrans, target_d=d)
    return build_groth
            
if __name__ == "__main__":
    
    import sys

    
    n = int(sys.argv[1])
    x = Gx.genset
    beta = Gx.beta
    for perm in Permutation.all_permutations(n):
        # When this file is executed as a module (-m), the local Permutation class
        # can differ from the one imported inside GrothendieckRing.new; pass one-line
        # notation explicitly to avoid isinstance mismatches.
        if perm.inv == 0:
            continue
        base_groth = Gx(perm)
        d = perm.max_descent
        zeroed_groth = base_groth.as_polynomial().subs(x[d], S.Zero)
        build_groth = pivot_transition(perm)
        diff = expand(build_groth.expand() - zeroed_groth.expand())
        assert diff == S.Zero, f"Failed for {perm} with build {build_groth.expand()} and zeroed {zeroed_groth.expand()}\nDifference: {diff}"
        print(f"Passed for {perm}")
