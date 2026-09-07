from schubmult import *
import itertools

import symengine
from joblib import Parallel, delayed
from sympy import together

from schubmult.symbolic import sympify_sympy


def check_pair(perm1, perm2):
    poly1 = (DGx(perm1) * DGx(perm2)).as_polynomial()
    poly2 = DGx(perm1).as_polynomial() * DGx(perm2).as_polynomial()
    diff = poly1 - poly2
    # exact (not probabilistic) zero test: clear denominators via together(), which
    # combines via GCD instead of guessing, then expand only the numerator through
    # symengine's C++ engine -- sympy's Expr.expand() never returns on these once
    # perm1/perm2 have several inversions, but symengine handles it in seconds
    numer, _ = together(sympify_sympy(diff), deep=False).as_numer_denom()
    assert symengine.expand(symengine.sympify(numer)) == 0, f"Failed for {perm1}, {perm2}"
    return perm1, perm2


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    # -1 uses all available cores
    num_processors = int(sys.argv[2]) if len(sys.argv) > 2 else -1
    perms = Permutation.all_permutations(n)
    pairs = list(itertools.product(perms, repeat=2))
    results = Parallel(n_jobs=num_processors, verbose=10)(delayed(check_pair)(perm1, perm2) for perm1, perm2 in pairs)
    for perm1, perm2 in results:
        print(f"Success {perm1} {perm2}")