import itertools

import symengine
from joblib import Parallel, delayed
from sympy import together


def check_pair(perm1, perm2):
    from schubmult.symbolic import sympify_sympy, efficient_subs, S, expand
    from schubmult import GeneratingSet, DGx, Gx
    poly1 = (DGx(perm1) * DGx(perm2)).as_polynomial()
    poly2 = DGx(perm1).as_polynomial() * DGx(perm2).as_polynomial()
    diff = poly1 - poly2
    t = GeneratingSet("t")
    y = GeneratingSet("y")
    z = GeneratingSet("z")
    s = GeneratingSet("s")
    subs_dict = {z[i]: (t[i] - S.One)/Gx._beta for i in range(50)}
    subs_dict |= {y[i]: s[i]/(S.One - Gx._beta * s[i]) for i in range(50)}
    subs_dict2 = {t[i]: S.One + Gx._beta * z[i] for i in range(50)}
    subs_dict2 |= {s[i]: y[i]/(S.One + Gx._beta * y[i]) for i in range(50)}
    # exact (not probabilistic) zero test: clear denominators via together(), which
    # combines via GCD instead of guessing, then expand only the numerator through
    # symengine's C++ engine -- sympy's Expr.expand() never returns on these once
    # perm1/perm2 have several inversions, but symengine handles it in seconds
    numer = together(efficient_subs(efficient_subs(diff,subs_dict),subs_dict2), deep=False)
    assert expand(numer) == 0, f"Failed for {perm1}, {perm2}"
    return perm1, perm2


if __name__ == "__main__":
    import sys
    from schubmult import Permutation
    n = int(sys.argv[1])
    # -1 uses all available cores
    num_processors = int(sys.argv[2]) if len(sys.argv) > 2 else -1
    perms = Permutation.all_permutations(n)
    pairs = list(itertools.product(perms, repeat=2))
    results = Parallel(n_jobs=num_processors, verbose=10)(delayed(check_pair)(perm1, perm2) for perm1, perm2 in pairs)
    for perm1, perm2 in results:
        print(f"Success {perm1} {perm2}")