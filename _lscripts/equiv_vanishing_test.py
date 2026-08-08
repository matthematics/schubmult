from schubmult import *
from schubmult.symbolic import *
from schubmult.symbolic.common_polys import *

if __name__ == "__main__":
    import sys
    import itertools
    from schubmult.abc import y, z
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    z_dict = {z[i]: y[i] for i in range(50)}
    for perm1, perm2 in itertools.combinations(perms, 2):
        prd = DSx(perm1)*DSx(perm2, "z")
        for k, v in prd.items():
            if expand(v) == 0:
                continue
            # if perm2 <= k:
            #     assert expand(v.subs({z[i]: y[i] for i in range(50)})) != 0, f"Fabulous fail {perm1.trimcode} {perm2.trimcode} {k.trimcode}"
            if expand((DSx(perm2) * DSx(perm1, "z")).get(k, 0)) != 0:
                assert expand(efficient_subs(v, z_dict)) != 0, f"Fabulous fail {perm1.trimcode} {perm2.trimcode} {k.trimcode}"
        print(f"Immaculate {perm1,perm2}")
        