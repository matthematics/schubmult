from schubmult import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    # KEY HW IDEAL

    def ideal_decomp(rc_elem):
        rc_part = r.zero
        ideal_part = r.zero
        for rc, coeff in rc_elem.items():
            top = r(rc).project()
            rc_part += coeff * top
            ideal_part += coeff * (r(rc) - top)
        return rc_part, ideal_part

    def ideal_mul(rc1, rc2):
        rc_part1, ideal_part1 = ideal_decomp(rc1)
        rc_part2, ideal_part2 = ideal_decomp(rc2)
        return rc_part1 * rc_part2, ideal_part1 * rc_part2 + rc_part1 * ideal_part2 + ideal_part1 * ideal_part2

    for perm1, perm2 in itertools.product(perms, repeat=2):
        for lengtht1 in range(max(1, len(perm1.trimcode)), n):
            for lengtht2 in range(max(1, len(perm2.trimcode)), n):
                for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, lengtht1), RCGraph.all_rc_graphs(perm2, lengtht2)):
                    ideal_prod = ideal_mul(r(rc1), r(rc2))
                    true_prod = r(rc1) * r(rc2)
                    assert ideal_prod[0].almosteq(ideal_decomp(true_prod)[0]), f"Failed for {rc1} and {rc2}, got {ideal_prod[0]-ideal_decomp(true_prod)[0]}"
        