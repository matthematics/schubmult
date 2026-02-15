from schubmult import *

if __name__ == "__main__":
    import sys
    import random

    n = int(sys.argv[1])
    max_iter = int(sys.argv[2])

    perms = Permutation.all_permutations(n)
    for perm in perms:
        for rc in RCGraph.all_hw_rcs(perm, len(perm.trimcode)):
            nilp = NilPlactic.from_word(rc.perm_word)
            orig_nilp = nilp
            num_slides = 0
            for _ in range(max_iter):
                new_nilp = NilPlactic.from_word(nilp.row_word)
                assert new_nilp == orig_nilp, f"NilPlactic changed after JDT slide: {nilp} -> {new_nilp}"
                outer_corners = list(nilp.iter_outer_corners)
                if outer_corners:
                    nilp = nilp.up_jdt_slide(*random.choice(outer_corners))
                else:
                    inner_corners = list(nilp.iter_inner_corners)
                    if inner_corners:
                        nilp = nilp.down_jdt_slide(*random.choice(inner_corners))
                    else:
                        print(f"NilPlactic has no corners to slide: {nilp} only did {num_slides} slides")
                num_slides += 1
        