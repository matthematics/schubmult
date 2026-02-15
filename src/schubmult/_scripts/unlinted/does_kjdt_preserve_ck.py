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
                # print("Current NilPlactic:")
                # print(nilp)
                # print("Original NilPlactic:")
                # print(orig_nilp)
                new_nilp = NilPlactic.from_word(nilp.row_word)
                assert new_nilp == orig_nilp, f"NilPlactic changed after JDT slide: {nilp} -> {new_nilp}"
                outer_corners = list(nilp.iter_outer_corners)
                if outer_corners:
                    outer_corner = random.choice(outer_corners)
                    # print("Sliding outer corner:", outer_corner, "on")
                    # print(nilp)
                    nilp = nilp.up_jdt_slide(*outer_corner)
                else:
                    inner_corners = list(nilp.iter_inner_corners)
                    if inner_corners:
                        inner_corner = random.choice(inner_corners)
                        # print("Sliding inner corner:", inner_corner, "on")
                        # print(nilp)
                        nilp = nilp.down_jdt_slide(*inner_corner)
                    else:
                        print(f"NilPlactic has no corners to slide: {nilp} only did {num_slides} slides")
                num_slides += 1
        