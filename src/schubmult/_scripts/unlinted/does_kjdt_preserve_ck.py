from schubmult import *

if __name__ == "__main__":
    import sys
    import random
    import itertools

    n = int(sys.argv[1])
    max_iter = int(sys.argv[2])

    perms = Permutation.all_permutations(n)

    skew_shapes = set()


    for perm in perms:
        for rc in RCGraph.all_hw_rcs(perm, len(perm.trimcode)):
            skew_shapes.add(NilPlactic.from_word(rc.perm_word).shape)
            nilp = NilPlactic.from_word(rc.perm_word)
            orig_nilp = nilp
            old_nilp = nilp
            num_slides = 0
            for _ in range(max_iter):
                # print("Current NilPlactic:")
                # print(nilp)
                # print("Original NilPlactic:")
                # print(orig_nilp)
                new_nilp = NilPlactic.from_word(nilp.row_word)
                assert new_nilp == orig_nilp, f"NilPlactic changed after JDT slide: {old_nilp} -> {nilp}"
                outer_corners = list(nilp.iter_outer_corners)
                if outer_corners:
                    outer_corner = random.choice(outer_corners)
                    print("Sliding outer corner:", outer_corner, "on")
                    print(nilp)
                    old_nilp = nilp
                    nilp = nilp.up_jdt_slide(*outer_corner)
                else:
                    inner_corners = list(nilp.iter_inner_corners)
                    if inner_corners:
                        inner_corner = random.choice(inner_corners)
                        print("Sliding inner corner:", inner_corner, "on")
                        print(nilp)
                        old_nilp = nilp
                        nilp = nilp.down_jdt_slide(*inner_corner)
                    else:
                        print(f"NilPlactic has no corners to slide: {nilp} only did {num_slides} slides")
                        break
                num_slides += 1
            nilp = nilp.rectify()
            assert nilp == orig_nilp, f"NilPlactic changed after rectify"

    
    for perm1, perm2 in itertools.product(perms, repeat=2):
        mulpy = Sx(perm1) * Sx(perm2)
        
        for w in mulpy:
            shape = tuple([a for a in RCGraph.principal_rc(w).to_highest_weight()[0].length_vector if a != 0])
            for rc in RCGraph.all_hw_rcs(perm2, len(perm2.trimcode)):
                shape2 = tuple([a for a in rc.length_vector if a != 0])
                tabs = NilPlactic.all_skew_ed_tableaux(shape, perm1, shape2)
                # skew_shapes.add(NilPlactic.from_word(rc.perm_word).shape)
                for tab in tabs:
                    nilp = tab
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
                                print(f"NilPlactic has no corners to slide: {nilp}\norig {orig_nilp} \nonly did {num_slides} slides")
                                break
                        num_slides += 1
                    nilp = nilp.rectify()
                    assert nilp == orig_nilp, f"NilPlactic changed after rectify"
        
        