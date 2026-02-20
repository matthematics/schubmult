def test_ed_insert_rsk():
    from schubmult import NilPlactic, RCGraph, Plactic
    rc = RCGraph([(3,1),(2,),(4,)])
    nilp, plac = NilPlactic.ed_column_insert_rsk(rc.perm_word, rc.compatible_sequence)
    
    hw_rc = nilp.hw_rc(len(rc))
    rc_hw, raise_seq = rc.to_highest_weight()
    assert hw_rc == rc_hw, f"HW RC mismatch: {hw_rc} vs {rc_hw}"
    
    high_plac = Plactic.yamanouchi([a for a in rc_hw.length_vector if a != 0])
    high_plac = high_plac.reverse_raise_seq(raise_seq)
    assert plac == high_plac, f"Plac mismatch: {plac} vs {high_plac}"

def test_does_kjdt_preserve_ck():
    import random
    from schubmult import NilPlactic, RCGraph, Permutation

    n = 6
    max_iter = 20
    seed = 12515

    random.seed(seed)

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
                num_slides += 1
            assert num_slides > 0, f"NilPlactic had no corners to slide: {nilp}"