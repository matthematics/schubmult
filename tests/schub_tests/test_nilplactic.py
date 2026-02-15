def test_ed_insert_rsk():
    from schubmult import NilPlactic, RCGraph, Plactic
    rc = RCGraph([(3,1),(2,),(4,)])
    nilp, plac = NilPlactic.ed_insert_rsk(rc.perm_word, rc.compatible_sequence)
    
    hw_rc = nilp.hw_rc(len(rc))
    rc_hw, raise_seq = rc.to_highest_weight()
    assert hw_rc == rc_hw, f"HW RC mismatch: {hw_rc} vs {rc_hw}"
    plac = plac.transpose()

    high_plac = Plactic.yamanouchi([a for a in rc_hw.length_vector if a != 0])
    high_plac = high_plac.reverse_raise_seq(raise_seq)
    assert plac == high_plac, f"Plac mismatch: {plac} vs {high_plac}"