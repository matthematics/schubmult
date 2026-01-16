def test_bpd_zeroing():
    from schubmult import RCGraph, BPD

    # Create an RC graph with a last row that can be zeroed
    rc = RCGraph([(3, 2),(3,),()])
    bpd = BPD.from_rc_graph(rc)
    result = BPD.from_rc_graph(RCGraph([(3,1),(2,)]))

    assert bpd.zero_out_last_row() == result

def test_bpd_zero_padding():
    from schubmult import RCGraph, BPD

    # Create an RC graph with a last row that can be zeroed
    result = {BPD.from_rc_graph(RCGraph([(3, 2),(3,),()])),BPD.from_rc_graph(RCGraph([(3,1),(2,),()])),BPD.from_rc_graph(RCGraph([(4,1),(3,),()])), BPD.from_rc_graph(RCGraph([(4,2),(3,),()]))}
    rc = RCGraph([(3,1),(2,)])
    bpd = BPD.from_rc_graph(rc)
    assert bpd.right_zero_act() == result

def test_bpd_product():
    from schubmult import RCGraph, BPD
    rc1 = RCGraph([(3,1), (2,)])
    rc2 = RCGraph([(2,1),(3,),()])

    bpd1 = BPD.from_rc_graph(rc1)
    bpd2 = BPD.from_rc_graph(rc2)
    expected_rc = {RCGraph(((6, 2), (5,), (4, 3), (5,), ())): 1, RCGraph(((5, 3), (6,), (4, 3), (5,), ())): 1, RCGraph(((7, 2), (6,), (4, 3), (5,), ())): 1, RCGraph(((6, 5), (6,), (4, 3), (5,), ())): 1, RCGraph(((7, 1), (6,), (4, 3), (5,), ())): 1, RCGraph(((3, 2), (3,), (4, 3), (5,), ())): 1, RCGraph(((3, 2), (6,), (4, 3), (5,), ())): 1, RCGraph(((7, 5), (6,), (4, 3), (5,), ())): 1, RCGraph(((3, 1), (2,), (4, 3), (5,), ())): 1, RCGraph(((6, 1), (5,), (4, 3), (5,), ())): 1, RCGraph(((3, 2), (5,), (4, 3), (5,), ())): 1, RCGraph(((6, 3), (5,), (4, 3), (5,), ())): 1, RCGraph(((5, 1), (2,), (4, 3), (5,), ())): 1, RCGraph(((4, 3), (5,), (4, 3), (5,), ())): 1, RCGraph(((5, 2), (6,), (4, 3), (5,), ())): 1, RCGraph(((6, 2), (3,), (4, 3), (5,), ())): 1, RCGraph(((5, 1), (3,), (4, 3), (5,), ())): 1, RCGraph(((6, 1), (3,), (4, 3), (5,), ())): 1, RCGraph(((6, 1), (2,), (4, 3), (5,), ())): 1, RCGraph(((7, 3), (6,), (4, 3), (5,), ())): 1, RCGraph(((5, 2), (3,), (4, 3), (5,), ())): 1}

    expected_bpd = {BPD.from_rc_graph(rc): coeff for rc, coeff in expected_rc.items()}
    prod = bpd1.product(bpd2)

    assert prod == expected_bpd

def test_gao_huang():
    from schubmult import RCGraph, BPD, Permutation

    n = 5
    perms = Permutation.all_permutations(n)
    
    for perm in perms:
        for length in range(len(perm.trimcode), n):
            for rc in RCGraph.all_rc_graphs(perm, length):
                rc_test = BPD.from_rc_graph(rc).to_rc_graph()
                assert rc == rc_test                
            for bpd in BPD.all_bpds(perm, length):
                bpd_test = BPD.from_rc_graph(bpd.to_rc_graph())
                assert bpd == bpd_test

            assert len(BPD.all_bpds(perm, length)) == len(RCGraph.all_rc_graphs(perm, length)), f"Error: Number of BPDs and RC graphs do not match for permutation {perm} with length {length}"