def test_rc_zeroing():
    from schubmult import RCGraph

    # Create an RC graph with a last row that can be zeroed
    rc = RCGraph([(3, 2),(3,),()])
    result = RCGraph([(3,1),(2,)])
    assert rc.zero_out_last_row() == result

def test_rc_zero_padding():
    from schubmult import RCGraph, BPD

    # Create an RC graph with a last row that can be zeroed
    result = {RCGraph([(3, 2),(3,),()]),RCGraph([(3,1),(2,),()]),RCGraph([(4,1),(3,),()]), RCGraph([(4,2),(3,),()])}
    rc = RCGraph([(3,1),(2,)])
    assert rc.right_zero_act() == result

def test_rc_product():
    from schubmult import RCGraph
    rc1 = RCGraph([(3,1), (2,)])
    rc2 = RCGraph([(2,1),(3,),()])

    expected = {RCGraph(((6, 2), (5,), (4, 3), (5,), ())): 1, RCGraph(((5, 3), (6,), (4, 3), (5,), ())): 1, RCGraph(((7, 2), (6,), (4, 3), (5,), ())): 1, RCGraph(((6, 5), (6,), (4, 3), (5,), ())): 1, RCGraph(((7, 1), (6,), (4, 3), (5,), ())): 1, RCGraph(((3, 2), (3,), (4, 3), (5,), ())): 1, RCGraph(((3, 2), (6,), (4, 3), (5,), ())): 1, RCGraph(((7, 5), (6,), (4, 3), (5,), ())): 1, RCGraph(((3, 1), (2,), (4, 3), (5,), ())): 1, RCGraph(((6, 1), (5,), (4, 3), (5,), ())): 1, RCGraph(((3, 2), (5,), (4, 3), (5,), ())): 1, RCGraph(((6, 3), (5,), (4, 3), (5,), ())): 1, RCGraph(((5, 1), (2,), (4, 3), (5,), ())): 1, RCGraph(((4, 3), (5,), (4, 3), (5,), ())): 1, RCGraph(((5, 2), (6,), (4, 3), (5,), ())): 1, RCGraph(((6, 2), (3,), (4, 3), (5,), ())): 1, RCGraph(((5, 1), (3,), (4, 3), (5,), ())): 1, RCGraph(((6, 1), (3,), (4, 3), (5,), ())): 1, RCGraph(((6, 1), (2,), (4, 3), (5,), ())): 1, RCGraph(((7, 3), (6,), (4, 3), (5,), ())): 1, RCGraph(((5, 2), (3,), (4, 3), (5,), ())): 1}

    prod = rc1.product(rc2)

    assert prod == expected