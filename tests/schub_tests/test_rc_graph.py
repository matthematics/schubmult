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

def test_from_reduced_compatible():
    from schubmult import RCGraph

    redc = [(2,3,2), (1,2,2)]
    expected = RCGraph([(2,),(3,2),()])
    
    assert RCGraph.from_reduced_compatible(*redc) == expected
    assert RCGraph.from_reduced_compatible(*redc, 6) == expected.resize(6)

def test_as_reduced_compatible():
    from schubmult import RCGraph

    rc = RCGraph([(2,),(3,2),()])
    expected = ((2,3,2), (1,2,2))

    assert rc.as_reduced_compatible() == expected

def test_loc_of_inversion():
    from schubmult import RCGraph

    rc = RCGraph([(2,),(3,2),()])
    expected = (1,2)

    assert rc.loc_of_inversion(3,4) == expected

    expected = (2,1)

    assert rc.loc_of_inversion(2,3) == expected

def test_index_of_inversion():
    from schubmult import RCGraph

    rc = RCGraph([(2,),(3,2),()])
    expected = 0

    assert rc.index_of_inversion(3,4) == expected

    expected = 2

    assert rc.index_of_inversion(2,3) == expected

def test_little_bump():
    from schubmult import RCGraph

    rc = RCGraph([(4,2,),(3,2),()])
    expected = RCGraph([(5,2),(3,2),(),(),()])

    assert rc.little_bump(2,5) == expected

    expected = RCGraph([(4,3),(4,2),(),()])

    assert rc.little_bump(3,4) == expected

def test_tableau_decomp():
    from schubmult import RCGraph

    rc = RCGraph(((5, 4, 3, 2, 1), (6, 5, 3, 2), (6, 4), (5,), ()))

    expected = (RCGraph([(5, 4, 3, 2, 1)]), RCGraph([(4, 3, 2, 1)]), RCGraph(((3, 1), (2,))), RCGraph([()]))

    assert rc.tableau_decomp() == expected