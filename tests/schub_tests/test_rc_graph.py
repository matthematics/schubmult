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

def test_computes_double_schub():
    from schubmult import RCGraph, Permutation, DSx
    from schubmult.symbolic import S, expand
    from schubmult.abc import x, y
    n = 5
    perms = Permutation.all_permutations(n)
    
    for perm in perms:
        rc_graphs = RCGraph.all_rc_graphs(perm, len(perm.trimcode))
        the_sum = sum([rc.polyvalue(x,y) for rc in rc_graphs])
        assert expand(the_sum - DSx(perm).expand()) == S.Zero


def test_little_bump_zero_equivalent():
    from schubmult import RCGraph, Permutation
    n = 6
    perms = Permutation.all_permutations(n)
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
            if rc.perm.inv == 0 or len(rc[-1]) != 0:
                continue
            assert rc.little_bump_zero() == rc.zero_out_last_row()


def test_tableau_decomp():
    from schubmult import RCGraph

    rc = RCGraph(((5, 4, 3, 2, 1), (6, 5, 3, 2), (6, 4), (5,), ()))

    expected = (RCGraph([(5, 4, 3, 2, 1)]), RCGraph([(4, 3, 2, 1)]), RCGraph(((3, 1), (2,))), RCGraph([()]))

    assert rc.tableau_decomp() == expected


def test_pull_out_empty_row():
    from schubmult import RCGraph, Permutation
    n = 5
    perms = Permutation.all_permutations(n)
    for perm in perms:
        lower_perms = {}
        for i in range(2, len(perm.trimcode)):
            lower_perms = {}
            for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
                if rc.perm.inv == 0 or len(rc[i-1]) != 0:
                    continue
                pullout = rc.pull_out_row(i)
                lower_perms[pullout.perm] = lower_perms.get(pullout.perm, set())
                lower_perms[pullout.perm].add(pullout)

        for perm2 in lower_perms:
            assert lower_perms[perm2] == RCGraph.all_rc_graphs(perm2, len(perm.trimcode) - 1), f"Error: missing RC graphs for permutation {perm2} from pull out of {perm} at row {i}."    