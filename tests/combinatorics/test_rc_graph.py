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

def test_full_crystal_is_key():
    from schubmult import RCGraph, Permutation, PolynomialAlgebra, Sx
    from schubmult.rings.polynomial_algebra import KeyPolyBasis
    from schubmult.symbolic import expand, S
    n = 5
    perms = Permutation.all_permutations(n)
    Key = PolynomialAlgebra(KeyPolyBasis(Sx.genset))

    for perm in perms:
        # num_extremal = 0
        # hw = set()
        for rc in RCGraph.all_hw_rcs(perm, len(perm.trimcode)):
            assert expand(Key(rc.extremal_weight).expand() - sum([rc0.polyvalue(Sx.genset) for rc0 in rc.full_crystal])) == S.Zero, f"Error: full crystal of RC graph {rc} does not match Key polynomial expansion of its length vector."
        #assert len(hw) == num_extremal, f"Error: number of extremal RC graphs for permutation {perm} does not match number of distinct highest weights obtained from them."          

    
def test_squash_decompose():
    """Decompose an n-row RC graph into a pair of n-row RC graph in S_n and an n-grass."""
    import random
    from schubmult import uncode
    from schubmult.combinatorics.rc_graph import RCGraph
    seed = 250

    random.seed(seed)
    
    perm = uncode([2,0,4,2,3])

    n = 5

    rcs = list(RCGraph.all_rc_graphs(perm, n))
    
    rc = random.choice(rcs)
    regular, grass = rc.squash_decomp()

    while len(grass.perm.descents()) == 0 or regular.perm.inv == 0:
        rc = random.choice(rcs)
        regular, grass = rc.squash_decomp()

    assert len(regular.perm) <= n

    assert grass.perm.descents() == {n - 1}

    assert regular.squash_product(grass) == rc, f"Error: decomposition of {rc} into {regular} and {grass} does not satisfy the expected product relation."


def test_left_squash_equals_squash_product_for_equal_length_full_grass():
    """For equal-length full Grassmannian RCs, left_squash and squash_product give identical results."""
    from schubmult import RCGraph, Permutation

    # Test across several full Grassmannian RCs of various lengths
    for n in range(2, 5):
        perms = Permutation.all_permutations(n)
        for perm in perms:
            length = n - 1
            rcs = list(RCGraph.all_rc_graphs(perm, length))
            # Filter to full Grassmannian RCs (single descent at last position)
            grass_rcs = [rc for rc in rcs if rc.perm.inv == 0 or rc.perm.descents() == {length - 1}]
            
            # Test all pairs of equal-length Grassmannian RCs
            for rc1 in grass_rcs:
                for rc2 in grass_rcs:
                    if len(rc1) != len(rc2):
                        continue
                    
                    # For equal-length full Grassmannian RCs, these should be the same
                    left_result = rc1.left_squash(rc2)
                    product_result = rc1.squash_product(rc2)
                    
                    assert left_result == product_result, (
                        f"For equal-length full Grassmannian RCs:\n"
                        f"  rc1={rc1}\n"
                        f"  rc2={rc2}\n"
                        f"  left_squash gives:    {left_result}\n"
                        f"  squash_product gives: {product_result}"
                    )