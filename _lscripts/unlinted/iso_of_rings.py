from schubmult import *

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print
    r = RCGraphRing()
    e = EGTensorRing()
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    for perm1, perm2 in itertools.product(perms, repeat=2):
        for len1, len2 in itertools.product(range(len(perm1.trimcode), n), range(len(perm2.trimcode), n)):
            for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, len1), RCGraph.all_rc_graphs(perm2, len2)):
                rc_elem = r(rc1) * r(rc2)
                eg_tensor_elem = e.from_rc_graph(rc1) * e.from_rc_graph(rc2)
                try:
                    assert rc_elem.almosteq(eg_tensor_elem.to_rc_graph_ring_element()), f"Failed for {rc1} and {rc2}"
                except AssertionError as ae:
                    print(ae)
                    pretty_print(rc_elem)
                    pretty_print(eg_tensor_elem.to_rc_graph_ring_element())
                    raise
    # rc = RCGraph([(4,1),(3,2),(),(4,)])
    # rc2 = RCGraph([(1,),(4,3,2),(4,),()])
    # rc_elem = r(rc) * r(rc2)
    # assert rc_elem.almosteq((e.from_rc_graph(rc) * e.from_rc_graph(rc2)).to_rc_graph_ring_element())