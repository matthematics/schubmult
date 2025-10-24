from sympy import init_printing, pretty_print

from schubmult.rings.rc_graph_ring import RCGraphRing
from schubmult.rings.rc_graph import RCGraph
from schubmult.rings.ck_ring import CoxeterKnuthKey

# IS THIS ASSOCIATIVE?
# need associativity
rc_ring = RCGraphRing()


if __name__ == "__main__":
    # test module functionality

    import sys

    from schubmult import FA, ASx, Permutation, SchubertBasis, WordBasis, uncode
    from schubmult.abc import x
    from schubmult.rings import MonomialBasis, PolynomialAlgebra, SchubertPolyBasis
    from schubmult.utils.perm_utils import artin_sequences

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3

    perms = Permutation.all_permutations(n)
    modfull = 0
    dct = {}
    deg = 6
    
    fail = False
    failures = []
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm):
            ck_key = CoxeterKnuthKey.from_rc_graph(rc)
            for i in range(1, n + 1):
                print("Inserting box", i, "into RCGraph with perm", perm)
                rc_insert = ck_key.monk_insert(i)
                pretty_print(rc_insert)
                if rc_insert is None:
                    fail = True
                    failures.append((rc, i))

    if fail:
        print("Some insertions failed")
        for rc, i in failures:
            print(f"Failed insertion of box {i} into RCGraph:")
            pretty_print(rc.normalize())
    else:
        print("All insertions succeeded")