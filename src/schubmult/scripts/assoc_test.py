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
        if perm.inv == 0:
            continue
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            ck_key = CoxeterKnuthKey.from_rc_graph(rc)
            for i in range(1, n):
                for box_val in range(1, i + 1):
                    print(f"Inserting box {box_val=} for {i} into RCGraph with perm {perm}")
                    rc_insert = ck_key.monk_insert(box_val, i)
                    print("Before:")
                    pretty_print(rc.normalize())
                    print("After:")
                    for rc2 in rc_insert or []:
                        assert rc2.perm.inv == perm.inv + 1, f"Inserted RCGraph has wrong permutation: got {rc2.perm.inv}, want {perm.inv + 1}"
                        pretty_print(rc2)
                        print(tuple([tuple(row) for row in rc2]))
                    if rc_insert is None:
                        fail = True
                        failures.append((rc, i, box_val))

    if fail:
        print("Some insertions failed")
        for rc, i, box_val in failures:
            print(f"Failed insertion of box {box_val} for {i} into RCGraph:")
            pretty_print(rc.normalize())
    else:
        print("All insertions succeeded")