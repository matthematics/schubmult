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
        
        rc_prin = RCGraph.principal_rc(perm, n - 1)
        for perm2 in perms:
            if perm2 == perm:
                continue
            for rc in RCGraph.all_rc_graphs(perm2, n - 1):
                if rc.to_highest_weight()[0].length_vector == rc_prin.to_highest_weight()[0].length_vector and rc.to_lowest_weight()[0].length_vector == rc_prin.length_vector:
                    fail = True
                    raise Exception(f"Failure in RCGraph insertion test {rc=} {rc_prin=}")
