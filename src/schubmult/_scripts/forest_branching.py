from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from schubmult.utils.tuple_utils import pad_tuple
from schubmult.rings.combinatorial.hw_rc_ring import HWRCGraphRing
from sympy import pretty_print

def snap_key(rc):
    return next(iter([rcc for rcc in RCGraph.all_rc_graphs(rc.perm, len(rc), weight=rc.extremal_weight)]))

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    Key = PolynomialAlgebra(KeyPolyBasis(Sx.genset))
    # KeyDual = FreeAlgebra(KeyBasis)

    for perm in perms:
        if perm.inv == 0:
            continue
        comp = pad_tuple(tuple(perm.trimcode), n - 1)
        forest_set = [rc for rc in RCGraph.all_rc_graphs(perm, n - 1) if rc.forest_weight == comp]

        for index in range(1, n - 1):
            fbranch = Forest(*comp).branch(index)
            tbranch = (Forest @ Forest).zero
            for rc in forest_set:
                rc1, rc2 = rc.vertical_cut(index)
                if rc1.forest_weight == rc1.length_vector and rc2.forest_weight == rc2.length_vector:
                    tbranch += Forest(*rc1.forest_weight) @ Forest(*rc2.forest_weight)
            assert fbranch.almosteq(tbranch), f"Failed for {perm} at index {index} with {len(forest_set)} RC graphs, got {tbranch} but expected {fbranch}"

        key_set = [rc for rc in RCGraph.all_rc_graphs(perm, n - 1) if rc.extremal_weight == comp]

        for index in range(1, n - 1):
            fbranch = Key(*comp).branch(index)
            tbranch = (Key @ Key).zero
            for rc in key_set:
                rc1, rc2 = rc.vertical_cut(index)
                if rc1.extremal_weight == rc1.length_vector and rc2.extremal_weight == rc2.length_vector:
                    tbranch += Key(*rc1.extremal_weight) @ Key(*rc2.extremal_weight)
            assert fbranch.almosteq(tbranch), f"Failed for {perm} at index {index} with {len(key_set)} RC graphs, got {tbranch} but expected {fbranch}"
