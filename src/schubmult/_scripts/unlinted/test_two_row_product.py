from schubmult import *
from schubmult.utils.perm_utils import weak_compositions


if __name__ == "__main__":
    import sys
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra, SchubertPolyBasis, Schub
    import itertools

    num_parts = int(sys.argv[1])
    max_part_size = int(sys.argv[2])

    comps = weak_compositions(num_parts, max_part_size)

    r = DualRCGraphRing()

    for comp1, comp2 in itertools.product(comps, repeat=2):
        perm1 = uncode(comp1)
        perm2 = uncode(comp2)
        produ = Schub(perm1,2) * Schub(perm2,2)
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, num_parts), RCGraph.all_rc_graphs(perm2, num_parts)):
            result = next(iter(r(rc1) * r(rc2)))
            assert produ.get((result.perm, 2), 0) != 0, f"Failure for comp {comp1}, {comp2}, got {result.perm=}, {produ=}, {result=}"