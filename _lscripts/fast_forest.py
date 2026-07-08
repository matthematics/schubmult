from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.utils.tuple_utils import pad_tuple
import argparse

if __name__ == "__main__":
    import itertools
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int)
    parser.add_argument("--old-way", action="store_true", default=False)
    args = parser.parse_args()

    n = args.n
    old_way = args.old_way
    ForestPoly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.trimcode) for perm in perms]
    if old_way:
        print("Eat a potato")
        exit(1)
    r = ForestRCGraphRing()
    for comp1, comp2 in itertools.product(comps, repeat=2):
        result = 0
        length = max(1, max(len(comp1), len(comp2)))
        if len(comp1) < length:
            comp11 = pad_tuple(comp1, length)
        else:
            comp11 = comp1[:length]
        if len(comp2) < length:
            comp22 = pad_tuple(comp2, length)
        else:
            comp22 = comp2[:length]
        for rc1, rc2 in itertools.product(r.forest_poly(comp11), r.forest_poly(comp22)):
            rc = next(iter((r._forest_brc_lookup(rc1) * r._forest_brc_lookup(rc2)).to_rc_graph_ring_element().resize(length)))
            #r._forest_brc_lookup(rc1) * self._forest_brc_lookup(rc2)
            if rc.forest_weight == rc.length_vector:
                result += ForestPoly(*rc.length_vector)
        # true_result = ForestPoly(*comp1) * ForestPoly(*comp2)
        # assert result.almosteq(true_result), f"Mismatch for {comp1=}, {comp2=}: {result=} vs {true_result=}"
        print(f"Hot baby {comp1} {comp2}")