from schubmult import *
from schubmult.rings.polynomial_algebra import *
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
    comps = [tuple(perm.pad_code(n-1)) for perm in perms]
    if old_way:
        print("Eat a potato")
        exit(1)
    r = ForestRCGraphRing()
    for comp1, comp2 in itertools.product(comps, repeat=2):
        result = 0
        for rc1, rc2 in itertools.product(r.forest_poly(comp1), r.forest_poly(comp2)):
            rc = next(iter((r._forest_brc_lookup(rc1) * r._forest_brc_lookup(rc2)).to_rc_graph_ring_element().resize(n - 1)))
            #r._forest_brc_lookup(rc1) * self._forest_brc_lookup(rc2)
            if rc.forest_weight == rc.length_vector:
                result += ForestPoly(*rc.length_vector)
        true_result = ForestPoly(*comp1) * ForestPoly(*comp2)
        assert result.almosteq(true_result), f"Mismatch for {comp1=}, {comp2=}: {result=} vs {true_result=}"
        print(f"Hot baby {comp1} {comp2}")