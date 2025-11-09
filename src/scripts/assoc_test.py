from sympy import init_printing, pretty_print

from schubmult import RCGraphRing

# IS THIS ASSOCIATIVE?
# need associativity
rc_ring = RCGraphRing()


if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    from symengine import S
    from sympy import pretty_print

    from schubmult import CrystalGraphTensor, MonomialBasis, Permutation, PolynomialAlgebra, RCGraph, SchubertBasis, SchubertPolyBasis, WordBasis, uncode
    from schubmult.abc import x
    #from schubmult.utils.perm_utils import artin_sequences

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    
    for perm1 in perms:
        touched = set()
        val0 = S.Zero
        val1 = S.Zero
        for rc1 in RCGraph.all_rc_graphs(perm1, n - 1):
            
            rc2 = rc1.transpose().resize(n-1).to_highest_weight()[0]
            if rc2 in touched:
                continue
            cut_spot = (n - 1)//2
            for rc02 in rc2.full_crystal:
                left, right = rc02.vertical_cut(cut_spot)
                val0 += left.transpose().polyvalue(x) * right.transpose().polyvalue(x)
                val1 += rc02.transpose().polyvalue(x)
                pretty_print(rc02.transpose())
                pretty_print(CrystalGraphTensor(left.transpose(),right.transpose()))
        print(f"{val0=}")
        print(f"{val1=}")

