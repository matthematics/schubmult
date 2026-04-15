from schubmult import *

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    r = BoundedRCFactorAlgebra()
    rc_ring = RCGraphRing()

    for perm1, perm2 in itertools.product(perms, repeat=2):
        if perm1.inv == 0 or perm2.inv == 0:
            continue
        length = max(len(perm1.trimcode), len(perm2.trimcode))
        test_prod = r.schub_elem(perm1, length) * r.schub_elem(perm2, length)
        real_prod = Sx(perm1) * Sx(perm2)

        for perm, coeff in real_prod.items():
            test_prod -= coeff * r.schub_elem(perm, length)
        #assert test_prod.to_rc_graph_ring_element().almosteq(rc_ring.zero), f"Failed product test for {perm1} and {perm2}: got {test_prod}"    
        assert test_prod.almosteq(r.zero), f"Failed product test for {perm1} and {perm2}: got {test_prod}"