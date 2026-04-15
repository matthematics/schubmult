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
        
        length = 100
        # partition1 = tuple((~(perm1.strict_mul_dominant(length))).trimcode)
        # partition2 = tuple((~(perm2.strict_mul_dominant(length))).trimcode)
        real_prod = Sx(perm1) * Sx(perm2)
        #length = max(len(perm) - 1 for perm in real_prod)
        #max(len(perm1), len(perm2)) - 1
        #test_prod = r.schub_elem(perm1, length, partition=partition1) * r.schub_elem(perm2, length, partition=partition2)
        test_prod = r.schub_elem(perm1, length, partition=tuple(range(n - 1, 0, -1))) * r.schub_elem(perm2, length, partition=tuple(range(n - 1, 0, -1))) 
        

        for perm, coeff in real_prod.items():
            test_prod -= coeff * r.schub_elem(perm, length)#, partition=tuple(range(n - 1, 0, -1)))
        #assert test_prod.to_rc_graph_ring_element().almosteq(rc_ring.zero), f"Failed product test for {perm1} and {perm2}: got {test_prod}"    
        #test_prod_dct = {k: v for k, v in test_prod.items() if all(len(rc.perm) <= length + 1 for rc, coeff in r(k).to_rc_graph_ring_element().items() if coeff != 0)}
        #assert test_prod.almosteq(r.zero), f"Failed product test for {perm1} and {perm2}: got {test_prod}"
        #test_elem = r.from_tensor_dict(test_prod_dct, length)
        test_elem = test_prod
        assert test_elem.almosteq(r.zero), f"Failed product test for {perm1} and {perm2}: got {test_elem}"
        #assert all(test_elem[k] == 0 for k in test_elem if all(len(rc.perm) <= length + 1 for rc, coeff in r(k).to_rc_graph_ring_element().items() if coeff != 0))), f"Failed product test for {perm1} and {perm2}: got {test_elem}"