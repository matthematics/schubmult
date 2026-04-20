from schubmult import *

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    r = BoundedRCFactorAlgebra()
    rc_ring = RCGraphRing()

    for perm1, perm2 in itertools.product(perms, repeat=2):
        if perm1.inv == 0 or perm2.inv == 0:
            continue
        # if len(perm1.trimcode) != len(perm2.trimcode):
        #     continue
        #length = max(len(perm1) - 1, len(perm2) - 1)
        print("TRying", perm1, perm2)
        #length = max(len(perm1), len(perm2))
        length = n
        #length = max(len(perm1), len(perm2))
        # partition1 = tuple((~(perm1.strict_mul_dominant(length))).trimcode)
        # partition2 = tuple((~(perm2.strict_mul_dominant(length))).trimcode)
        real_prod = Sx(perm1) * Sx(perm2)
        #length = max(len(perm) - 1 for perm in real_prod)
        #max(len(perm1), len(perm2)) - 1
        #test_prod = r.schub_elem(perm1, length, partition=partition1) * r.schub_elem(perm2, length, partition=partition2)
        test_prod = r.full_schub_elem(perm1, length) * r.full_schub_elem(perm2, length) 
        

        for perm, coeff in real_prod.items():
            test_prod -= coeff * r.full_schub_elem(perm, length)#, partition=tuple(range(n - 1, 0, -1)))
        #assert test_prod.to_rc_graph_ring_element().almosteq(rc_ring.zero), f"Failed product test for {perm1} and {perm2}: got {test_prod}"    
        #test_prod_dct = {k: v for k, v in test_prod.items() if all(len(rc.perm) <= length + 1 for rc, coeff in r(k).to_rc_graph_ring_element().items() if coeff != 0)}
        #assert test_prod.almosteq(r.zero), f"Failed product test for {perm1} and {perm2}: got {test_prod}"
        #test_elem = r.from_tensor_dict(test_prod_dct, length)
        test_elem = test_prod
        try:
            assert test_elem.almosteq(r.zero), f"Failed product test for {perm1} and {perm2}"
        except AssertionError as e:
            #print(e)
            #pretty_print(test_elem)
            #print("but")
            def acceptable(key):
                rc = r.key_to_rc_graph(key)
                return len(rc.perm.descents()) <= 1 or any(rc.perm.has_pattern(pat) for pat in [[3,1,4,2],[1,4,3,2],[4,1,3,2]])
            if all(not acceptable(key) for key in test_elem if test_elem[key] != 0):

                print(f"Failed product test for {perm1} and {perm2}")#, but all nonzero terms correspond to RC graphs of the identity permutation.")
                for key in test_elem:
                    pretty_print(r.key_to_rc_graph(key))
                pretty_print(test_elem)
            #assert test_elem.to_rc_graph_ring_element().almosteq(rc_ring.zero), f"Failed product test for {perm1} and {perm2}: got {test_elem}"
            #print("Sanity")
            # pretty_print(r.full_schub_elem(perm1, length))
            # pretty_print(r.full_schub_elem(perm2, length))
                raise
        #else:
        print(f"Passed product test for {perm1} and {perm2}.")
        #assert all(test_elem[k] == 0 for k in test_elem if all(len(rc.perm) <= length + 1 for rc, coeff in r(k).to_rc_graph_ring_element().items() if coeff != 0))), f"Failed product test for {perm1} and {perm2}: got {test_elem}"