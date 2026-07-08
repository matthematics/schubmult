from schubmult import *

r = RCGraphRing()
br = BoundedRCFactorAlgebra()

def ring_elem_to_schubert(ring_elem):
    return Sx.from_dict({rc.perm: coeff for rc, coeff in ring_elem.items() if rc.is_principal})

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm1, perm2 in itertools.product(perms, repeat=2):
        real_product = (br.full_schub_elem(perm1, n) * br.full_schub_elem(perm2, n)).to_rc_graph_ring_element()
        the_product = Sx(perm1) * Sx(perm2)
        assert ring_elem_to_schubert(real_product).almosteq(the_product), f"Mismatch for {perm1=}, {perm2=}: {real_product=} vs {the_product=}"
        try_product = (br.from_rc_graph_ring_element(r.schub(perm1, n), n) * br.from_rc_graph_ring_element(r.schub(perm2, n), n))
        if  not ring_elem_to_schubert(try_product.to_rc_graph_ring_element()).almosteq(the_product):
            print("The real test")
            test_product = try_product.straighten().to_rc_graph_ring_element()
            assert ring_elem_to_schubert(test_product).almosteq(the_product), f"Mismatch for {perm1=}, {perm2=}: {test_product=} vs {the_product=}"
            print("Yay")
