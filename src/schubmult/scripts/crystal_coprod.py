if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.rc_graph_ring import RCGraphRing, tensor_to_highest_weight
    from sympy import pretty_print
    import sys
    ring = RCGraphRing()
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    tring = ring @ ring
    for perm in perms:
        cprd = ASx(perm).change_basis(WordBasis).coproduct()
        elem = tring.zero
        for (w1, w2), coeff in cprd.items():
            elem += coeff*tring.ext_multiply(ring.from_free_algebra_element(FA(*w1)),ring.from_free_algebra_element(FA(*w2)))
        elem = tensor_to_highest_weight(elem)
        elem = elem.ring.from_dict({(k1, k2): v for (k1, k2), v in elem.items() if k1.perm.bruhat_leq(perm) and k2.perm.bruhat_leq(perm)})
        print(f"Coprod of {perm.trimcode=}")
        pretty_print(elem)
        # elem = ring.from_free_algebra_element(ASx(perm, len(perm.trimcode)))
        # pretty_print(elem)
    