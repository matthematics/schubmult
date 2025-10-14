if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.rc_graph_ring import RCGraphRing
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
        pretty_print(elem)
        # elem = ring.from_free_algebra_element(ASx(perm, len(perm.trimcode)))
        # pretty_print(elem)
    