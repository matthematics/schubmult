if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.rc_graph_ring import RCGraphRing
    from sympy import pretty_print
    import sys
    ring = RCGraphRing()
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        elem = ring.from_free_algebra_element(ASx(perm, len(perm.trimcode)))
        print(elem)
        pretty_print(elem)
    