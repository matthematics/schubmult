if __name__ == "__main__":
    import sys

    from sympy import pretty_print

    from schubmult import Permutation, RCGraph, RCGraphRing

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    rc_ring = RCGraphRing()
    addup = {}
    for perm in perms:

        for rc in RCGraph.all_rc_graphs(perm, n-1):

            for d in rc.perm.descents():
                if len(rc[d]) > 0 and rc[d][-1] == d + 1 and rc.right_root_at(d + 1, 1) == (d+1, d+2):
                    rcd = rc.exchange_property(d + 1)
                    addup[(perm, d+1)]  = addup.get((perm,d+1),rc_ring.zero) + rc_ring(rcd)

    for (perm, desc), val in addup.items():
        print(f"{perm.trimcode=} {desc=}")
        pretty_print(val)
