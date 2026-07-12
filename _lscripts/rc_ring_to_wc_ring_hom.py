from schubmult import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    w = WCGraphRing()
    r = RCGraphRing()
    perms = Permutation.all_permutations(n)
    def morphism(wc_elem):
        return r.from_dict({RCGraph(wc): coeff for wc, coeff in wc_elem.items() if wc.is_reduced})
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1), RCGraph.all_rc_graphs(perm2)):
            for len1, len2 in itertools.product(range(len(rc1), n+1), range(len(rc2), n+1)):
                wc1 = w(WCGraph(rc1.resize(len1)))
                wc2 = w(WCGraph(rc2.resize(len2)))
                rc_prod = r(rc1.resize(len1)) * r(rc2.resize(len2))
                wc_prod = wc1 * wc2
                rc_try_prod = morphism(wc_prod)
                if not rc_prod.almosteq(rc_try_prod):
                    print(f"Mismatch for {rc1}, {rc2}: {rc_prod} != {rc_try_prod}")
        print(f"Done with perm1={perm1}, perm2={perm2}")
        