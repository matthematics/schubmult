from schubmult import *
from sympy import pretty_print, S

r = RCGraphRing()

def bpd_auto(rc):
    if isinstance(rc, RCGraphRingElement):
        ret = r.zero
        for rc2, coeff in rc.items():
            ret += coeff * r(bpd_auto(rc2))
        return ret
    if len(rc) == 0:
        return RCGraph()
    return BPD.from_rc_graph(rc.transpose()).transpose().to_rc_graph().resize(len(rc))


if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    for perm1, perm2 in itertools.product(perms, repeat=2):
        # if perm.is_dominant:
        #     continue
        if perm2.inv != 1 or perm1.inv != 1:
            continue
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n - 1), RCGraph.all_rc_graphs(perm2,n - 1)):
            

            try:
                r(rc1) % r(rc2)
            except Exception as e:
                continue
            print(f"Testing bob")
            pretty_print(rc1)
            pretty_print(rc2)
            pretty_print(r(rc1)%r(rc2))
            print("Testing joe")
            auto_rc1 = bpd_auto(rc1)
            auto_rc2 = bpd_auto(rc2)
            pretty_print(r(auto_rc2) % r(auto_rc1))
            assert (r(auto_rc2) % r(auto_rc1)).almosteq(bpd_auto(r(rc1) % r(rc2)))
