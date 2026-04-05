from schubmult import *

def transition_rc(rc):
    ring = RCGraphRing()
    if rc.perm.inv == 0:
        return (rc,)
    if len(rc[-1]) == 0:
        return (*transition_rc(rc.zero_out_last_row()), rc)
    lower_rc = rc.extend(1).exchange_property(len(rc)).zero_out_last_row()
    return (*transition_rc(lower_rc), rc) 


if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1), RCGraph.all_rc_graphs(perm2)):
            cp1 = r(rc1).coproduct()
            cp2 = r(rc2).coproduct()
            prd_cp = cp1 * cp2
            cp_prd = (r(rc1) * r(rc2)).coproduct()
            pretty_print(prd_cp)
            pretty_print(cp_prd)
            assert prd_cp.almosteq(cp_prd), f"Failed for {rc1}, {rc2}: {prd_cp-cp_prd}"
            print("yay")