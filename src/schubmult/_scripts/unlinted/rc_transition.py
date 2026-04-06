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
    from schubmult.rings.combinatorial.hw_rc_ring import HWRCGraphRing
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = HWRCGraphRing()
    rr = RCGraphRing()
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for rc1, rc2 in itertools.product(RCGraph.all_hw_rcs(perm1, n - 1), RCGraph.all_hw_rcs(perm2, n - 1)):
            cp1 = r(rc1).coproduct()
            cp2 = r(rc2).coproduct()
            prd_cp_low = cp1 * cp2
            prd_cp = (r@r).zero
            for (key1, key2), coeff in prd_cp_low.items():
                prd_cp += coeff * r._snap_highest_weight(r(key1)) @ r._snap_highest_weight(r(key2))
            cp_prd_low = (r(rc1) * r(rc2)).coproduct()
            cp_prd = (r@r).zero
            for (key1, key2), coeff in cp_prd_low.items():
                cp_prd += coeff * r._snap_highest_weight(r(key1)) @ r._snap_highest_weight(r(key2))
            pretty_print(prd_cp)
            pretty_print(cp_prd)
            assert prd_cp.almosteq(cp_prd), f"Failed for {rc1}, {rc2}: {prd_cp-cp_prd}"
            print("yay")