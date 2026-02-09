from schubmult import *
from sympy import pretty_print
import numpy as np

def vexgen(rc):
    r = RCGraphRing()
    if rc.perm.is_vexillary:
        return rc
    downrow = len(rc) - 1
    while downrow > 0 and not rc.rowrange(downrow).perm.is_vexillary:
        downrow -= 1
    downrow += 1
    cut0, cut1 = rc.vertical_cut(downrow)


if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm):
            if rc.perm.is_vexillary:
                continue
            vex_rc = rc
            while not vex_rc.perm.is_vexillary:
                if len(vex_rc[-1]) == 0:
                    if len(vex_rc.perm.trimcode) < len(vex_rc):
                        vex_rc = RCGraph(vex_rc.shiftup(len(vex_rc) - len(vex_rc.perm.trimcode)))
                    vex_rc = vex_rc.zero_out_last_row().resize(len(rc))
                else:
                    vex_rc = vex_rc.resize(len(rc) + 1)
                    if len(vex_rc.perm.trimcode) < len(vex_rc):
                        vex_rc = RCGraph(vex_rc.shiftup(len(vex_rc) - len(vex_rc.perm.trimcode))) 
                    vex_rc = vex_rc.zero_out_last_row().resize(len(rc))
            assert vex_rc.perm.is_vexillary, f"Failed for {perm} with vexillary RC {vex_rc}"
            print("RC")
            pretty_print(rc)
            print("Vexillary RC")
            pretty_print(vex_rc)