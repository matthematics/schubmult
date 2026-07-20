from schubmult import *
from schubmult.abc import *

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        if perm.max_descent != 2:
            continue    
        for p in range(1, perm.max_descent):
            try_sum = 0
            real_sum = 0
            for wc in WCGraph.all_wc_graphs(perm, perm.max_descent):
                real_sum += wc.polyvalue(x, y, beta=1, minus_convention=True)
                wc_up, wc_down = wc.vertical_cut(p)
                cutt = wc.resize(p).resize(len(wc))
                left_sum = 0
                working_rc = wc_up
                for i in range(len(wc_up[0]) -1, -1, -1):
                    (r1, c1) = wc_up.left_to_right_inversion_coords(i)
                    (r2, c2) = cutt.left_to_right_inversion_coords(i)
                    the_diff = (y[c1] + y[c1])*(1 + y[c2])**(-1)
                    working_rc = working_rc.toggle_ref_at(r2, c2)
                    left_sum += the_diff * working_rc.polyvalue(x, y, beta=1, minus_convention=True)
                    working_rc = working_rc.toggle_ref_at(r1, c1)    
                try_sum += left_sum * wc_down.polyvalue(x[p:], y, beta=1, minus_convention=True)
            diff = (try_sum - real_sum).expand()
            try:
                assert diff == 0, f"WCGraph for {perm} does not match: {diff=}, \n{try_sum=}\n {real_sum=}"
            except AssertionError as e:
                print(e)
                continue
            print(f"WCGraph for {perm} verified for p={p}.")