from schubmult import *
from schubmult.abc import e, x
from schubmult.symbolic import prod
from schubmult.utils.schub_lib import kdown_perms

def elem_monom_formula(n):
    """???"""
    perms = Permutation.all_permutations(n)
    # seen = {}
    #build_polys = {}
    w0 = Permutation.w0(n)
    for perm in perms:
        # rc_t_by_weight = {}
        # for rc_t in RCGraph.all_rc_graphs(~perm, n - 1):
        #     # if rc.extremal_weight != RCGraph.principal_rc(perm, n - 1).length_vector:
        #     #     continue
        #     wt = rc_t.length_vector
        #     rc_t_by_weight[wt] = rc_t_by_weight.get(wt, 0) + 1
        buildup = 0
        #for wt, ct in rc_t_by_weight.items():
        n = len(perm)
        w0 = Permutation.w0(n)
        stack = [(1, 1, [], perm * w0)]
        while stack:
            index, sgn, path, permu = stack.pop()
            if index == n and permu.inv == 0:
                buildup += Sx([]) * (sgn * prod([e(n - i - path[i - 1], n - i, x[1:]) for i in range(1, n) if path[i - 1]<=n-i]))
                continue
            monoperm = Permutation.w0(2 * n - index)
            for new_perm, diff, new_sgn in kdown_perms(permu, monoperm, n - index, index):
                # if n - index - diff < 0:
                #     continue
                if len(new_perm) <= len(permu):
                    stack.append((index + 1, sgn * new_sgn, path + [diff], new_perm))
        if buildup != Sx(perm):
            print(f"buildup {buildup} != Sx({perm})")
        else:
            print("Successfully verified", perm)

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    elem_monom_formula(n)