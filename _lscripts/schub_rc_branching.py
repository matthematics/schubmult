from schubmult import *
from schubmult.mult.double import schubmult_double

def branch_that_double_schub(perm, genset, coeff_genset, p):
    if p <= 0:
        raise ValueError("p must be a positive integer")
    if p > perm.max_descent:
        raise ValueError("p must be less than or equal to the maximum descent of the permutation")
    result = 0
    n = len(perm)
    w0up = (Permutation.w0(n) * Permutation.w0(n - p))
    w0down = Permutation.w0(n - p)
    w0 = Permutation.w0(n)
    for rc in RCGraph.all_rc_graphs(perm, n - 1):
        rc_up, rc_down = rc.vertical_cut(p)
        #coeff = rc.resize(p).resize(len(rc)).polyvalue([-yy for yy in coeff_genset],[-yy for yy in coeff_genset])
        #schubmult_double({rc_up.perm * w0up: 1}, rc_down.perm * w0down,[-yy for yy in coeff_genset[:20]],[-yy for yy in coeff_genset[:20]]).get(perm * w0, 0)
        coeff = 1
        result += coeff * rc_up.polyvalue(genset, coeff_genset) * rc_down.polyvalue(genset[p:], coeff_genset)
    return result

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    ring = DSx([]).ring
    for perm in perms:
        real_result = DSx(perm).as_polynomial().expand()
        for p in range(1, perm.max_descent + 1):
            branch_result = branch_that_double_schub(perm, ring.genset, ring.coeff_genset, p).expand()
            if (real_result - branch_result).expand() != 0:
                print(f"Discrepancy for permutation {perm} at p={p}:")
                print(f"Real result: {real_result}")
                print(f"Branch result: {branch_result}")
                print(f"Difference: {(real_result - branch_result).expand()}")
            else:
                print(f"Permutation {perm} at p={p} matches.")