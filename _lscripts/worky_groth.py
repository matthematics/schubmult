from schubmult import *
from schubmult.combinatorics.pipe_dream import PipeDream
from schubmult.mult.groth_double import *
import numpy as np

if __name__ == "__main__":
    import sys
    from schubmult.abc import y
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    poly_results = {perm: DSx([]).ring.zero for perm in perms}
    w0 = Permutation.w0(n)
    # for perm in perms:
    #     for rc in RCGraph.all_rc_graphs(perm * w0, n):
    #         cpd = PipeDream.from_rc_graph(rc).co_pipe_dream()
    #         keys = tuple(poly_results.keys())
    #         cperm = cpd.perm * w0
    #         result = 1
    #         for (_, col) in [rc.left_to_right_inversion_coords(i) for i in range(w0.inv - perm.inv)]:
    #             result *= (1 + Gx._beta * y[col])
    #         if perm.bruhat_leq(cperm):
    #             #crosses = np.argwhere(cpd._grid == cpd.CROSS)
    #             poly_results[cperm] += (Gx._beta**(cperm.inv - perm.inv)) * result * DSx(perm)
    for perm in perms:#, poly in poly_results.items():
        poly = DSx([]).ring.from_dict(dgroth_to_dschub_positive(perm, y, Gx._beta))
        assert (DGx(perm).as_polynomial() - poly.as_polynomial()).expand() == 0, f"Mismatch for permutation {perm}, expected {DGx(perm).as_polynomial()}, got {poly.as_polynomial()}"
        print(f"Funky spinach {perm}")