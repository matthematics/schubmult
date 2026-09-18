"""Test squash product coincides with Fomin-Kirillov"""

from schubmult import *

if __name__ == "__main__":
    from schubmult.symbolic import expand, S
    from schubmult.rings.combinatorial.dual_rc_graph_ring import DualRCGraphRing
    from schubmult.abc import x
    from schubmult.utils.perm_utils import weak_compositions
    from schubmult.visualization import draw_pipe_dream_tikz
    from sympy import pretty_print, pretty, Mul
    import itertools
    import sys

    n = int(sys.argv[1])
    #m = int(sys.argv[2])
    r = DualRCGraphRing()
    #comps = weak_compositions(n, m)
    perms = Permutation.all_permutations(n)
    prod_dict = {}
    for perm in perms:
        result = Sx(perm)
        rc_result = 0
        for rc in RCGraph.all_rc_graphs(perm, n):
            for p in range(n + 1):
                elem_rc_elem = r.schub(uncode([0] * (n - p) + [1] * (p + 1)), n)
                for elem_rc in elem_rc_elem:
                    rc_result += r(rc) * r(elem_rc)