
if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.crystal_graph import CrystalGraph, CrystalGraphTensor
    from schubmult.utils.perm_utils import elem_sym_perms
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n-1):
            for k in range(1, n):
                for p in range(1, k + 1):
                    result = None
                    # for cut in range(p, len(rc) + 1):
                    #     result = rc.vertical_cut(cut)[0].monk_crystal_mul(p, min(cut,k), prev_result=None)
                    result = rc.monk_crystal_mul(p, k, warn=False)
                    print(f"Success {p=} {k=}")
                    pretty_print(rc)
                    print("Result:")
                    pretty_print(result)
                # if len(results) != 1:
                    #     print("AMBIGUOUS")
                    #     raise AssertionError
                # result = rc.flat_elem_sym_mul(k)
                # if result:
                #     print(f"Success {k=}")
                #     pretty_print(rc)
                #     print("Result:")
                #     pretty_print(result)
                # assert result is not None
                