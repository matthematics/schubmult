if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.crystal_graph import CrystalGraph, CrystalGraphTensor
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms}

    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n-1):
            for k in range(1, n):
                try:
                    result = rc.flat_elem_sym_mul(k)
                    print(f"Result for {perm} and k={k}:")
                    pretty_print(rc)
                    print(f"Multiplied by k={k}:")
                    pretty_print(result)
                except ValueError:
                    print(f"No result for {perm} and k={k}")