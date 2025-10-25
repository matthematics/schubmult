if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.crystal_graph import CrystalGraph, CrystalGraphTensor
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    dominant_graphs = {RCGraph.principal_rc(Permutation([]).swap(i, i+1),n-1) for i in range(n - 1)}

    for perm in perms:
        hw_rcs = set()
        for rc0 in RCGraph.all_rc_graphs(perm, n-1):
             hw_rcs.add(rc0.to_highest_weight()[0])
        for rc in hw_rcs:
            for dom in dominant_graphs:
                tens = CrystalGraphTensor(rc, dom)
                hw_set = tens.all_highest_weights()
                print("We have")
                pretty_print(rc.full_crystal)
                print("Descent")
                pretty_print(dom)
                print("Highest weights:")
                for hw in hw_set:
                    pretty_print(hw)
                print("=============")