from sympy import init_printing, pretty_print




if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    from schubmult.abc import x
    from symengine import S
    from sympy import pretty_print

    from schubmult import Permutation, RCGraph, RCGraphRing, RootTableau
    # from schubmultutils.perm_utils import artin_sequences

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    # sKEW DIV DIFF WEIGHT
    # is dual pieri Cauchy?
    dom_perms = {perm.minimal_dominant_above() for perm in perms}
    for w in perms:
        print(f"{w.trimcode}")
        for rc in RCGraph.all_rc_graphs(w, n-1):
            rt = RootTableau.from_rc_graph(rc)
            pretty_print(rt)