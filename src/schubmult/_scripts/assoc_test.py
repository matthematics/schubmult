from sympy import init_printing, pretty_print




if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    from schubmult.abc import x
    from symengine import S
    from sympy import pretty_print

    from schubmult import Permutation, RCGraph, RCGraphRing
    # from schubmultutils.perm_utils import artin_sequences

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    # sKEW DIV DIFF WEIGHT
    # is dual pieri Cauchy?
    dom_perms = {perm.minimal_dominant_above() for perm in perms}
    for w in perms:
        for u in dom_perms:
            if not u.bruhat_leq(w):
                continue
            for v in perms:
                for rc in RCGraph.all_rc_graphs(v, n - 1):
                    for (rc1, rc2) in rc.dualpieri(u, w):
                        print(f"we gots the {v.trimcode} rc")

                        pretty_print(rc)
                        print(f"We iz gots the results of acting with da {w.trimcode}/{u.trimcode} dual pieri:")
                        pretty_print((rc1, rc2))