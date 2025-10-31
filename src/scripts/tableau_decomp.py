

if __name__ == "__main__":
    import sys
    from schubmult.schub_lib.rc_graph import *
    from schubmult import *
    from sympy import pretty_print
    perms = Permutation.all_permutations(int(sys.argv[1]))
    for perm in perms:
        print(f"{perm.trimcode} tableaux:")
        for rc in RCGraph.all_rc_graphs(perm):
            pretty_print(rc.tableau_decomp())