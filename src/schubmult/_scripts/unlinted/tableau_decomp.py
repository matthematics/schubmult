

if __name__ == "__main__":
    import sys

    from schubmult import *
    from sympy import pretty_print

    from schubmult import *
    perms = Permutation.all_permutations(int(sys.argv[1]))
    for perm in perms:
        print(f"{perm.trimcode} tableaux:")
        for rc in RCGraph.all_rc_graphs(perm):
            pretty_print(rc.tableau_decomp())
