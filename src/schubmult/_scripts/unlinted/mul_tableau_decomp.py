from schubmult import *
from schubmult.symbolic import *

if __name__ == "__main__":
    import sys
    from sympy import pretty_print

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    for perm in perms:
        if len(perm.descents()) <= 1:
            continue
        #for rc in RCGraph.all_lw_rcs(perm, len(perm.trimcode)):
        if True:
            rc = RCGraph.principal_rc(perm)
            print("RC:")
            pretty_print(rc)
            elem = prod([r(rcc) for rcc in rc.tableau_decomp()])
            #pretty_print(elem)
            counts = {}
            for rcc, coeff in elem.items():
                if not rcc.is_lowest_weight:
                    continue
                counts[len(rcc.perm.descents())] = counts.get(len(rcc.perm.descents()), 0) + 1
            print("Descents: ", len(rc.perm.descents()))
            print("Num descents in product:")
            for k in sorted(counts.keys()):
                print(f"  {k}: {counts[k]}")
            # try:
            #     assert rc == max(elem.keys(), key=lambda x: tuple(x.perm.trimcode)), f"Failed for RC graph of {perm}"
            # except AssertionError as e:
            #     print(e)
            print(elem)
            #    sys.exit(1)