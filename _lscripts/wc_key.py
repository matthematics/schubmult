import sys
from collections import defaultdict

from schubmult import IncreasingTableau, Permutation, WCGraph, x
from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra

KeyPoly = PolynomialAlgebra(KeyPolyBasis(x))

n = int(sys.argv[1])

failures = []
for w in Permutation.all_permutations(n):
    all_wc = list(WCGraph.all_wc_graphs(w, n - 1))
    classes = defaultdict(list)

    for wc in all_wc:
        ptab, _ = IncreasingTableau.hecke_column_insert_rsk(wc.compatible_sequence, wc.perm_word)
        ptab_key = tuple(tuple(row) for row in ptab._grid.tolist())
        classes[ptab_key].append(wc)

    for ptab, gs in classes.items():
        poly = sum(
            (g.polyvalue(x, beta=1, prop_beta=True) for g in gs),
            0,
        )
        key_exp = KeyPoly.from_expr(poly)
        nonzero = [(kk, vv) for kk, vv in key_exp.items() if vv != 0]

        print("perm =", w)
        print("P-tableau =", ptab)
        print("sum poly =", poly)
        print("key expansion =", nonzero)
        print()

        if any(vv < 0 for _, vv in nonzero):
            failures.append((w, ptab, len(gs), nonzero))

print("failures =", len(failures))
for failure in failures:
    print("FAIL", failure)