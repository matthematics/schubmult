from schubmult import *
from schubmult.symbolic import S, expand

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    ring = DGx([]).ring
    for perm in perms:
        if perm.inv == 0:
            continue
        persch = ring._as_schub(DGx(Permutation.w0(n))).divdiff_perm((~perm)*Permutation.w0(n)).as_polynomial()
        persch2 = S.Zero
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            the_prod = S.One
            for row in range(1, n):
                for col in range(1, n + 1 - row):
                    if rc.has_element(row, col):
                        the_prod *= (ring.genset[row] + ring.coeff_genset[col] + ring._beta * ring.genset[row] * ring.coeff_genset[col])
                    else:
                        the_prod *= (1 + ring._beta * ring.coeff_genset[col])

            persch2 += the_prod
        assert (persch - persch2).expand() == 0, f"Failed for {perm}: {(persch-persch2).expand()}"
        print("Hot poatto")
