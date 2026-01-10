from schubmult import *

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    ring = RCGraphRing()
    for perm in perms:
        rc = RCGraph.principal_rc(perm, len(perm.trimcode))
        cp = ring.coproduct_on_basis(rc)
        assert all(rc.is_potential_coproduct(rc1, rc2) for (rc1, rc2), coeff in cp.items() if coeff != 0), f"Failed on {rc}"
        print("Paint")