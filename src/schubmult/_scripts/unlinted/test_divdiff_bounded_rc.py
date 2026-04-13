import sys

from schubmult import Permutation, RCGraphRing
from schubmult.rings.combinatorial.bounded_rc_factor_algebra import BoundedRCFactorAlgebra

if __name__ == "__main__":
    n = int(sys.argv[1])

    g = BoundedRCFactorAlgebra()
    perms = Permutation.all_permutations(n)

    for perm in perms:
        elem = g.schub_elem(perm, n + 2)
        # for key, coeff in elem.items():
        #     if coeff == 0:
        #         continue
        #     single_term = g.from_dict({key: coeff})
        for descent in perm.descents():
            test1 = elem.divdiff_descs(descent).to_rc_graph_ring_element()
            test2 = elem.to_rc_graph_ring_element().divdiff(descent)
            if not test1.almosteq(test2):
                print(f"FAIL: perm={perm}, elem={elem}, descent={descent}")
                print(f"  divdiff_descs then to_rc: {test1}")
                print(f"  to_rc then divdiff:       {test2}")
                print(f"{elem=}")
                sys.exit(1)
        print(f"PASS: perm={perm}")

    print(f"All permutations passed for n={n}")
