import sys

from schubmult import Permutation, RCGraphRing, CrystalGraphTensor
from schubmult.rings.combinatorial.bounded_rc_factor_algebra import BoundedRCFactorAlgebra

if __name__ == "__main__":
    n = int(sys.argv[1])

    g = BoundedRCFactorAlgebra()
    perms = Permutation.all_permutations(n)
    length = n + 2
    for perm in perms:
        if perm.inv == 0:
            continue
        elem = g.schub_elem(perm, length)
        # for key, coeff in elem.items():
        #     if coeff == 0:
        #         continue
        #     single_term = g.from_dict({key: 1})
            #test_elem = single_term
        used = set()
        for key, coeff in elem.items():
            test_elem = g.zero
            hw_key = CrystalGraphTensor(*[key[i].resize(length) for i in range(len(key))]).to_highest_weight()[0]
            if hw_key in used:
                continue
            used.add(hw_key)
            for the_key in hw_key.truncated_crystal(length=length):
                test_elem += coeff * g(g.make_key(CrystalGraphTensor(*[the_key[i].normalize() for i in range(len(the_key))]), length))
            for descent in perm.descents():
                test1 = test_elem.divdiff_descs(descent + 1).to_rc_graph_ring_element()
                test2 = test_elem.to_rc_graph_ring_element().divdiff(descent + 1)
                if not test1.almosteq(test2):
                    diff = test1 - test2
                    if any(rc.perm == perm.swap(descent, descent + 1) for rc, v in diff.items() if v != 0):
                        print(f"FAIL: perm={perm}, elem={test_elem}, descent={descent}")
                        print(f"  divdiff_descs then to_rc: {test1}")
                        print(f"  to_rc then divdiff:       {test2}")
                        print(f"{next(iter(test_elem.to_rc_graph_ring_element())).perm=}")
                        print(f"{test1 - test2=}")
                        sys.exit(1)
            print(f"PASS: perm={perm}")

    print(f"All permutations passed for n={n}")
