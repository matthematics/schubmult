import math
import sys

from schubmult import ASx, Sx, uncode
from schubmult.rings import FA, WordBasis
from schubmult.rings.free_algebra_basis import SchubertBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, TensorModule
from schubmult.utils.perm_utils import artin_sequences


def product_too_big(perm1, perm2, n):
    return any(len(perm) > n for perm in (Sx(perm1)*Sx(perm2)).keys())

def main():
    n = int(sys.argv[1])

    unit_rc_module = RCGraphModule({RCGraph(): 1})

    # 100% positive!

    unit_tensor_rc_module = TensorModule.ext_multiply(unit_rc_module, unit_rc_module)

    solution_module = TensorModule()

    aseqs = artin_sequences(n - 1)

    # THIS IS THE CORRECT COPRODUCT
    for seq in aseqs:
        perm = uncode(seq)
        elem = ASx(perm, n - 1).change_basis(WordBasis)

        for word, coeff in elem.items():
            modmod = FA(*word).coproduct() * unit_tensor_rc_module
            for word_double, coeff_double in elem.items():
                mod = FA(*word_double) * unit_rc_module
                for rc, coeff2 in mod.items():
                    assert rc.length_vector() == word_double
                    solution_module += coeff * coeff_double * coeff2 * TensorModule.ext_multiply(1 * rc, modmod)

    coprods_interim = {}
    coprods = {}
    coprods_length_vector = {}
    coprods_length_vector2 = {}
    coprods_interim_rc = {}
    products = {}
    products_rc = {}
    weight_sum = {}
    for (rc, (rc1, rc2)), coeff in solution_module.items():
        perm = rc.perm
        perm1, perm2 = rc1.perm, rc2.perm
        if len(perm1) > n or len(perm2) > n or len(perm) > n:
            continue
        # assert coeff >= 0 NOPE
        # coprods_interim[rc.perm] = coprods_interim.get(rc.perm, TensorModule()) + coeff * TensorModule.ext_multiply(1*rc1, ASx(perm2, len(rc2)))
        # coprods_interim_rc[rc.perm] = coprods_interim.get(rc.perm, TensorModule()) + coeff * TensorModule.ext_multiply(Sx(rc1.perm), 1 * rc2)
        # coprods_length_vector[rc.length_vector()] = coprods_length_vector.get(rc.length_vector(), 0) + coeff * (FA@FA)((rc1.length_vector(), rc2.length_vector()))
        # coprods_length_vector2[rc.length_vector()] = coprods_length_vector2.get(rc.length_vector(), 0) + coeff * (ASx@ASx)(((rc1.perm,len(rc1)), (rc2.perm,len(rc2))))
        products[(rc1.perm, rc2.perm)] = products.get((rc1.perm, rc2.perm), RCGraphModule()) + coeff * RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(rc.perm, n - 1), 1))

    for (perm1, perm2), module in products.items():
        for rc, coeff in module.items():
            weight_sum[rc.length_vector()] = weight_sum.get(rc.length_vector(), 0) + coeff
        #products_rc[(rc1, rc2)] = products_rc.get((rc1, rc2), 0) + coeff * Sx(perm)

    assert all(v > 0 for v in weight_sum.values()), f"{weight_sum}" # total sums all >= 0
    

    # for perm, module in coprods_interim.items():
    #     #assert all(v >=0 for v in module.values()), "NOPE"
    #     for (rc, key), coeff in module.items():
    #         coprods[perm] = coprods.get(perm, 0) + coeff * (ASx@ASx)(((rc.perm, len(rc)),key))


    # for perm, module in coprods_interim_rc.items():
    #     assert all(v >=0 for v in module.values()), "NOPE"
    #     for (perm1, rc2), coeff in module.items():
    #         coprods[perm] = coprods.get(perm, 0) + coeff * (ASx@ASx)(((perm1, n-1),(rc2.perm, len(rc2))))

    # for length_vector, sx_elem in coprods_length_vector2.items():
    #     print(sx_elem)
    #     #assert all(v >= 0 for v in fa_elem.values()), "NOPE"
    #     elem = FA(*length_vector).change_basis(SchubertBasis)
    #     for key, v in elem.items():
    #         assert v * coeff >= 0, "NOPE"
    #         coprods[key] = coprods.get(key, 0) + coeff * v * sx_elem

    # for length_vector, fa_elem in coprods_length_vector.items():
    #     print(fa_elem)
    #     #assert all(v >= 0 for v in fa_elem.values()), "NOPE"
    #     elem = FA(*length_vector).change_basis(SchubertBasis)
    #     for key, v in elem.items():
    #         assert v * coeff >= 0, "NOPE"
    #         coprods[key] = coprods.get(key, 0) + coeff * v * 1

    num_successes = 0

    for perm, elem in coprods.items():
        print(f"{perm.trimcode}")
        print(elem)
        check = ASx(perm, n - 1).coproduct()
        print(check)
        diff = elem - check
        print(diff)
        assert all(v == 0 for v in diff.values()), f"Failed check on {perm}, diff = {diff}"
        num_successes += 1


    for (perm1, perm2), elem in products.items():
        if product_too_big(perm1, perm2, n):
            continue
        print(f"{perm1.trimcode}, {perm2.trimcode}")
        check = Sx(perm1)*Sx(perm2)
        sumup = 0
        for rc, coeff in elem.items():
            # diff = elem - check
            # print(diff)
            assert check.get(rc.perm, 0) == coeff
            if rc.is_principal:
                sumup += coeff * Sx(rc.perm)
        assert sumup == check
        print(f"Success {sumup}")
        num_successes += 1

    # for (rc1, rc2), elem in products_rc.items():
    #     perm1, perm2 = rc1.perm, rc2.perm
    #     if product_too_big(perm1, perm2, n):
    #         continue
    #     print(f"{perm1.trimcode}, {perm2.trimcode}")
    #     print(elem)
    #     check = Sx(perm1)*Sx(perm2)
    #     print(check)
    #     # diff = elem - check
    #     # print(diff)
    #     assert elem == check
    #     num_successes += 1

    

    #assert num_successes == math.factorial(n), f"Only {num_successes} successes, not the full group"

if __name__ == "__main__":
    main()
