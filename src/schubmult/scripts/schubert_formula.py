import sys

from schubmult import ASx, Permutation, uncode
from schubmult.rings import FA, Sx, WordBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, TensorModule


def vector_sum(v1, v2):
    return tuple([a + b for a, b in zip(v1, v2)])


def main():
    n = int(sys.argv[1])

    unit_rc_module = RCGraphModule({RCGraph(): 1})

    # 100% positive!
    # aseqs = artin_sequences(n - 1)
    # perm_modules = {}
    # ring = ASx @ ASx

    def asxt(rc):
        return (rc.perm, len(rc))

    unit_tensor_rc_module = TensorModule.ext_multiply(unit_rc_module, unit_rc_module)

    perm_modules3 = {}

    # this is an inner product of RCs
    perms = Permutation.all_permutations(n)
    # additionally extract the coefficient of a specific rc graph by length vector
    for perm in perms:
        perm_words = ASx(perm, n - 1).change_basis(WordBasis)
        # mod = ASx(perm, n - 1) * unit_rc_module
        alternative_magic_coeffs = {}
        for rc4 in RCGraph.all_rc_graphs(perm, n - 1):
            extra_coeff = perm_words.get(rc4.length_vector(), 0)
            magic_coeffs = {}

            # PUT THIS BACK IF NEEDED
            # multiplier = 0
            # for rc3, coeff3 in mod.items(): #mod is ASx rc graph module, coeff3 ia an ASx perm coeff
            #     if len(rc3.perm) > n:
            #         continue
            #     if rc3.perm != rc4.perm:
            #         continue
            #     multiplier += coeff3
            # assert multiplier == 1
            # this inner products it
            # sumtiplier = 0
            # RC graph bijection

            # the coproduct of the ASx RC graphs
            # selects a term, termwise
            for (
                aseq,
                coeff,
            ) in perm_words.items():  # ASx perm word coeff is coeff
                mod2 = FA(*aseq).coproduct() * unit_tensor_rc_module
                # coeff is perm_words for some coeff3
                for (rc1, rc2), coeff1 in mod2.items():  # positive coproduct coeff is coeff1
                    if len(rc1.perm) > n or len(rc2.perm) > n:
                        continue
                    assert coeff1 == 1
                    perm1, perm2 = rc1.perm, rc2.perm
                    magic_coeffs[(perm1, perm2)] = magic_coeffs.get((perm1, perm2), 0) + coeff * coeff1
                    alternative_magic_coeffs[(perm1, perm2, perm, rc4.length_vector())] = (
                        alternative_magic_coeffs.get((perm1, perm2, perm, rc4.length_vector()), 0)
                        + coeff
                        * coeff1
                        * extra_coeff
                    )
                    # print(f"{perm1.trimcode} {perm2.trimcode} += {multiplier*coeff*coeff1}")
                    # print(RCGraphTensor(rc1, rc2))
                    # print("inner product")
                    # print(rc4)
                    # sumtiplier += multiplier * coeff
                # this specifically pulls out a single coefficient
            for key, magic_coeff in magic_coeffs.items():
                assert magic_coeff >= 0
                #assert alternative_magic_coeffs.get((*key,rc4.length_vector()), 0) == magic_coeff * 
                perm_modules3[key] = perm_modules3.get(key, RCGraphModule()) + magic_coeff * rc4

    for (perm1, perm2, perm, seq), val in alternative_magic_coeffs.items():
        assert val == sum([v for rc, v in perm_modules3[(perm1, perm2)].items() if rc.length_vector() == seq and rc.perm == perm])

    for (perm1, perm2), elem in perm_modules3.items():
        product = Sx(perm1) * Sx(perm2)
        if any(len(perm) > n for perm in product.keys()):
            continue
        print(perm1.trimcode, perm2.trimcode)
        print(elem)
        val = Sx([]).ring.zero
        for rc, coeff in elem.items():
            assert coeff == product.get(rc.perm, 0)
            val[rc.perm] = coeff
        assert val == product
        print(f"Success {perm1.trimcode} {perm2.trimcode}")
    exit()


if __name__ == "__main__":
    main()
