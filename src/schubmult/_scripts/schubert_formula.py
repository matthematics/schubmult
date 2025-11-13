import sys

from schubmult import x
from schubmult import FA, Sx, WordBasis
from schubmult import RCGraph, RCGraphModule, TensorModule

from schubmult import ASx, Permutation


def vector_sum(v1, v2):
    return tuple([a + b for a, b in zip(v1, v2)])


def main():
    n = int(sys.argv[1])

    unit_rc_module = RCGraphModule({RCGraph(): 1})

    # 100% positive!

    unit_tensor_rc_module = TensorModule.ext_multiply(unit_rc_module, unit_rc_module)

    by_coefficient_dict = {}
    by_schub_dict = {}
    by_term_dict = {}

    # this is an inner product of RCs
    perms = Permutation.all_permutations(n)

    for perm in perms:
        perm_words = ASx(perm, n - 1).change_basis(WordBasis)
        for rc4 in RCGraph.all_rc_graphs(perm, n - 1):
            extra_coeff = perm_words.get(rc4.length_vector(), 0)
            alternative_magic_coeffs = {}
            term_coeffs = {}
            magic_coeffs = {}
            for (
                aseq,
                coeff,
            ) in perm_words.items():  # ASx perm word coeff is coeff
                mod2 = FA(*aseq).coproduct() * unit_tensor_rc_module
                do_sum = False
                if rc4.length_vector() == aseq:
                    do_sum = True
                for (rc1, rc2), coeff1 in mod2.items():  # positive coproduct coeff is coeff1
                    if len(rc1.perm) > n or len(rc2.perm) > n:
                        continue
                    assert coeff1 == 1
                    perm1, perm2 = rc1.perm, rc2.perm

                    # this works by taking term by term inner products with the RCs in the negative terms
                    # there's an extra inner sum with FA(*seq) * unit_rc_module times the coeff
                    # rc's from the pos mod, <ASx(rc4),ASx(rc3)> = delta perm(rc3), perm(rc4)
                    # rc's with the same permutation, sign of coeff, take sum of schub(rc3) adds up to schub(perm)
                    magic_coeffs[(perm1, perm2)] = magic_coeffs.get((perm1, perm2), 0) + coeff * coeff1

                    # This is the same as above for the principal one but 0 otherwise
                    alternative_magic_coeffs[(perm1, perm2)] = alternative_magic_coeffs.get((perm1, perm2), 0) + coeff * coeff1 * extra_coeff

                    if do_sum:
                        # this picks out the specific coefficient of the term
                        # we get the right weight but from different Schubs
                        assert coeff == extra_coeff
                        term_coeffs[(perm1, perm2)] = term_coeffs.get((perm1, perm2), 0) + coeff * coeff1 * extra_coeff

            # term by term and schub by schub
            for key, magic_coeff in magic_coeffs.items():
                assert magic_coeff >= 0
                by_coefficient_dict[key] = by_coefficient_dict.get(key, RCGraphModule()) + magic_coeff * rc4

            # schub by schub
            for key, amagic_coeff in alternative_magic_coeffs.items():
                assert amagic_coeff == 0 or (amagic_coeff == magic_coeffs.get(key, 0) and rc4.is_principal)
                by_schub_dict[key] = by_schub_dict.get(key, 0) + amagic_coeff * Sx(rc4.perm)

            # term by term only
            for key, term_coeff in term_coeffs.items():
                assert term_coeff >= 0
                by_term_dict[key] = by_term_dict.get(key, 0) + term_coeff * Sx(rc4.polyvalue(x))

            # BIJECTION

    for (perm1, perm2), elem in by_coefficient_dict.items():
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
        assert by_schub_dict.get((perm1, perm2), 0) == product
        assert by_term_dict.get((perm1, perm2), 0) == product
        print(f"Success {perm1.trimcode} {perm2.trimcode}")
    exit()


if __name__ == "__main__":
    main()
