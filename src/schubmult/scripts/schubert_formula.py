import sys

from schubmult import ASx, uncode
from schubmult.rings import FA, Sx, WordBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, TensorModule
from schubmult.utils.perm_utils import artin_sequences


def vector_sum(v1, v2):
    return tuple([a + b for a, b in zip(v1, v2)])


def main():
    n = int(sys.argv[1])

    unit_rc_module = RCGraphModule({RCGraph(): 1})

    # 100% positive!
    aseqs = artin_sequences(n - 1)
    perm_modules = {}
    ring = ASx @ ASx

    def asxt(rc):
        return (rc.perm, len(rc))

    unit_tensor_rc_module = TensorModule.ext_multiply(unit_rc_module, unit_rc_module)

    perm_modules3 = {}

    for aseq0 in aseqs:
        perm = uncode(aseq0)
        perm_words = ASx(perm, n - 1).change_basis(WordBasis)
        mod = ASx(perm, n - 1) * unit_rc_module
        for aseq1 in aseqs:
            posmod = FA(*aseq1) * unit_rc_module
            for rc3, coeff2 in mod.items():
                if len(rc3.perm) > n:
                    continue
                for rc4, coeff4 in posmod.items():
                    if rc3.perm != rc4.perm:
                        continue
                    for aseq in perm_words:
                        coeff = perm_words.get(aseq, 0)
                        if coeff == 0:
                            continue
                        mod2 = FA(*aseq).coproduct() * unit_tensor_rc_module
                        for (rc1, rc2), coeff1 in mod2.items():
                            if len(rc1.perm) > n or len(rc2.perm) > n:
                                continue
                            perm_modules[perm] = perm_modules.get(perm, 0) + coeff1 * ring((asxt(rc1), asxt(rc2)))
                            perm1, perm2 = rc1.perm, rc2.perm
                            perm_modules3[(perm1, perm2)] = perm_modules3.get((perm1, perm2), RCGraphModule()) + coeff * coeff1 * coeff2 * coeff4 * rc4

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
