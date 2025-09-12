import math
import sys

from schubmult import ASx, uncode
from schubmult.rings import FA, WordBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, TensorModule
from schubmult.utils.perm_utils import artin_sequences


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
                    solution_module += coeff * coeff_double * coeff2 * TensorModule.ext_multiply(1 * rc, modmod)

    coprods = {}
    for (rc, (rc1, rc2)), coeff in solution_module.items():
        perm = rc.perm
        if len(rc1.perm) > n or len(rc2.perm) > n or len(perm) > n:
            continue
        # assert coeff >= 0 NOPE
        coprods[perm] = coprods.get(perm, 0) + coeff * (ASx @ ASx).ext_multiply(ASx(rc1.perm, len(rc1)), ASx(rc2.perm, len(rc2)))

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

    assert num_successes == math.factorial(n), f"Only {num_successes} successes, not the full group"

if __name__ == "__main__":
    main()
