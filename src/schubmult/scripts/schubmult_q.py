import sys

import numpy as np
import sympy
from symengine import sympify

from schubmult import GeneratingSet, Permutation, check_blocks, factor_out_q_keep_factored, permtrim, q_vector, schubmult_q, schubmult_q_fast, uncode
from schubmult.perm_lib import (
    longest_element,
)
from schubmult.utils.argparse import schub_argparse
from schubmult.utils.perm_utils import (
    count_less_than,
    is_parabolic,
    omega,
)

q_var = GeneratingSet("q")


def _display_full(coeff_dict, args, formatter):
    raw_result_dict = {}
    ascode = args.ascode
    Permutation.print_as_code = ascode
    parabolic_index = []
    start = 0
    for i in range(len(args.parabolic)):
        end = start + int(args.parabolic[i])
        parabolic_index += list(range(start + 1, end))
        start = end
    parabolic = len(parabolic_index) != 0

    if parabolic:
            max_len = parabolic_index[-1] + 1
            # parabolic_index += list(range(parabolic_index[-1] + 2, max_len))
            w_P = longest_element(parabolic_index)
            # max_len = len(w_P)
            w_P_prime = Permutation([1, 2])
            coeff_dict_update = {}
            for w_1 in coeff_dict.keys():
                val = coeff_dict[w_1]
                q_dict = factor_out_q_keep_factored(val)
                for q_part in q_dict:
                    qv = q_vector(q_part)
                    w = w_1
                    good = True
                    parabolic_index2 = []
                    for i in range(len(parabolic_index)):
                        if omega(parabolic_index[i], qv) == 0:
                            parabolic_index2 += [parabolic_index[i]]
                        elif omega(parabolic_index[i], qv) != -1:
                            good = False
                            break
                    if not good:
                        continue
                    w_P_prime = longest_element(parabolic_index2)
                    if not check_blocks(qv, parabolic_index):
                        continue
                    w = (w * w_P_prime) * w_P
                    if not is_parabolic(w, parabolic_index):
                        continue

                    w = permtrim(w)
                    if len(w) > max_len:
                        continue
                    new_q_part = np.prod(
                        [q_var[index + 1 - count_less_than(parabolic_index, index + 1)] ** qv[index] for index in range(len(qv)) if index + 1 not in parabolic_index],
                    )
                    try:
                        new_q_part = int(new_q_part)
                    except Exception:
                        pass
                    q_val_part = q_dict[q_part]
                    coeff_dict_update[w] = coeff_dict_update.get(w, 0) + new_q_part * q_val_part
            coeff_dict = coeff_dict_update

    coeff_perms = list(coeff_dict.keys())
    coeff_perms.sort(key=lambda x: (x.inv, *x))

    for perm in coeff_perms:
        val = sympify(coeff_dict[perm]).expand()
        if val != 0:
            raw_result_dict[perm] = val
            if formatter:
                print(f"{sympy.sstr(perm)!s}  {formatter(val)}")
    return raw_result_dict


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        args, formatter = schub_argparse(
            "schubmult_q",
            "Compute products of quantum Schubert polynomials",
            argv=argv[1:],
            quantum=True,
        )

        perms = args.perms

        for perm in perms:
            try:
                for i in range(len(perm)):
                    perm[i] = int(perm[i])
            except Exception as e:
                print("Permutations must have integer values")
                raise e

        ascode = args.ascode
        pr = args.pr
        parabolic_index = [int(s) for s in args.parabolic]
        parabolic = len(parabolic_index) != 0
        slow = args.slow

        if parabolic and len(perms) != 2:
            print("Only two permutations supported for parabolic.")
            exit(1)

        if ascode:
            for i in range(len(perms)):
                perms[i] = uncode(perms[i])
        else:
            perms = [Permutation(perm) for perm in perms]

        # if parabolic:
        #     for i in range(len(parabolic_index)):
        #         index = parabolic_index[i] - 1
        #         if sg(index, perms[0]) == 1 or sg(index, perms[1]) == 1:
        #             print(
        #                 "Parabolic given but elements are not minimal length coset representatives.",
        #             )
        #             exit(1)

        coeff_dict = {perms[0]: 1}

        if not slow:
            for perm in perms[1:]:
                coeff_dict = schubmult_q_fast(coeff_dict, perm)
        else:
            for perm in perms[1:]:
                coeff_dict = schubmult_q(coeff_dict, perm)

        # if mult:
        #     mul_exp = sympify(mulstring)
        #     coeff_dict = mult_poly(coeff_dict, mul_exp)

        if pr or formatter is None:
            raw_result_dict = _display_full(coeff_dict, args, formatter)
        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv))
