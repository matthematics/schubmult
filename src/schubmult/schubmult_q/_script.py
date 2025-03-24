import sys

import numpy as np
from symengine import sympify

from schubmult._base_argparse import schub_argparse
from schubmult.perm_lib import (
    check_blocks,
    count_less_than,
    inv,
    is_parabolic,
    longest_element,
    mulperm,
    omega,
    permtrim,
    q_var,
    q_vector,
    sg,
    trimcode,
    uncode,
)
from schubmult.schubmult_q._funcs import (
    schubmult,
    schubmult_db,
)
from schubmult.schubmult_q_double import factor_out_q_keep_factored


def _display_full(coeff_dict, args, formatter):
    raw_result_dict = {}
    ascode = args.ascode
    parabolic_index = [int(s) for s in args.parabolic]
    parabolic = len(parabolic_index) != 0

    if parabolic:
        w_P = longest_element(parabolic_index)
        w_P_prime = [1, 2]
        coeff_dict_update = {}
        for w_1 in coeff_dict:
            val = coeff_dict[w_1]
            q_dict = factor_out_q_keep_factored(val)
            for q_part in q_dict:
                qv = q_vector(q_part)
                w = [*w_1]
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
                w = permtrim(mulperm(mulperm(w, w_P_prime), w_P))
                if not is_parabolic(w, parabolic_index):
                    continue

                w = tuple(permtrim(w))

                new_q_part = np.prod(
                    [
                        q_var[index + 1 - count_less_than(parabolic_index, index + 1)] ** qv[index]
                        for index in range(len(qv))
                        if index + 1 not in parabolic_index
                    ],
                )

                try:
                    new_q_part = int(new_q_part)
                except Exception:
                    pass
                q_val_part = q_dict[q_part]
                coeff_dict_update[w] = coeff_dict_update.get(w, 0) + new_q_part * q_val_part
        coeff_dict = coeff_dict_update

    coeff_perms = list(coeff_dict.keys())
    coeff_perms.sort(key=lambda x: (inv(x), *x))

    for perm in coeff_perms:
        val = sympify(coeff_dict[perm]).expand()
        if val != 0:
            if ascode:
                raw_result_dict[tuple(trimcode(perm))] = val
                if formatter:
                    print(f"{trimcode(perm)!s}  {formatter(val)}")
            else:
                raw_result_dict[perm] = val
                if formatter:
                    print(f"{perm!s}  {formatter(val)}")
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

        if parabolic:
            for i in range(len(parabolic_index)):
                index = parabolic_index[i] - 1
                if sg(index, perms[0]) == 1 or sg(index, perms[1]) == 1:
                    print(
                        "Parabolic given but elements are not minimal length coset representatives.",
                    )
                    exit(1)

        coeff_dict = {tuple(permtrim([*perms[0]])): 1}

        if not slow:
            for perm in perms[1:]:
                coeff_dict = schubmult_db(coeff_dict, tuple(permtrim([*perm])))
        else:
            for perm in perms[1:]:
                coeff_dict = schubmult(coeff_dict, tuple(permtrim([*perm])))

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
