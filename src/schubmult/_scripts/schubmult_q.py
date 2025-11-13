import sys

from schubmult.symbolic import sstr, sympify
from schubmult.utils.argparse import schub_argparse

from schubmult import GeneratingSet, Permutation, apply_peterson_woodward, schubmult_q, schubmult_q_fast, uncode

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
        coeff_dict = apply_peterson_woodward(coeff_dict, parabolic_index)

    coeff_perms = list(coeff_dict.keys())
    coeff_perms.sort(key=lambda x: (x.inv, *x))

    for perm in coeff_perms:
        val = sympify(coeff_dict[perm]).expand()
        if val != 0:
            raw_result_dict[perm] = val
            if formatter:
                print(f"{sstr(perm)!s}  {formatter(val)}")
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
