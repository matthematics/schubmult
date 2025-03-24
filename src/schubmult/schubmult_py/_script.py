import sys

from symengine import sympify

from schubmult._base_argparse import schub_argparse
from schubmult.perm_lib import (
    inv,
    inverse,
    permtrim,
    theta,
    trimcode,
    uncode,
)
from schubmult.schubmult_py._funcs import mult_poly, schub_coprod, schubmult


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        args, formatter = schub_argparse(
            "schubmult_py", "Compute products of ordinary Schubert polynomials", argv=argv[1:],
        )

        mult = args.mult
        mulstring = args.mulstring

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
        coprod = args.coprod
        raw_result_dict = {}
        if coprod:
            if ascode:
                perms[0] = tuple(permtrim(uncode(perms[0])))
            pos = [*perms[1]]
            pos.sort()
            mperm = tuple(permtrim(perms[0]))

            coeff_dict = schub_coprod(mperm, pos)

            if pr or formatter is None:
                for firstperm, secondperm in coeff_dict:
                    val = coeff_dict[(firstperm, secondperm)]
                    if val != 0:
                        if ascode:
                            if formatter is None:
                                raw_result_dict[(tuple(trimcode(firstperm)), tuple(trimcode(secondperm)))] = val
                            else:
                                print(f"{val} {trimcode(firstperm)} {trimcode(secondperm)}")
                        elif formatter is None:
                            raw_result_dict[(firstperm, secondperm)] = val
                        else:
                            print(f"{val} {firstperm} {secondperm}")
        else:
            if ascode:
                for i in range(len(perms)):
                    perms[i] = tuple(permtrim(uncode(perms[i])))

            perms.sort(reverse=True, key=lambda x: sum(theta(inverse(x))) - inv(x))

            coeff_dict = {tuple(permtrim([*perms[0]])): 1}

            for perm in perms[1:]:
                coeff_dict = schubmult(coeff_dict, tuple(permtrim([*perm])))
            if mult:
                mul_exp = sympify(mulstring)
                coeff_dict = mult_poly(coeff_dict, mul_exp)

            if pr or formatter is None:
                for perm, val in coeff_dict.items():
                    if val != 0:
                        if ascode:
                            raw_result_dict[tuple(trimcode(perm))] = val
                            if formatter:
                                print(f"{val}  {trimcode(perm)}")
                        else:
                            raw_result_dict[tuple(perm)] = val
                            if formatter:
                                print(f"{val}  {perm}")
        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
