import sys

from schubmult.symbolic import sstr, sympify
from schubmult.utils.argparse import schub_argparse

from schubmult import Permutation, mult_poly_py, schub_coprod_py, schubmult_py, theta, uncode


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        args, formatter = schub_argparse(
            "schubmult_py",
            "Compute products of ordinary Schubert polynomials",
            argv=argv[1:],
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
        Permutation.print_as_code = ascode
        pr = args.pr
        coprod = args.coprod
        raw_result_dict = {}
        if coprod:
            if ascode:
                perms[0] = uncode(perms[0])
            pos = [*perms[1]]
            pos.sort()
            mperm = Permutation(perms[0])

            coeff_dict = schub_coprod_py(mperm, pos)

            if pr or formatter is None:
                for firstperm, secondperm in coeff_dict:
                    val = coeff_dict[(firstperm, secondperm)]
                    if val != 0:
                        if formatter is None:
                            raw_result_dict[(firstperm, secondperm)] = val
                        else:
                            print(f"{val} {sstr(firstperm)} {sstr(secondperm)}")
        else:
            if ascode:
                for i in range(len(perms)):
                    perms[i] = Permutation(uncode(perms[i]))
            else:
                perms = [Permutation(perm) for perm in perms]
            perms.sort(reverse=True, key=lambda x: sum(theta(~x)) - x.inv)

            coeff_dict = {Permutation([*perms[0]]): 1}

            for perm in perms[1:]:
                coeff_dict = schubmult_py(coeff_dict, Permutation(perm))
            if mult:
                mul_exp = sympify(mulstring)
                coeff_dict = mult_poly_py(coeff_dict, mul_exp)

            if pr or formatter is None:
                for perm, val in coeff_dict.items():
                    if val != 0:
                        raw_result_dict[perm] = val
                        if formatter:
                            print(f"{val}  {sstr(perm)}")
        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
