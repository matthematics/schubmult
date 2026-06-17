import sys

from symengine import Symbol

from schubmult import Permutation, uncode
from schubmult.abc import x
from schubmult.symbolic import sstr, sympify
from schubmult.symbolic.poly.schub_poly import groth_dict_to_poly, groth_mul_full, to_groth
from schubmult.symbolic.poly.variables import ZeroGeneratingSet
from schubmult.utils.argparse import schub_argparse


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        args, formatter = schub_argparse(
            "grothmult_py",
            "Compute products of Grothendieck polynomials",
            argv=argv[1:],
            coprod=False,
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
        raw_result_dict = {}
        zz = ZeroGeneratingSet()
        beta = Symbol("\u03B2")

        if ascode:
            perms = [Permutation(uncode(perm)) for perm in perms]
        else:
            perms = [Permutation(perm) for perm in perms]
        perms.sort(reverse=True, key=lambda perm: sum((~perm).theta()) - perm.inv)

        coeff_dict = {Permutation([*perms[0]]): 1}

        for perm in perms[1:]:
            coeff_dict = groth_mul_full(coeff_dict, Permutation(perm), x, zz, beta)
        if mult:
            mul_exp = sympify(mulstring)
            coeff_dict = to_groth(groth_dict_to_poly(coeff_dict, x, zz, beta) * mul_exp, x, zz, beta)

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