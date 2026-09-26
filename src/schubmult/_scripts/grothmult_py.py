"""``grothmult_py`` console script: products of (ordinary) Grothendieck polynomials.

Example:
    grothmult_py 3 1 2 - 2 1 3
    grothmult_py --code 2 0 - 1 0
"""

import sys

from schubmult import Gx, Permutation, uncode
from schubmult.abc import x
from schubmult.mult.groth import grothmult_py, mult_poly_groth
from schubmult.symbolic import sympify
from schubmult.utils.argparse import schub_argparse


def main(argv=None):
    """Entry point for the ``grothmult_py`` console script.

    Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
    their Grothendieck polynomials via `schubmult.mult.groth.grothmult_py`, and prints
    the resulting coefficient dictionary ``{Permutation: coefficient}`` (coefficients are
    Laurent polynomials in the K-theory parameter ``beta``). Returns the raw result dict
    when the caller passes a ``None`` formatter (e.g. from tests); otherwise prints and
    returns ``None``.
    """
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
        beta = Gx._beta

        if ascode:
            perms = [Permutation(uncode(perm)) for perm in perms]
        else:
            perms = [Permutation(perm) for perm in perms]
        perms.sort(reverse=True, key=lambda perm: sum((~perm).theta()) - perm.inv)

        coeff_dict = {Permutation([*perms[0]]): 1}

        for perm in perms[1:]:
            coeff_dict = grothmult_py(coeff_dict, Permutation(perm), beta)
        if mult:
            mul_exp = sympify(mulstring)
            coeff_dict = mult_poly_groth(coeff_dict, mul_exp, x, beta)

        if pr or formatter is None:
            for perm, val in coeff_dict.items():
                if val != 0:
                    raw_result_dict[perm] = val
                    if formatter:
                        print(f"{val}  {str(perm)}")

        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))