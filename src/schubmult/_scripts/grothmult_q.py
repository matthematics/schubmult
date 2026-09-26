"""``grothmult_q`` console script: products of quantum Grothendieck polynomials.

Conjectural quantum K-theory Pieri rule; ``--parabolic`` and ``--mult`` are not yet
supported.

Example:
    grothmult_q 3 1 2 - 2 1 3
"""

import sys

from schubmult import Gx, Permutation, uncode
from schubmult.mult.groth_quantum import grothmult_q
from schubmult.symbolic import expand, sympify
from schubmult.utils.argparse import schub_argparse


def _display_full(coeff_dict, args, formatter):
    raw_result_dict = {}
    Permutation.print_as_code = args.ascode
    coeff_perms = sorted(coeff_dict.keys(), key=lambda x: (x.inv, *x))
    for perm in coeff_perms:
        val = expand(sympify(coeff_dict[perm]))
        if val != 0:
            raw_result_dict[perm] = val
            if formatter:
                print(f"{str(perm)!s}  {formatter(val)}")
    return raw_result_dict


def main(argv=None):
    """Entry point for the ``grothmult_q`` console script.

    Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
    their quantum Grothendieck polynomials via `schubmult.mult.groth_quantum.grothmult_q`,
    and prints the resulting coefficient dictionary ``{Permutation: coefficient}``
    (coefficients are polynomials in the quantum parameters ``q_i`` and the K-theory
    parameter ``beta``). ``--parabolic`` and ``--mult`` are not yet supported and cause
    an early exit. Returns the raw result dict when the caller passes a ``None``
    formatter (e.g. from tests); otherwise prints and returns ``None``.
    """
    if argv is None:
        argv = sys.argv
    try:
        args, formatter = schub_argparse(
            "grothmult_q",
            "Compute products of quantum Grothendieck polynomials (conjectural quantum K-Pieri rule)",
            argv=argv[1:],
            quantum=True,
        )

        perms = args.perms

        if args.parabolic:
            print("Parabolic quantum Grothendieck products are not supported yet.")
            return 1
        if args.mult:
            print("--mult is not supported for quantum Grothendieck polynomials.")
            return 1

        ascode = args.ascode
        pr = args.pr
        beta = Gx._beta

        if ascode:
            perms = [uncode(perm) for perm in perms]
        else:
            perms = [Permutation(perm) for perm in perms]

        coeff_dict = {perms[0]: 1}
        for perm in perms[1:]:
            coeff_dict = grothmult_q(coeff_dict, perm, beta)

        if pr or formatter is None:
            raw_result_dict = _display_full(coeff_dict, args, formatter)
        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
