"""``grothmult_q`` console script: products of quantum Grothendieck polynomials.

Quantum K-theory Pieri rule (Naito--Sagaki); ``--mult`` is not yet supported.

Example:
    grothmult_q 3 1 2 - 2 1 3
    grothmult_q 1 3 2 4 - 2 4 1 3 --parabolic 2 2
"""

import sys

from schubmult import Gx, Permutation, uncode
from schubmult.mult.groth_quantum import grothmult_q
from schubmult.mult.groth_quantum_double import apply_kato
from schubmult.symbolic import expand, sympify
from schubmult.utils.argparse import schub_argparse
from schubmult.utils.perm_utils import is_parabolic


def _parabolic_blocks(blocks):
    """``--parabolic`` block sizes -> (1-indexed generators of ``W_P``, ambient ``n``)."""
    parabolic_index = []
    start = 0
    for size in blocks:
        end = start + int(size)
        parabolic_index += list(range(start + 1, end))
        start = end
    return parabolic_index, start


def _check_parabolic_inputs(perms, parabolic_index, n):
    """Print a message and return ``False`` unless every input indexes a Schubert class of ``G/P``:
    a minimal coset representative of ``S_n / W_P`` (increasing on each block)."""
    for perm in perms:
        if perm.inv > 0 and (max(perm) > n or not is_parabolic(perm, parabolic_index)):
            print(f"{perm} is not a minimal coset representative for --parabolic blocks (must be increasing on each block of positions and lie in S_{n}).")
            return False
    return True


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
    parameter ``beta``). With ``--parabolic g1 g2 ...`` (block sizes of a parabolic subgroup
    ``S_{g1} x S_{g2} x ...``), the result is projected to ``QK(G/P)`` by the ``beta``-homogeneous
    form of Kato's ring homomorphism (`schubmult.mult.groth_quantum_double.apply_kato`); the input
    permutations must then be minimal coset representatives (increasing on each block), since only
    those index Schubert classes of ``G/P``. ``--mult`` is not yet supported and causes an early
    exit. Returns the raw result dict when the caller passes a ``None`` formatter (e.g. from
    tests); otherwise prints and returns ``None``.
    """
    if argv is None:
        argv = sys.argv
    try:
        args, formatter = schub_argparse(
            "grothmult_q",
            "Compute products of quantum Grothendieck polynomials",
            argv=argv[1:],
            quantum=True,
        )

        perms = args.perms

        if args.mult:
            print("--mult is not supported for quantum Grothendieck polynomials.")
            return 1

        ascode = args.ascode
        pr = args.pr
        beta = Gx._beta
        parabolic_index, n = _parabolic_blocks(args.parabolic)

        if ascode:
            perms = [uncode(perm) for perm in perms]
        else:
            perms = [Permutation(perm) for perm in perms]

        if parabolic_index and not _check_parabolic_inputs(perms, parabolic_index, n):
            return 1

        coeff_dict = {perms[0]: 1}
        for perm in perms[1:]:
            coeff_dict = grothmult_q(coeff_dict, perm, beta)

        if parabolic_index:
            coeff_dict = apply_kato(coeff_dict, parabolic_index, n=n)

        if pr or formatter is None:
            raw_result_dict = _display_full(coeff_dict, args, formatter)
        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
