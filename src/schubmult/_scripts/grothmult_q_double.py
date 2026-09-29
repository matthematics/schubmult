"""``grothmult_q_double`` console script: products of quantum double Grothendieck polynomials.

Conjectural quantum K-theory Pieri rule; ``--display-positive`` and
``--nil-hecke``/``--nil-hecke-apply`` are not yet supported.

Example:
    grothmult_q_double 3 1 2 - 2 1 3 --mixed-var
    grothmult_q_double 1 3 2 4 - 2 4 1 3 --parabolic 2 2
"""

import sys

from schubmult import GeneratingSet, Permutation, uncode
from schubmult._scripts.grothmult_q import _check_parabolic_inputs, _parabolic_blocks
from schubmult.mult.groth_double import normalize_coeff
from schubmult.mult.groth_quantum_double import apply_kato, grothmult_q_double
from schubmult.symbolic import S
from schubmult.utils.argparse import schub_argparse


def main(argv=None):
    """Entry point for the ``grothmult_q_double`` console script.

    Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
    their quantum double Grothendieck polynomials via
    `schubmult.mult.groth_quantum_double.grothmult_q_double`, and prints the resulting
    coefficient dictionary ``{Permutation: coefficient}`` in the ``y``/``z`` coefficient
    variables and quantum parameters ``q_i``. With ``--simplify`` coefficients are printed in
    cancelled form ``numer / prod (1 + beta*y_i)**e``. With ``--parabolic g1 g2 ...`` (block
    sizes of a parabolic subgroup) the result is projected to ``QK_T(G/P)`` by the
    ``beta``-homogeneous form of Kato's ring homomorphism
    (`schubmult.mult.groth_quantum_double.apply_kato`); inputs must be minimal coset
    representatives. ``--display-positive``, ``--nil-hecke``, and ``--mult`` are not yet supported
    and cause an early exit.
    Returns the raw result dict when the caller passes a ``None`` formatter (e.g. from
    tests); otherwise prints and returns ``None``.
    """
    if argv is None:
        argv = sys.argv

    try:
        var2 = GeneratingSet("y")
        var3 = GeneratingSet("z")
        sys.setrecursionlimit(1000000)

        args, formatter = schub_argparse(
            "grothmult_q_double",
            "Compute coefficients of products of quantum double Grothendieck polynomials in the same or different sets of coefficient variables (conjectural quantum K-Pieri rule)",
            argv=argv[1:],
            yz=True,
            quantum=True,
            coprod=False,
        )

        for flag, given in (
            ("--display-positive", args.display_positive),
            ("--nil-hecke", args.nilhecke is not None),
            ("--nil-hecke-apply", args.nilhecke_apply is not None),
            ("--mult", args.mult),
        ):
            if given:
                print(f"{flag} is not supported for grothmult_q_double yet.")
                return 1

        perms = args.perms
        ascode = args.ascode
        Permutation.print_as_code = ascode
        same = args.same
        pr = args.pr

        if same:
            var3 = var2

        if ascode:
            perms = [uncode(perm) for perm in perms]
        else:
            for i in range(len(perms)):
                if len(perms[i]) < 2 and (len(perms[i]) == 0 or perms[i][0] == 1):
                    perms[i] = Permutation([])
                perms[i] = Permutation(perms[i])

        parabolic_index, n = _parabolic_blocks(args.parabolic)
        if parabolic_index and not _check_parabolic_inputs(perms, parabolic_index, n):
            return 1

        coeff_dict = {perms[0]: 1}
        for perm in perms[1:]:
            coeff_dict = grothmult_q_double(coeff_dict, perm, var2, var3)

        if parabolic_index:
            coeff_dict = apply_kato(coeff_dict, parabolic_index, n=n)

        if args.simplify:
            coeff_dict = {perm: normalize_coeff(val, var2) for perm, val in coeff_dict.items()}

        coeff_perms = [perm for perm, val in coeff_dict.items() if val != S.Zero]
        coeff_perms.sort(key=lambda x: (-abs(perms[0].inv + perms[1].inv - x.inv), *x))
        width = max([len(str(perm)) for perm in coeff_perms]) if coeff_perms else 0

        raw_result_dict = {}
        for perm in coeff_perms:
            val = coeff_dict[perm]
            raw_result_dict[perm] = val
            if pr and formatter:
                print(f"{str(perm)!s:>{width}}  {formatter(val)}", flush=True)

        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
