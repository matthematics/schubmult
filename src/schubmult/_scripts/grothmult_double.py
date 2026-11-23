"""``grothmult_double`` console script: products of double Grothendieck polynomials.

Example:
    grothmult_double 3 1 2 - 2 1 3 --display-positive
    grothmult_double --code 2 0 - 1 0 --mixed-var
"""

import sys

from schubmult import Gx, Permutation, uncode
from schubmult.abc import beta
from schubmult.mult.groth_double import grothmult_double, mult_poly_groth_double, normalize_coeff
from schubmult.mult.positivity import groth_posify_many
from schubmult.symbolic import S, expand, sympify, sympify_sympy
from schubmult.utils.argparse import grothendieck_cli_variables, schub_argparse
from schubmult.utils.logging import get_logger

logger = get_logger(__name__)


def _limit_memory():
    """Cap this process's address space so runaway cases die with MemoryError
    instead of taking down the host (WSL is especially fragile under OOM).
    No-op where ``resource`` is unavailable (Windows)."""
    try:
        import resource

        import psutil

        cap = max(int(psutil.virtual_memory().available * 0.7), 1024**3)
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        if soft == resource.RLIM_INFINITY or soft > cap:
            resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    except Exception:  # noqa: BLE001
        logger.debug("could not set memory limit", exc_info=True)


def _denominator_budget(den, var2):
    """Atom multiplicities ``{y_symbol: exp}`` of a denominator ``prod (1 + beta*y_i)**e``.

    Keys are symengine symbols (consistent with GeneratingSet identity).
    """
    from sympy import factor_list

    bs = sympify_sympy(beta)
    budget = {}
    _, factors = factor_list(sympify_sympy(den))
    for fac, mult in factors:
        syms = [s for s in fac.free_symbols if s != bs]
        if len(syms) != 1 or var2.index(syms[0]) == -1:
            raise ValueError(f"unexpected denominator factor {fac}")
        key = sympify(syms[0])
        budget[key] = budget.get(key, 0) + int(mult)
    return budget


def main(argv=None):
    """Entry point for the ``grothmult_double`` console script.

    Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
    their double Grothendieck polynomials via `schubmult.mult.groth_double.grothmult_double`,
    and prints the resulting coefficient dictionary ``{Permutation: coefficient}`` in the
    ``y``/``z`` coefficient variables and the K-theory parameter ``beta``. With
    ``--multiplicative-variables``, use ``beta=-1`` and ``1-1/y_i`` coefficient
    coordinates (also ``1-1/z_i`` for mixed variables). With
    ``--display-positive`` (only supported together with ``--mixed-var``), coefficients
    are rewritten in the positive FGL basis of ``z_i - y_j``, ``(1 + y_j)^{-1}``, and
    ``(1 + z_i)^{-1}`` via an integer program (see `schubmult.mult.positivity.groth_posify`);
    the integer programs of different coefficients are solved concurrently, up to ``--jobs`` at
    a time (default 1; see `schubmult.mult.positivity.groth_posify_many`). Returns the raw
    result dict when the caller passes a ``None`` formatter (e.g. from tests); otherwise
    prints and returns ``None``.
    """
    if argv is None:
        argv = sys.argv

    try:
        sys.setrecursionlimit(1000000)

        args, formatter = schub_argparse(
            "grothmult_double",
            "Compute coefficients of products of double Grothendieck polynomials in the same or different sets of coefficient variables",
            argv=argv[1:],
            yz=True,
            coprod=False,
        )

        if args.display_positive and args.same:
            print("--display-positive is only supported with --mixed-var for grothmult_double")
            return 1

        if args.display_positive:
            _limit_memory()

        mult = args.mult
        mulstring = args.mulstring

        perms = args.perms

        ascode = args.ascode
        Permutation.print_as_code = ascode
        pr = args.pr
        beta_value = sympify(args.beta) if args.beta is not None else Gx._beta

        var2, var3 = grothendieck_cli_variables(args)

        if ascode:
            perms = [uncode(perm) for perm in perms]
        else:
            for i in range(len(perms)):
                if len(perms[i]) < 2 and (len(perms[i]) == 0 or perms[i][0] == 1):
                    perms[i] = Permutation([])
                perms[i] = Permutation(perms[i])

        coeff_dict = {perms[0]: 1}

        for perm in perms[1:]:
            coeff_dict = grothmult_double(coeff_dict, perm, var2, var3, beta=beta_value)

        if mult:
            mul_exp = sympify(mulstring)
            coeff_dict = mult_poly_groth_double(coeff_dict, mul_exp, var2, var3, beta=beta_value)

        if args.simplify and not args.display_positive:
            if args.multiplicative_variables:
                coeff_dict = {perm: sympify(sympify_sympy(val).cancel()) for perm, val in coeff_dict.items()}
            else:
                coeff_dict = {perm: normalize_coeff(val, var2, beta=beta_value) if beta_value != 0 else expand(val) for perm, val in coeff_dict.items()}

        positive_coeff_dict = coeff_dict
        if args.display_positive and args.beta is not None:
            # The positivity solver needs the symbolic beta grading before specialization.
            positive_coeff_dict = {perms[0]: 1}
            for perm in perms[1:]:
                positive_coeff_dict = grothmult_double(positive_coeff_dict, perm, var2, var3, beta=Gx._beta)
            if mult:
                positive_coeff_dict = mult_poly_groth_double(positive_coeff_dict, mul_exp, var2, var3, beta=Gx._beta)

        # sort/filter up front so posified coefficients can stream out as each
        # one finishes (the LPs can take a long time); structural zero test only
        coeff_perms = [perm for perm, val in coeff_dict.items() if val != S.Zero]
        coeff_perms.sort(key=lambda x: (-abs(perms[0].inv + perms[1].inv - x.inv), *x))
        width = max([len(str(perm)) for perm in coeff_perms]) if coeff_perms else 0

        raw_result_dict = {}
        posified = None
        if args.display_positive:
            posified = groth_posify_many(
                [positive_coeff_dict[perm] for perm in coeff_perms],
                var2,
                var3,
                args.msg,
                max_workers=args.jobs,
                progress=lambda i: print(f"posify {i + 1}/{len(coeff_perms)}: {coeff_perms[i]!s}", file=sys.stderr),
            )
        try:
            for perm in coeff_perms:
                val = coeff_dict[perm]
                if posified is not None:
                    try:
                        # groth_posify verifies its own output exactly
                        val = next(posified)
                        if args.beta is not None:
                            val = sympify(val).subs(Gx._beta, beta_value)
                    except Exception:
                        import traceback

                        traceback.print_exc()
                        print(f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=}")
                        return 1
                raw_result_dict[perm] = val
                if pr and formatter:
                    print(f"{str(perm)!s:>{width}}  {formatter(val)}", flush=True)
        finally:
            # stops the solver processes of coefficients that are no longer needed
            if posified is not None:
                posified.close()

        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
