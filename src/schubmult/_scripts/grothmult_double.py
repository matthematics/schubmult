import sys

from schubmult import GeneratingSet, Permutation, uncode
from schubmult.mult.groth_double import grothmult_double, mult_poly_groth_double
from schubmult.symbolic import sstr, sympify, expand
from schubmult.utils.argparse import schub_argparse
from schubmult.utils.logging import get_logger

logger = get_logger(__name__)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        var2 = GeneratingSet("y")
        var3 = GeneratingSet("z")
        sys.setrecursionlimit(1000000)

        args, formatter = schub_argparse(
            "grothmult_double",
            "Compute coefficients of products of double Grothendieck polynomials in the same or different sets of coefficient variables",
            argv=argv[1:],
            yz=True,
            coprod=False,
        )

        if args.display_positive:
            print("--display-positive is not supported for grothmult_double")
            return 1

        mult = args.mult
        mulstring = args.mulstring

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

        coeff_dict = {perms[0]: 1}

        for perm in perms[1:]:
            coeff_dict = grothmult_double(coeff_dict, perm, var2, var3)

        if mult:
            mul_exp = sympify(mulstring)
            coeff_dict = mult_poly_groth_double(coeff_dict, mul_exp, var2, var3)

        raw_result_dict = {}
        if pr or formatter is None:
            width = max([len(sstr(perm)) for perm in coeff_dict]) if coeff_dict else 0
            coeff_perms = list(coeff_dict.keys())
            coeff_perms.sort(key=lambda x: (-abs(perms[0].inv + perms[1].inv - x.inv), *x))

            for perm in coeff_perms:
                val = coeff_dict[perm]
                if expand(val) != 0:
                    raw_result_dict[perm] = val
                    if formatter:
                        print(f"{sstr(perm)!s:>{width}}  {formatter(val)}")

        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
