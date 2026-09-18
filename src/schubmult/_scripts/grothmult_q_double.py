import sys

from schubmult import GeneratingSet, Permutation, uncode
from schubmult.mult.groth_quantum_double import grothmult_q_double
from schubmult.symbolic import S, sstr
from schubmult.utils.argparse import schub_argparse


def main(argv=None):
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
            ("--parabolic", args.parabolic),
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

        coeff_dict = {perms[0]: 1}
        for perm in perms[1:]:
            coeff_dict = grothmult_q_double(coeff_dict, perm, var2, var3)

        coeff_perms = [perm for perm, val in coeff_dict.items() if val != S.Zero]
        coeff_perms.sort(key=lambda x: (-abs(perms[0].inv + perms[1].inv - x.inv), *x))
        width = max([len(sstr(perm)) for perm in coeff_perms]) if coeff_perms else 0

        raw_result_dict = {}
        for perm in coeff_perms:
            val = coeff_dict[perm]
            raw_result_dict[perm] = val
            if pr and formatter:
                print(f"{sstr(perm)!s:>{width}}  {formatter(val)}", flush=True)

        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
