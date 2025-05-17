import sys  # noqa: F401
from argparse import SUPPRESS, ArgumentParser, RawDescriptionHelpFormatter

from schubmult.symbolic import init_printing, latex, pretty, sstr, sympify_sympy

# from sympy import init_printing
from schubmult.utils.logging import init_logging

# Indexed.sympify_sympystr = lambda x, p: f"{p.doprint(x.args[0])}_{x.args[1]}"


def schub_argparse(prog_name, description, argv, quantum=False, yz=False):
    parser = ArgumentParser(
        prog=prog_name,
        description=description,
        epilog=f"""Example:
        {prog_name} 5 1 7 3 2 6 4 - 2 1 6 3 5 4 {"" if not yz else "[ --display-positive]"}
            or equivalently
        {prog_name} --code 4 0 4 1 0 1 - 1 0 3 0 1 {"" if not yz else "[ --display-positive]"}
        {"    or alternatively" if not quantum else ""}
        {prog_name if not quantum else ""} {"--coprod --code 2 0 3 0 1 - 2 4 " if not quantum else ""} {"" if not yz or quantum else "[ --display-positive]"}
    """,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "perms",
        nargs="+",
        action="append",
        help="Space-delimited permutations separated by hyphens, e. g. 3 4 1 2 - 5 1 2 4 3",
        metavar="hyphen-separated list of perms",
    )
    parser.add_argument("-", nargs="+", action="append", dest="perms", help=SUPPRESS)

    parser.add_argument(
        "-np",
        "--no-print",
        action="store_false",
        default=True,
        dest="pr",
        help="Compute the result but do not print it",
    )

    parser.add_argument(
        "--code",
        action="store_true",
        default=False,
        dest="ascode",
        help="Permutations represented by the Lehmer code",
    )

    parser.add_argument(
        "--mult",
        nargs="+",
        required=False,
        default=None,
        help="Some additional terms in the ring to multiply by",
    )

    if not quantum:
        parser.add_argument(
            "--coprod",
            action="store_true",
            default=False,
            help="Compute the coproduct (different syntax: one permutation, then a hyphen, the variable index positions to split on)",
        )

    if yz:
        parser.add_argument(
            "--display-positive",
            action="store_true",
            default=False,
            dest="display_positive",
            help="Display the result in terms of the positive roots, or if mixed variable attempt to display the result as a positive algebraic combination of terms of the form y_i - z_j",
        )

        parser.add_argument(
            "--optimizer-message",
            action="store_true",
            default=False,
            dest="msg",
            help="Display debug output during integer optimization for --display-positive",
        )

        parser.add_argument(
            "--down",
            action="store_true",
            default=False,
            help="Reverse multiplication",
        )

        parser.add_argument(
            "-nc",
            "--no-check",
            action="store_false",
            default=True,
            dest="check",
            help="Do not check if positive result matches the original",
        )

        parser.add_argument(
            "--mixed-var",
            action="store_false",
            dest="same",
            help="Used mixed variables y and z",
        )

        parser.add_argument(
            "--expand",
            action="store_true",
            default=False,
            dest="expa",
            help="Expand the output rather than leaving it as originally computed (slow)",
        )

    if quantum:
        parser.add_argument(
            "--parabolic",
            nargs="+",
            required=False,
            default=[],
            help="Generators of the parabolic subgroup to compute quantum coeffs for",
        )

        parser.add_argument(
            "--basic-pieri",
            action="store_true",
            default=False,
            dest="slow",
            help="Do not apply conjectural computation optimization to quantum",
        )
        if yz:
            parser.add_argument(
                "--nil-hecke",
                type=int,
                required=False,
                default=None,
                dest="nilhecke",
                metavar="N",
                help="Substitute up to N of Fomin-Gelfand-Postnikov commuting difference operators",
            )

            parser.add_argument(
                "--nil-hecke-apply",
                type=int,
                required=False,
                default=None,
                dest="nilhecke_apply",
                metavar="N",
                help="Substitute commuting difference operators for perm1, then apply to Schub indexed by perm2",
            )

    disp_mode_list = ["basic", "pretty", "latex", "raw"]

    if not yz and not quantum:
        disp_mode_list = ["basic", "raw"]

    parser.add_argument(
        "--display-mode",
        type=str,
        required=False,
        choices=disp_mode_list,
        default="basic",
        dest="disp_mode",
        help="Method of displaying the output. Default basic",
    )

    parser.add_argument(
        "-g",
        action="store_true",
        dest="gen",
        help=SUPPRESS,
    )

    parser.add_argument(
        "-debug",
        action="store_true",
        dest="debug",
        help=SUPPRESS,
    )

    parser.add_argument(
        "--secret",
        action="store_true",
        dest="secret",
        help=SUPPRESS,
    )

    args = parser.parse_args(argv)
    args.mulstring = ""

    if args.mult is not None:
        args.mulstring = " ".join(args.mult)
        args.mult = True
    else:
        args.mult = False

    for perm in args.perms:
        try:
            for i in range(len(perm)):
                perm[i] = int(perm[i])
        except Exception as e:
            print("Permutations must have integer values")
            raise e

    if args.gen:
        import json

        argv.pop(argv.index("-g"))
        args.__dict__["cmd_line"] = [prog_name, *argv]
        del args.__dict__["gen"]
        cmd = " ".join(args.cmd_line)
        cmd = cmd.replace("--", "").replace(" - ", "T").replace(" ", "_")
        with open(f"{cmd}.json", "w") as js:
            json.dump(args.__dict__, js, ensure_ascii=False, indent=1)
        exit(0)

    init_printing()

    if args.disp_mode == "latex":
        formatter = (  # noqa: E731
            lambda bob, width=None: latex(sympify_sympy(bob), order="old").replace("\\left", "").replace("\\right", "")
        )
    elif args.disp_mode == "pretty":
        # pretty we need to keep centered
        formatter = (  # noqa: E731
            lambda bob, width=None: pretty(sympify_sympy(bob), order="old")
            if width is None
            else pretty(sympify_sympy(bob), order="rev-lex" if args.same else "none", use_unicode=False).replace("\n", "\n" + " ".join(["" for i in range(width)]))
        )
    elif args.disp_mode == "basic":
        formatter = lambda bob, width=None: sstr(sympify_sympy(bob), order="old")  # , order="rev-lex" if args.same else "none")  # noqa: E731
    elif args.disp_mode == "raw":
        formatter = None
    init_logging(debug=args.debug)
    return args, formatter
