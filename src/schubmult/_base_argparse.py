from argparse import ArgumentParser

def schub_argparse(quantum=False, double=False, yz=False):
    parser = ArgumentParser()

    parser.add_argument("perms", nargs="+", action="append")
    parser.add_argument("-", nargs="+", action="append", dest="perms")

    parser.add_argument(
        "-np", "--no-print", action="store_false", default=True, dest="pr"
    )

    parser.add_argument("--code", action="store_true", default=False, dest="ascode")

    parser.add_argument("--mult", nargs="+", required=False, default=None)

    parser.add_argument("--coprod", action="store_true", default=False)

    if yz:
        parser.add_argument(
            "--display-positive",
            action="store_true",
            default=False,
            dest="display_positive",
        )

        parser.add_argument(
            "--optimizer-message", action="store_true", default=False, dest="msg"
        )

        parser.add_argument(
            "--down",
            action="store_true",
            default=False,
        )

        parser.add_argument(
            "-nc", "--no-check", action="store_false", default=True, dest="check"
        )

        parser.add_argument(
            "--same-var", action="store_true", default=False, dest="same"
        )

        parser.add_argument("--expand", action="store_true", default=False, dest="expa")

        parser.add_argument("--norep", action="store_true", default=False)

    if quantum:
        parser.add_argument("--parabolic", nargs="+", required=False, default=[])

        parser.add_argument(
            "--nil-hecke", type=int, required=False, default=None, dest="nilhecke"
        )

        parser.add_argument(
            "--nil-hecke-apply",
            type=int,
            required=False,
            default=None,
            dest="nilhecke_apply",
        )

        parser.add_argument(
            "--basic-pieri", action="store_true", default=False, dest="slow"
        )

    args =  parser.parse_args()
    args.mulstring = ""

    if args.mult is not None:
        args.mulstring = " ".join(args.mult)
        args.mult = True

    for perm in args.perms:
        try:
            for i in range(len(perm)):
                perm[i] = int(perm[i])
        except Exception as e:
            print("Permutations must have integer values")
            raise e
    return args    