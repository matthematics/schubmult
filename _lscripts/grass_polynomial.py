#!/usr/bin/env python3
"""Command line front end for Grassmannian polynomial computations.

Usage::

    python grass_polynomial.py a1 a2 ... ak - b1 b2 ... bk

where both ``a1 <= a2 <= ... <= ak`` and ``b1 <= b2 <= ... <= bk`` are weakly
increasing lists of nonnegative integers of the same length.
"""

import sys
from argparse import SUPPRESS, ArgumentParser, RawDescriptionHelpFormatter


def _build_parser():
    parser = ArgumentParser(
        prog="grass_polynomial",
        description="Compute a Grassmannian polynomial from two weakly increasing lists of nonnegative integers.",
        epilog="""Example:
    grass_polynomial 0 2 3 - 1 1 4
""",
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "groups",
        nargs="+",
        action="append",
        help="Two hyphen-separated weakly increasing lists of nonnegative integers of equal length, e. g. 0 2 3 - 1 1 4",
        metavar="list1 - list2",
    )
    # Lets argparse split the command line on the bare hyphen.
    parser.add_argument("-", nargs="+", action="append", dest="groups", help=SUPPRESS)
    parser.add_argument(
        "-np",
        "--no-print",
        action="store_false",
        default=True,
        dest="pr",
        help="Compute the result but do not print it",
    )
    return parser


def _parse_list(parser, raw, label):
    try:
        vals = tuple(int(a) for a in raw)
    except ValueError:
        parser.error(f"{label} must consist of integers, got {' '.join(raw)}")

    if any(a < 0 for a in vals):
        parser.error(f"{label} must consist of nonnegative integers, got {' '.join(map(str, vals))}")

    if any(a > b for a, b in zip(vals, vals[1:])):
        parser.error(f"{label} must be weakly increasing, got {' '.join(map(str, vals))}")

    return vals


def parse_args(argv=None):
    """Parse ``argv`` and return ``(list1, list2, args)``.

    Both returned tuples are weakly increasing tuples of nonnegative ints of the
    same length.
    """
    if argv is None:
        argv = sys.argv[1:]
    parser = _build_parser()
    args = parser.parse_args(argv)

    groups = args.groups
    if len(groups) != 2:
        parser.error("expected exactly one hyphen separating the two lists")

    list1 = _parse_list(parser, groups[0], "the first list")
    list2 = _parse_list(parser, groups[1], "the second list")

    if len(list1) != len(list2):
        parser.error(f"the two lists must have the same length, got {len(list1)} and {len(list2)}")

    return list1, list2, args


def main(argv=None):
    list1, list2, args = parse_args(argv)

    result = None  # TODO: compute the Grassmannian polynomial from list1 and list2

    if args.pr:
        print(f"list1 = {list1}")
        print(f"list2 = {list2}")
        print(result)

    
    perm1 = lambda k: uncode([a * k for a in list1])
    perm2 = lambda k: uncode([a * k for a in list2])
    result = {perm1(1): Sx.one}
    block_size = perm2(1).mul_dominant()[0] - 1
    #return result
    genset = [0, Sx(uncode([1]))] + [Sx(uncode([0]* (i - 1) + [1])) - Sx(uncode([0] * i + [1])) for i in range(30)]
    zagong = [0 for _ in range(30)]
    MAX_T = 3
    t = Symbol("t")

    for k in range(1, MAX_T + 1):
        new_result = {uncode([a//k for a in key.trimcode]): (t**k) * v for key, v in schubmult_double({perm1(k): Sx.one}, perm2(k), zagong, genset).items() if all(a % k == 0 for a in key)}

if __name__ == "__main__":
    main()
