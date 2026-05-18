import argparse

from schubmult import Permutation, RCGraph, uncode


def check_disjoint_vs_squash_omega(n: int) -> tuple[bool, int, int]:
    """Compare omega_invariant[1] for disjoint_union vs squash_product.

    Loop over minimal-length RC graphs rc. For each rc, loop over elementary
    symmetric RC graphs of matching length and last descent, and verify:

        rc.disjoint_union(elem_sym_rc).omega_invariant[1]
        ==
        rc.squash_product(elem_sym_rc).omega_invariant[1]
    """
    checked = 0
    mismatches = 0

    perms = Permutation.all_permutations(n)
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm):
            length = len(rc)
            last_desc = len(rc.perm.trimcode)

            if last_desc <= 0:
                continue

            for p in range(1, last_desc + 1):
                elem_perm = uncode([0] * (last_desc - p) + [1] * p)
                for elem_sym_rc in RCGraph.all_rc_graphs(elem_perm, length):
                    checked += 1
                    disjoint = rc.disjoint_union(elem_sym_rc)
                    squashed = rc.squash_product(elem_sym_rc)
                    if disjoint.omega_invariant[1] != squashed.omega_invariant[1]:
                        mismatches += 1
                        print(
                            "Mismatch:",
                            f"perm={perm}",
                            f"rc={rc}",
                            f"p={p}",
                            f"elem_perm={elem_perm}",
                            f"elem_sym_rc={elem_sym_rc}",
                            f"disjoint_union={disjoint}",
                            f"squash_product={squashed}",
                            f"disjoint omega_invariant[1]={disjoint.omega_invariant[1]}",
                            f"squash omega_invariant[1]={squashed.omega_invariant[1]}",
                            sep="\n  ",
                            flush=True,
                        )

    return mismatches == 0, checked, mismatches


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check omega_invariant[1] for rc.disjoint_union(elem_sym_rc) vs rc.squash_product(elem_sym_rc)"
    )
    parser.add_argument("n", type=int, help="Permutation size")
    args = parser.parse_args()

    ok, checked, mismatches = check_disjoint_vs_squash_omega(args.n)
    if ok:
        print(f"All checks passed for n={args.n}. Cases checked: {checked}")
        raise SystemExit(0)

    print(f"Found {mismatches} mismatches for n={args.n} across {checked} checked cases")
    raise SystemExit(1)


if __name__ == "__main__":
    main()
