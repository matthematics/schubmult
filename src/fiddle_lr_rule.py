import sys

from schubmult import ASx, Permutation
from schubmult.rings.rc_graph_module import RCGraph, try_lr_module

# THE ASX TABLEAU TELL YOU WHAT TO DO


def main():
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    rc_graphs = {perm: RCGraph.all_rc_graphs(perm, n - 1) for perm in perms}

    rc_graphs_by_weight = {}

    for perm, rcs in rc_graphs.items():
        for rc in rcs:
            dct = rc_graphs_by_weight.get(perm, {})
            st = dct.get(rc.length_vector(), set())
            st.add(rc)
            dct[rc.length_vector()] = st
            rc_graphs_by_weight[perm] = dct

    for perm in perms:
        try_mod = try_lr_module(perm)
        elem = try_mod.asdtype(ASx @ ASx)

        check = ASx(perm).coproduct()
        try:
            assert all(v == 0 for v in (elem - check).values())
        except AssertionError:
            print(f"Fail on {perm}")
            raise

        print(f"Success {perm.trimcode}")
        print(try_mod)


if __name__ == "__main__":
    main()
