from schubmult import *


if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    for perm in perms:
        for k in range(1, len(perm.trimcode)):
            for p in range(1, k + 1):
                for bpd in BPD.all_bpds(perm, len(perm.trimcode)):
                    rc = bpd.to_rc_graph()
                    elem_perm = ~uncode([0] * (k - p) + [p])
                    for gg in RCGraph.all_rc_graphs(elem_perm, len(perm.trimcode)):
                        da_bottom, top = rc.vertical_cut(k)
                        new_rc = RCGraph([*da_bottom.squash_product(gg.rowrange(0,k))])
                        da_bpd = BPD.from_rc_graph(new_rc).resize(len(perm))
                        grid = bpd.resize(len(perm))._grid.copy()
                        grid[:k, :da_bpd.cols] = da_bpd._grid[:k, :]
                        spom_bpd = BPD(grid)
                        spom_bpd.rebuild()
                        print(f"{perm} , k={k} , p={p} : {spom_bpd}")
                        assert spom_bpd.is_valid
                        