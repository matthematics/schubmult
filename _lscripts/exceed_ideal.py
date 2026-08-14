from schubmult import *

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    bw = BoundedWCFactorAlgebra()
    seen = set()
    exideal = set()
    base_keys = set()
    for perm in perms:
        ge = bw.full_groth_elem(perm, n)
        for key in ge.keys():
            if key in seen:
                continue
            else:
                seen.add(key)
            wc = bw.key_to_wc_graph(key).normalize()
            if wc.excess > 0:
                exideal.add(key)
            else:
                base_keys.add(key)


    for key1 in exideal:
        wc1 = bw.key_to_wc_graph(key).normalize()
        for key2 in base_keys:
            pain = bw(key) * bw(key2)
            spain = next(iter(pain.keys()))
            wc2 = bw.key_to_wc_graph(key2).normalize()
            wcc = bw.key_to_wc_graph(spain).normalize()
            if wcc.excess > wc1.excess + wc2.excess:
                #print(f"Excess not additive for {key} * {key2}: {wc1.excess=} + {wc2.excess=} != {wcc.excess=}")
                raise ValueError(f"Excess spanish for {key} * {key2}:\n{wc1=} excess {wc1.excess}\n{wc2=} excess {wc2.excess}\n{wcc=} excess {wcc.excess}")
                # if wcc.excess < wc1.excess + wc2.excess:
                #     raise ValueError(f"Excess not superadditive for {key} * {key2}: {wc1.excess=} + {wc2.excess=} != {wcc.excess=}")