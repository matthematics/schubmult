def crystal_isomorphic(c1, c2, cutoff=None):
    hw_1, _ = c1.to_highest_weight(length=cutoff)
    hw_2, _ = c2.to_highest_weight(length=cutoff)

    stack = [(c1, c2)]
    if cutoff is None:
        cutoff = c1.crystal_length()
    if hw_1.crystal_weight[:cutoff] != hw_2.crystal_weight[:cutoff]:
        return False
    while len(stack) > 0:
        c1_test, c2_test = stack.pop()
        for i in range(1, cutoff):
            c1_test0 = c1_test.lowering_operator(i)
            if c1_test0 is not None:
                c2_test0 = c2_test.lowering_operator(i)
                if c2_test0 is None:
                    return False
                stack.append((c1_test0, c2_test0))
    return True

if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.crystal_graph import CrystalGraph, CrystalGraphTensor
    from schubmult.schub_lib.schub_lib import elem_sym_perms
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms}

    hw_rcs = {}
    for perm in Permutation.all_permutations(n + 1):
        hw_rcs[perm] = set()
        for rc in RCGraph.all_rc_graphs(perm):
            hw_rcs[perm].add(rc.to_highest_weight()[0])
    
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n-1):
            for k in range(1, n):
                for p in range(1, k + 1):
                    monk_rc = next(iter(RCGraph.all_rc_graphs(Permutation([]).swap(k-1, k), n-1, weight=tuple([*([0] * (p - 1)),1,*([0] * (n - 1 - p))]))))
                    tensor = CrystalGraphTensor(rc, monk_rc)
                    good = False
                    up_perms = [pperm for pperm, L in elem_sym_perms(perm, 1, k) if L == 1]
                    results = set()
                    lv = [*rc.length_vector]
                    lv[p-1] += 1
                    for up_perm in up_perms:
                        for rc2 in RCGraph.all_rc_graphs(up_perm, length=n-1, weight=tuple(lv)):
                            if crystal_isomorphic(tensor, rc2, cutoff=k):
                                good = True
                                results.add(rc2)
                                break

                    if good:
                        print(f"Success {p=} {k=}")
                        pretty_print(rc)
                        print("Result:")
                        pretty_print(results)
                    else:
                        print(f"FAIL {p=} {k=}")
                        pretty_print(rc)
                    assert good
                    if len(results) != 1:
                        print("AMBIGUOUS")
                