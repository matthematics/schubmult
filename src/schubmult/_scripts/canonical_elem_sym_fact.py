from schubmult import *

# coordinates RC

def elem_factor_from_rc(rc):
    from schubmult.utils.schub_lib import elem_sym_perms_op
    # find maximal downard path
    # target_perm = None
    descent = len(rc.perm.trimcode)
    if descent == 0:
        return {1: ()}
    target_perm = rc.perm
    while len(target_perm.trimcode) >= descent:
        results = elem_sym_perms_op(rc.perm, descent, descent)
        new_results = [permperm for permperm in results if len(permperm[0].trimcode) < descent]
        if len(new_results) == 0:
            descent += 1
            results = elem_sym_perms_op(rc.perm, descent, descent)
        target_perm = min(results, key=lambda p: (len(p[0].trimcode), p[0].inv), default=target_perm)[0]
    #except ValueError:
    #print(f"{rc.perm.trimcode}: {target_perm.trimcode}")
    if target_perm.inv == 0:
        return {descent: rc.resize(descent).length_vector}
    # awful slow way
    for old_rc in RCGraph.all_rc_graphs(target_perm, descent):
        weight_diff = tuple(rc.length_vector[i] - old_rc.length_vector[i] for i in range(descent))
        try:
            elem_sym_rc = next(iter(RCGraph.all_rc_graphs(uncode([0] * (descent - sum(weight_diff)) + [1] * sum(weight_diff)), descent, weight=weight_diff)))
        except StopIteration:
            continue
        try_rc = old_rc.squash_product(elem_sym_rc)
        if try_rc == rc.resize(len(try_rc)):
            return elem_factor_from_rc(old_rc) | {descent: weight_diff}
    #return {rc.perm: target_perm}
    raise ValueError(f"Could not find element factorization for RC graph {rc} with target permutation {target_perm}")

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            # dct = elem_factor_from_rc(rc)
            # new_rc = RCGraph([()])
            # assert sum([len(val) for key, val in dct.items()]) == rc.perm.inv
            # for elem_factor_index in range(1, n):
            #     new_rc = new_rc.resize(elem_factor_index)
            #     if elem_factor_index in dct:
            #         lst = dct[elem_factor_index]
            #         weight = [1 if i + 1 in lst else 0 for i in range(elem_factor_index)]
            #         elem_sym_rc = next(iter(RCGraph.all_rc_graphs(uncode([0] * (elem_factor_index - len(lst)) + [1] * len(lst)), elem_factor_index, weight=tuple(weight))))
            #         new_rc = new_rc.squash_product(elem_sym_rc)
            # if new_rc == rc:
            #     print("Yay")
            # else:
            #     print("Mismatch for permutation", perm)
            #     print("Original RC graph:", rc)
            #     print("Reconstructed RC graph:", new_rc)
            print(f"{repr(rc)}: {elem_factor_from_rc(rc)}")
