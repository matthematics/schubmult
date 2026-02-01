from schubmult import *
from schubmult.symbolic import *
from schubmult.abc import *
from schubmult.utils.schub_lib import elem_sym_perms_op
def is_decomposable(w):
    for i in range(1, len(w) - 1):
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)), *list(range(i + 2, len(w))))
        if coset.inv == 0 and set(w_J.code[: i + 1]) != {0} and set(w_J.code[i + 2 :]) != {0}:
            return True
    return False

if __name__ == "__main__":
    import sys
    
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        # for i in perm.descents(zero_indexed=False):
        lower_perms = {}
        for i in range(2, len(perm.trimcode)):
            summ = S.Zero
            lower_perms = {}
            for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
                if len(rc[i-1]) != 0:
                    continue
                pullout = rc.pull_out_row(i)
                lower_perms[pullout.perm] = lower_perms.get(pullout.perm, set())
                lower_perms[pullout.perm].add(pullout)

        for perm2 in lower_perms:
            assert lower_perms[perm2] == RCGraph.all_rc_graphs(perm2, len(perm.trimcode) - 1), f"Error: missing RC graphs for permutation {perm2} from pull out of {perm} at row {i}."
