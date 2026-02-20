from schubmult import *
import numpy as np

def tableau_lift(nilp):
    grd = np.full((len(nilp.perm), len(nilp.perm)), dtype=int, fill_value=None)
    for col in len(nilp.shape[0]):
        for row in range(len(nilp.shape) - 1, -1, -1):
            if nilp[row, col] is not None:
                pass
            

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    for perm in Permutation.all_permutations(n):
        for rc in RCGraph.all_rc_graphs(perm):
            nilp = NilPlactic().ed_insert(*rc.perm_word)