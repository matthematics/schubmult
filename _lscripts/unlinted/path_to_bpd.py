from schubmult import *
from schubmult.utils.schub_lib import elem_sym_perms_op
import numpy as np

if __name__ == "__main__":
    import sys
    perm = Permutation([int(a) for a in sys.argv[1:]])
    n = len(perm)
    path_dir = list(range(n-1, 0, -1))
    perms = Permutation.all_permutations(n)
    the_big_one = (perm) * Permutation.w0(n)
    paths = [[the_big_one]]
    for stinkbat in reversed(path_dir):
        new_paths = []
        for path in paths:
            last_perm = path[-1]
            for perm1, fatbat in elem_sym_perms_op(last_perm, stinkbat, stinkbat):
                if stinkbat == n - 1 and perm1.inv != 0:
                    continue
                new_paths.append(path + [perm1])
        paths = new_paths
    print(f"Permutation: {perm}, #paths: {len(paths)}")
    def pathweight(path):
        spit = list(reversed(path))
        weight = [0] * n
        for i in range(1, len(spit)):
            for j in range(n - i):
                if spit[i][j] == spit[i-1][j]:
                    weight[spit[i][j] - 1] += 1
        return tuple(weight)

    for path in paths:
        print(f"BANGFAT DUGONG {pathweight(path)}")
        spanky = np.array([[row[p] for p in range(n)] for row in reversed(path)]).transpose()
        print("\n".join(" ".join([str(a) for a in row]) for row in spanky))