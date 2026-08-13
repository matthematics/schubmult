from schubmult import *
from schubmult.symbolic import expand

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    max_factors = {}
    maximizers = {}
    for k in range(1, n):
        max_factors[k] = {}
        maximizers[k] = {}
        for u in perms:
            
            last_count = 1
            factors = []
            for p in range(1, k + 1):
                product = DSx(u) * DSx(uncode([0]*(k - p) + [1]*p), "z")
                count = len({k for k,v in product.items() if expand(v) != 0})
                #assert count <= n, "" + str(u) + " " + str(k) + " " + str(count)
                factors = count / last_count
                if factors > max_factors[k].get(p, 0):
                    max_factors[k][p] = factors
                    maximizers[k][p] = u
                last_count = count

    for k in range(1, n):
        print(f"n={n} k={k} max factors: " + ", ".join(f"{p}: {max_factors[k][p]:.4f}" for p in range(1, k + 1)))
        print(f"n={n} k={k} maximizers: " + ", ".join(f"{p}: {maximizers[k][p]}" for p in range(1, k + 1)))