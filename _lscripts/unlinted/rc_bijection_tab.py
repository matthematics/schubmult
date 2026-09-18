from schubmult import *

if __name__ == "__main__":
    import sys
    
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    checkset = set()
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            tup = (rc,)
            old_tup = tup
            while not tup[0].is_full_grass:
                nepple = tup[0].squash_decomp()
                tup = (nepple[0].resize(len(nepple[0]) - 1), nepple[1], *tup[1:])
                if len(tup) == len(old_tup):
                    raise ValueError(f"Squash decomposition did not reduce size: {tup}")
                old_tup = tup
            if tup in checkset:
                raise ValueError(f"Duplicate bijection found: {tup}")
            print("sar")
            checkset.add(tup)
    print(f"Found {len(checkset)} unique bijections.")