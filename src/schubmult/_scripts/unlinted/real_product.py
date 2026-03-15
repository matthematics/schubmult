from schubmult import *

def full_factorization(rc):
    if len(rc) == 1:
        return [rc]
    base, grass = rc.squash_decomp()
    return full_factorization(base.resize(len(base) - 1)) + [grass]

def is_full_grass(rc):
    if rc.perm.inv == 0 or rc.perm.descents() == {len(rc) - 1}:
        return True
    return False

def grass_left_multiply(grass_rc, rc2):
    # full
    assert is_full_grass(grass_rc)
    if len(rc2.perm) < len(grass_rc):
        return grass_rc.left_squash(rc2.resize(len(grass_rc)))
    if len(grass_rc) == len(rc2):
        return grass_rc.left_squash(rc2)
    if len(grass_rc) < len(rc2):
        base, grass = rc2.squash_decomp()
        return grass_left_multiply(grass_rc, base.resize(len(base.perm) - 1)).resize(len(rc2)).squash_product(grass)
    return grass_rc.left_squash(rc2.resize(len(grass_rc)))

def recursive_multiply(rc1, rc2):
    assert len(rc1) == len(rc2)
    
    factorization = full_factorization(rc1)
    result = rc2
    for factor in reversed(factorization):
        result = grass_left_multiply(factor, result)
    return result

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    for perm1, perm2 in itertools.product(perms, repeat=2):
        prodo = Sx(perm1) * Sx(perm2)
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n), RCGraph.all_rc_graphs(perm2, n)):
            if is_full_grass(rc1) and is_full_grass(rc2):
                faffle = rc1.squash_product(rc2)
            elif is_full_grass(rc1):
                faffle = rc1.left_squash(rc2)
            elif is_full_grass(rc2):
                faffle = rc1.squash_product(rc2)
            else:
                faffle = recursive_multiply(rc1, rc2)
            assert prodo.get(faffle.perm, 0) != 0, f"Failure for {rc1}, {rc2}, got {faffle.perm}, {prodo=}, {faffle=}"