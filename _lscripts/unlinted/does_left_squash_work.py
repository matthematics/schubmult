from schubmult import *
from schubmult.utils.schub_lib import all_grassmannian_rc_graphs

if __name__ == "__main__":
    import sys
    import itertools

    r = DualRCGraphRing()

    n = int(sys.argv[1])
    m = int(sys.argv[2])

    perms = Permutation.all_permutations(n)

    grasses = all_grassmannian_rc_graphs(n, m)

    #for perm1, perm2 in itertools.produc
    # t(perms, repeat=2):
    for grass_rc1, grass_rc2 in itertools.product(grasses, repeat=2):
        # prod = Sx(rc.perm) * Sx(grass_rc2.perm)
        # result = grass_rc1.left_squash(rc)
        # assert prod.get(result.perm, 0) != 0, f"Failure for {rc}, {grass_rc1}, {grass_rc2}, got {result.perm}, {prod=}, {result=}"
        result1 = grass_rc1.squash_product(grass_rc2)
        result2 = grass_rc1.left_squash(grass_rc2)
        assert result1 == result2, f"Failure for {grass_rc1}, {grass_rc2}, got {result1=}, {result2=}"
        print("potato bangfat")

