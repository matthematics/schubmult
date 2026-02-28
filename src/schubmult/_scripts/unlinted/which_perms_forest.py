from schubmult import *
from schubmult.rings.free_algebra import *
from schubmult.utils.perm_utils import weak_compositions

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    Forest = FreeAlgebra(ForestBasis)
    comps = weak_compositions(n-1, n - 1)
    r = RCGraphRing()
    for comp in comps:
        forr = r.from_free_algebra_element(Forest(*comp))
        new_bucket = r.from_dict({**forr})
        for rc, coeff in forr.items():
            if coeff < 0:
                try:
                    rc2 = next(iter([rc2 for rc2, coeff2 in new_bucket.items() if rc2.perm == rc.perm and coeff2 > 0]))
                except StopIteration:
                    continue
                del new_bucket[rc]
                new_bucket += coeff * r(rc2)
                # if new_bucket[rc2] == 0:
                #     del new_bucket[rc2]
        print(f"Newbubck bucket for {comp}: {new_bucket}") 