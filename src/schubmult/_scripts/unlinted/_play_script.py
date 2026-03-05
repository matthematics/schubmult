from schubmult import *

if __name__ == "__main__":
    import itertools
    import sys
    from schubmult.utils.perm_utils import weak_compositions
    from schubmult.rings.polynomial_algebra import *
    n = int(sys.argv[1])
    
    weaks = weak_compositions(n, n - 1)

    def num_descs(comp):
        return sum(1 for i in range(len(comp) - 1) if comp[i] > comp[i + 1])
    for comp in sorted(weaks):
        if num_descs(comp) > 1:
            continue
        pinto = Sx.from_expr(Forest(*comp).expand())
        if len(pinto) <= 1:
            print(comp)
            print(pinto)
            print(uncode(comp).is_vexillary)



    