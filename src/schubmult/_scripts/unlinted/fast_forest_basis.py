from schubmult import *
from schubmult.abc import *
from schubmult.rings.schubert.schubert_ring import SingleSchubertRing
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic import sympify

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    
    w0 = Permutation.w0(n)

    for perm in perms:
        pass