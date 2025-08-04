from schubmult import *
from schubmult.rings import *
from schubmult.rings.free_algebra_basis import JBasis, NElementaryBasis
from sage.all import Compositions, Integer, QuasiSymmetricFunctions, ZZ

perms = Permutation.all_permutations(5)

for perm in perms:
    print(f"{perm.code}: {ASx(perm).kill_zero().change_basis(JBasis)}")