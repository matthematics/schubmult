from schubmult import *
from schubmult.rings import *
from schubmult.rings.free_algebra_basis import JBasis, NElementaryBasis
from schubmult.symbolic import *
from schubmult.abc import *
from sage.all import Compositions, Integer, QuasiSymmetricFunctions, ZZ

J = FreeAlgebra(JBasis)

Elputz = FreeAlgebra(NElementaryBasis)

monom_to_j_dict = {}

for i in range(1, 7):
    for c in Compositions(Integer(i)):
        cc = tuple([int(a) for a in c])
        c = cc
        elem = FA(*cc).change_basis(JBasis)
        # add_elem = FA.zero
        # for (perm, n), coeff in elem.items():
        #     cd = perm.code
        #     if len(cd) < n or 0 in cd[:n]:
        #         continue
        #     add_elem[tuple(cd[:n])] = coeff
        # monom_to_j_dict[c] = add_elem
        monom_to_j_dict[c] = elem
        

# transpose

basis_dict = {}


QSym = QuasiSymmetricFunctions(ZZ)



M = QSym.M()
F = QSym.F()
E = QSym.E()
QS = QSym.QS()
YQS = QSym.YQS()
DOING = QS.dual()

for k, v in monom_to_j_dict.items():    
    for a, b in v.items():
        basis_dict[a] = basis_dict.get(a, QSym.zero()) + b*M[*k]


for comp in basis_dict:
    print(f"S({comp}) = {basis_dict[*comp]}")
    print(f"QS({comp}) = {M(QS[*comp])}")
    print()
    print(f"S({comp}) = {QS(basis_dict[*comp])}")
    print(f"S({comp}) = {YQS(basis_dict[*comp])}")
    print(f"S({comp}) = {F(basis_dict[*comp])}")
    print(f"S({comp}) = {E(basis_dict[*comp])}")
    print()
    #fingus = J(uncode([int(a) for a in comp]), len(comp)).change_basis(WordBasis)
    fingus = J(*[int(a) for a in comp]).change_basis(WordBasis)
    fingus.kill_zero()
    finker = fingus.bcoproduct()
    #finker2 = (FA@FA).zero
    # for (tup1, tup2), v in finker.items():
    #     finker2[(tuple([a - 1 for a in tup1]), tuple([a - 1 for a in tup2]))] = v
    #finker2 = FreeAlgebraBasis.change_tensor_basis(finker, SchubertBasis, SchubertBasis)
    finker2 = FreeAlgebraBasis.change_tensor_basis(finker, JBasis, JBasis)
    stink_elem = finker2
    # for ((perm1, n1), (perm2, n2)), coeff in finker2.items():
    #     cd1 = perm1.code
    #     if len(cd1) < n1 or 0 in cd1[:n1]:
    #         continue
    #     cd2 = perm2.code
    #     if len(cd2) < n2 or 0 in cd2[:n2]:
    #         continue
    #     stink_elem[(tuple(cd1[:n1]),tuple(cd2[:n2]))] = coeff
    print(f"{stink_elem}")
    print(f"{DOING[*comp].coproduct()}")
    biffdank = sympify(str(basis_dict[*comp].expand(10)))
    subs_dict = {Symbol(f"x{i}"): x[i + 1] for i in range(10)}
    biffdank = efficient_subs(biffdank, subs_dict)
    #print(f"{biffdank=}")
    Permutation.print_as_code = True
    print(Sx(biffdank))

