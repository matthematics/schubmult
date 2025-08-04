from schubmult import *
from schubmult.rings import *
from schubmult.symbolic import *
Permutation.print_as_code=False
Spoing = PolynomialAlgebra(basis=SchubertPolyBasis(5))
phantasm = Spoing(uncode([0,2,0,2,3])).change_basis(SepDescPolyBasis(numvars=5,k=4))
print(f"{phantasm=}")
phantasm2 = Spoing(uncode([0,2,0,2,3])).change_basis(ElemSymPolyBasis(numvars=5))
print(f"{phantasm2=}")
Fong = FreeAlgebra(basis=SeparatedDescentsBasis(3))
print(f"{Fong.from_dict({(uncode([0, 1, 2, 1]),uncode([2]),5): S.One}).coproduct()=}")