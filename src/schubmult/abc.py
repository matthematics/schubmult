from schubmult.rings.symmetric_polynomials.complete_sym import CompleteSym
from schubmult.rings.symmetric_polynomials.elem_sym import ElemSym
from schubmult.rings.variables import GeneratingSet

e = lambda *x: ElemSym(*x)
h = lambda *x: CompleteSym(*x)
x = GeneratingSet("x")
y = GeneratingSet("y")
z = GeneratingSet("z")
