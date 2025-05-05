from schubmult.rings.symmetric_polynomials.complete_sym import CompleteSym
from schubmult.rings.symmetric_polynomials.elem_sym import ElemSym
from schubmult.rings.variables import GeneratingSet


def e(*x):
    return ElemSym(*x)
def h(*x):
    return CompleteSym(*x)
x = GeneratingSet("x")
y = GeneratingSet("y")
z = GeneratingSet("z")
