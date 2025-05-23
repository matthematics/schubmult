from symengine import var

from schubmult.rings.variables import GeneratingSet
from schubmult.symmetric_polynomials import CompleteSym, E, E_q, ElemSym, FactorialCompleteSym, FactorialElemSym, H, QFactorialElemSym, e, h  # noqa: F401

var("x_(1:100)")
var("y_(1:100)")
var("z_(1:100)")
var("q_(1:100)")

x = GeneratingSet("x")
y = GeneratingSet("y")
z = GeneratingSet("z")
q = GeneratingSet("q")
