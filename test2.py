from schubmult import *
from schubmult.abc import *
coeff = (DSx([4,1,5,2,3],elem_sym=True)*DSx([5,1,4,2,3],"z"))[Permutation([6,2,5,1,3,4])]
coeff2 = coeff.expand(mul=True, func=False)
from sympy import sympify
f = coeff2.args[0].args[1]
print(f)
g = f.split_out_vars([y[5]], [z[2],z[3]])
print(f"{f=}")
print(f"{g=}")
import sympy
print(f"{sympy.expand(f-g, func=True)}")