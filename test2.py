from schubmult import *
from schubmult.rings.elem_sym import ElemSym
from schubmult.schub_lib.double import schubmult_double_from_elems

def elem_sym(p, k, varl1, varl2):
    return ElemSym(p, k, varl1, varl2)

print(schubmult_double_from_elems({Permutation([4,3,1,2]): 1}, Permutation([3,1,4,2]), GeneratingSet("y"), GeneratingSet("z"), elem_sym))