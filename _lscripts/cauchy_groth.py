from schubmult import *
from schubmult.rings.schubert.grothendieck_ring import GrothendieckRing
from schubmult.abc import x, y
from schubmult.symbolic.common_polys import grothendieck_poly
from schubmult.symbolic import S, is_of_func_type
from functools import cache
from schubmult.symbolic.symmetric_polynomials import CompleteSym_base, ElemSym_base, FactorialElemSym, coeffvars, degree, genvars, numvars, split_out_vars

#$$F_k(x;y) = \sum_{i=0}^{k} f_{k-i}(x)(y_1^i + y_1^{i+1})$$

def _elem_perm(p, k):
    return uncode([0] * (k - p) + [1] * p)

rx = GrothendieckRing(x, beta=1)
ry = GrothendieckRing(y, beta=1)

@cache
def _elem_sym(p, k):
    return rx(_elem_perm(p, k)).expand()

@cache
def _elem_sym_s(p, k):
    return Sx(_elem_perm(p, k)).expand()


def elem_sym_by_cauchy(p, k):
    from schubmult.abc import x, y, z, H
    from schubmult.symbolic.poly.schub_poly import elem_sym_poly
    #from sympy import expand, expand_func, sympify, prod, Add, Mul
    from symengine import Add, Mul, Pow
    # sm = 0
    space_set_potato = [x[i]/(1 + x[i]) for i in range(1, 10)]
    ny = [-y[i] for i in range(1, 10)]
    
    ones = [S.NegativeOne for i in range(1, 10)]

    E = FactorialElemSym
    if p > k or k < 0:
        return S.Zero
    if p == 0:
        return S.One
    if p == k:
        return E(k, k, x[1:], ones) * E(k, k, space_set_potato, ny)
    
    # space_set = [z[i] for i in range(1, 10)]
    # nx = [-x[i] for i in range(1, 10)]
    
    fat_donkey = elem_sym_by_cauchy(p + 1, k)

    

    def _recurse(poly):
        if isinstance(poly, Add):
            return sum([_recurse(arg) for arg in poly.args])
        elif isinstance(poly, Mul):
            return prod([_recurse(arg) for arg in poly.args])
        elif isinstance(poly, Pow):
            return _recurse(poly.base) ** poly. exp
        elif is_of_func_type(poly, ElemSym_base) and len(poly.coeffvars) == p:# and y.index(-poly.coeffvars[0]) != -1:
            print("APSGJK", poly, poly.degree, poly.numvars, poly.coeffvars, poly.genvars)
            return (1 + y[p]) * E(poly.degree - 1, poly.numvars, space_set_potato, ny)
        else:
            return poly

    result = _recurse(fat_donkey) - fat_donkey
    return result
        



def elem_sym_by_cauchy2(k):
    sm = 0
    # sm = rx(_elem_perm(k, k)).expand()
    #z = 1 / y[i]
    for i in range(k + 1):
        sm += (1 + y[1]) * _elem_sym(k - i, k) * (y[1]/(1 + y[1]))

    return sm


#(1 + x_1)d_1 - 1

if __name__ == "__main__":
    #from symengine import expand
    #from sympy import expand, sympify
    from schubmult.symbolic import expand_func
    import sys

    n = int(sys.argv[1])
    for k in range(1, n):
        for p in range(1, k + 1):
            poly1 = elem_sym_by_cauchy(p, k)
            poly2 = grothendieck_poly(_elem_perm(p, k), x, y, beta=1).expand()
            try:
                test_poly = expand_func(poly1).simplify()
                assert (test_poly - poly2).expand() == 0, f"Failed for k={k} {p=}: {test_poly=}, {poly2=}"
                print(f"Success for k={k} {p=}")
            except AssertionError as e:
                print(e)
            