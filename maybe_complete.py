
from schubmult import *
from schubmult.symbolic import *
from schubmult.abc import *
from schubmult.rings.variables import ZeroGeneratingSet
from schubmult.rings.schub_poly import schubpoly_classical_from_elems
import sympy

def elem_to_complete_expr(p, k, var1, var2):
    if k < 0 or p > k:
        return S.Zero
    if p == 0:
        return S.One
    if p == 1:
        return h(p, k, var1, var2)
    return Add(*[(S.NegativeOne**(i+1))*h(i, k, var1, var2) * elem_to_complete_expr(p - i, k, var1, var2) for i in range(1, p + 1)])


def elem_to_complete_expr2(p, k, var1, var2):
    if k < 0 or p > k:
        return S.Zero
    if p == 0:
        return S.One
    if p == 1:
        return H(p, k, var1, var2)
    return sympy_Add(*[(S.NegativeOne**(i+1))*H(i, k, var1, var2) * elem_to_complete_expr(p - i, k, var1, var2) for i in range(1, p + 1)])

def elem_top(*args):
    #print(f"{args=}")
    global _cache
    if args[0] in _cache:
        return _cache[args[0]](*args[1:])
    return S.Zero

_cache = {
        0: lambda *x: S.One,
        1: lambda *x: H(1, *x),
          }

def elem_full_func(*bob):
    global _cache
    print(f"{bob=}")
    if bob[0] in _cache:
        return _cache[bob[0]](*bob[1:])
    p = bob[0]
    a = (S.NegativeOne**p) * schubpoly_classical_from_elems(uncode(([0]*(k-1))+[p]),x,y,elem_func = elem_top) + (S.NegativeOne**(p+1))* H(p, *bob[1:])
    _cache[p] = lambda *x: elem_full_func(x[0], *x[1:])
    return a

elem_symbs = []

for i in range(20):
    elem_symbs += [[S.One, *[S.Zero if j > i else sympy.Symbol(f"e_{j}_{i}") for j in range(1, 20)]]]

def elem_func(i, j, varl1, varl2):
    if i == 0 and j>=0:
        return S.One
    if i > j or i < 0:
        return S.Zero
    spark = E(i, j, varl1, varl2)
    return sympy.Symbol(f"E({','.join([str(arg) for arg in spark.args])})")

# p, k = 4, 6

# print(schubpoly_classical_from_elems(uncode(([0]*(k-1))+[p]),x,y,elem_func = elem_func))

# print(elem_func(2, 2, x[1:], y[1:]))
# print(sympy.solve(schubpoly_classical_from_elems(uncode(([0]*(k-1))+[4]),x,y,elem_func = elem_func) - sympy.Symbol("H"), elem_func(4, 6, x[1:], y[1:])))
# p = 3
# k = 4
# A = DSx(uncode(([0]*(k-1))+[p])).in_CEM_basis()
# #A = schubpoly_classical_from_elems
# B = H(p, k, x[1:], y[1:])
# C = elem_full_func(p,k,x[1:],y[1:])
# DD = E(p,k,x[1:],y[1:])
# #print(f"{A=}")
# print(f"{C=}")
# #print(f"{C=}")
# #print(f"{expand(A-B, func=True)}")
# print(f"{expand(C-DD,func=True)=}")
# #print(f"{=}")

#complete hom groebner


def double_elem_to_complete_expr(p, k, var1, var2):
    if k < 0 or p > k:
        return S.Zero
    if p == 0:
        return S.One
    if p == 1:
        return h(p, k, var1, var2)
    return sympify_sympy(Add(*[(S.NegativeOne**(i+1))*H(i, k, var1,z[1:]) * elem_to_complete_expr(p - i, k, var1, var2) for i in range(1, p + 1)]))
Permutation.print_as_code = True
p, k = 3, 4
A = E(p, k, x[1:], y[1:])
B = double_elem_to_complete_expr(p, k, x[1:],y[1:])
print(f"{A=}")
print(f"{DSx(expand_func(B))=}")
print(f"{expand(A-B,func=True)}")

