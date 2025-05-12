
from schubmult.symbolic import *
from schubmult.abc import *

def elem_to_complete_expr(p, k, var1, var2):
    if k < 0 or p > k:
        return S.Zero
    if p == 0:
        return S.One
    if p == 1:
        return h(p, k, var1, var2)
    return Add(*[(S.NegativeOne**(i+1))*h(i, k, var1, var2) * elem_to_complete_expr(p - i, k, var1, var2) for i in range(1, p + 1)])

A = elem_to_complete_expr(2, 3, x[1:], y[1:])
B = e(2, 3, x[1:], y[1:])

print(f"{A=}")
print(f"{B=}")
print(f"{expand(A-B, func=True)}")