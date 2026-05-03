from schubmult import *
from schubmult.symbolic import *
from schubmult.abc import *
from schubmult.utils.schub_lib import has_bruhat_descent
from sympy import Rational
from schubmult.utils.tuple_utils import pad_tuple
import math

def dual_schub(perm1, perm2):
    if perm1 == perm2:
        return 1
    if not perm1.bruhat_leq(perm2):
        return 0
    sm = 0
    perm_stack = [perm2]
    seen = set()
    while perm_stack:
        perm = perm_stack.pop()
        seen.add(perm)
        for i in range(len(perm)):
            for j in range(i + 1, len(perm)):
                if has_bruhat_descent(perm, i, j):
                    test_perm = perm.swap(i, j)
                    if perm1.bruhat_leq(test_perm):
                        
                        #if test_perm not in seen:
                        sm +=  (S.NegativeOne**(test_perm.inv - perm1.inv) ) * (x[i + 1] - x[j + 1]) * dual_schub(perm1, test_perm)
                            #perm_stack.append(test_perm)
    return sm*Rational(1,math.factorial(perm2.inv - perm1.inv))

if __name__ == "__main__":
    import sys
    from schubmult.rings.polynomial_algebra import *
    n = int(sys.argv[1])    
    perms = Permutation.all_permutations(n)
    import sympy
    for perm2 in perms:
        if perm2.inv == 0:# or perm2[0] != 1:
            continue
        ds = dual_schub(Permutation([]), perm2)
        #print(f"dual_schub(1, {perm2}) = {ds}")
        piss = PA.from_expr(ds)
        piss2 = 0
        for k, coeff in piss.items():
            k = pad_tuple(k, n)
            piss2 += coeff * prod([sympy.factorial(a) for a in k])*PA(*k)
        piss2 = piss2 * math.factorial(perm2.inv)
        # print(f"piss2 = {piss2}")
        # print(f"{}")
        the_dual = PA.from_dict(ASx(perm2, n).change_basis(WordBasis))
        # if not all((v  == piss2.get(k, 0)) or (v == -piss2.get(k, 0)) for k, v in the_dual.items()):
        #     print(f"Failed for {perm2}, got {piss2} but expected {the_dual}")
        q, r = sympy.poly(sympy.sympify(piss2.expand()), *[sympy.sympify(a) for a in x[1:n + 1]]).div(sympy.poly(sympy.sympify(the_dual.expand()), *[sympy.sympify(a) for a in x[1:n + 1]]))
        if r != S.Zero:
            print(f"Failed for {perm2}, got {piss2} but expected {the_dual}")
        print(f"Success {perm2} scalar = {q.as_expr()=}")
        if not q.as_expr().is_Rational:
            print(f"Failed for {perm2}, got non-rational coefficient {q.as_expr()} in quotient, expected rational.")