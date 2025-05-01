from schubmult.utils.logging import init_logging
init_logging(False)
from schubmult.poly_lib.elem_sym import pull_out_vars
from schubmult import *
from schubmult.abc import *
from sympy import *


def remove_zeros(poly):
    if expand(poly, func=True) == S.Zero:
        return S.Zero
    if isinstance(poly, Add):
        return Add(*[remove_zeros(arg) for arg in poly.args])
    if isinstance(poly, Mul):
        if any(expand(arg, func=True) == S.Zero for arg in poly.args):
            return S.Zero
        else:
            return poly
    if isinstance(poly, Pow):
        return Pow(remove_zeros(poly.args[0]), poly.args[1])
    return poly

coeffs = (DSx([4,1,5,2,3],elem_sym=True)*DSx([5,1,4,2,3],"z"))
import random 
bob = None
it = iter(coeffs.values())
for i in range(random.randrange(1,len(coeffs.values()))):
    try:
        a = next(it)
    except Exception:
        break
    if not isinstance(a, int) and not isinstance(a, Integer):
        bob = a

bob = remove_zeros(simplify(bob))

print(bob)

varlist = [sympify(arg) for arg in y[1:6]]
varlist2 = [sympify(arg) for arg in z[1:10]]



import random

bob2 = bob

print(bob2)

final = bob2
# for j in range(20):
#     random.shuffle(varlist2)
#     bob3 = split_out_vars(bob2, random.sample(varlist, 3), varlist2)
#     bob3 = expand(remove_zeros(simplify(bob3.expand(mul=True, func=False))),func=False, mul=True)
#     if bob3.count(S.NegativeOne) < bob2.count(S.NegativeOne):
#         bob2 = bob3
#     # try:
#     #     if len(bob3.match(S.NegativeOne)) < len(bob2.match(S.NegativeOne)):
#     #         bob2 = bob3
#     # except TypeError:
#     #     if not bob3.match(S.NegativeOne):
#     #         final = bob3 
#     #     elif not bob2.match(S.NegativeOne):
#     #         final = bob2


# for j in range(20):
#     random.shuffle(varlist2)
#     bob3 = split_out_vars(bob2, random.sample(varlist, 2), varlist2)
#     bob3 = expand(remove_zeros(simplify(bob3.expand(mul=True, func=False))),func=False, mul=True)
#     if bob3.count(S.NegativeOne) < bob2.count(S.NegativeOne):
#         bob2 = bob3
negf = lambda x: isinstance(x, Integer) and x < 0
boing = True
final = None
nomoves = 0
while nomoves < 20:
    res = None
    changed = True
    while changed:
        symbs = [v for v in (res if res else bob2).free_symbols if y.index(v) != -1]
        varlist2 = [v for v in (res if res else bob2).free_symbols if z.index(v) != -1]
        random.shuffle(varlist2)
        v = symbs
        if len(symbs) > 1:
            v = random.sample(symbs, 1 if boing else random.randrange(1,len(symbs)))
        bob3 = expand(remove_zeros(simplify(split_out_vars((res if res else bob2), v, varlist2))))
        if bob3 == (res if res else bob2):
            changed = False
            break
        if not res:
            res = bob3
        if bob3.count(S.NegativeOne) <= 1.7 * res.count(S.NegativeOne):
            res = bob3
        # else:
        #     print(f"{bob3.count(S.NegativeOne)=} {res.count(S.NegativeOne)=}")
    if (res and not final) or (res and res.count(negf) < final.count(negf)):
        final = res
        print(f"{final=} {final.count(negf)=}")
        if final.count(negf) == 0:
            print("Success")
            break
    elif res:
        bob2 = simplify(expand(bob2))
    nomoves += 1




# #print(simplify(expand(split_out_vars(split_out_vars(bob, [y[1]], [z[2],z[1]]),[y[4]],[z[1],z[2]]))))
# bob2 = bob

# for min_degree in range(5, 0, -1):
#     for j in range(1,6):
#         for i in range (1,6):
#             print(f"i: {i}, j: {j}")
#             bob2 = pull_out_vars(bob2, y[i], z[j], min_degree=min_degree)
#             bob2 = expand(bob2, func=False)
#             bob2 = remove_zeros(simplify(bob2))
#             print(f"bob2: {bob2}")
#             print(f"bob-bob2: {expand(bob-bob2, func=True)}")
# #bob2 = simplify(expand(pull_out_vars(pull_out_vars(pull_out_vars(bob, y[1], z[1]),y[2],z[1]),y[4],z[1]), mul=True))
# print(expand(remove_zeros(simplify(bob2)),mul=True,func=False))
# print(remove_zeros(simplify(expand(remove_zeros(simplify(bob2)),mul=False,func=True))))
print(expand(bob-final, func=True))
# print(expand(remove_zeros(simplify(bob2)),func=False,mul=True))
print(final)
print(simplify(expand_func(final)))
print(f"{final.count(negf)=}")
print(f"{final=}")
final2 = None
while not final2 or final2 != final:
    if not final2:
        final2 = final
    final = final2
    final2 = elem_sym_unify(simplify(expand_mul(final2)))
print(f"{final2=}")

# split split out one

# coeff2 = coeff.expand(mul=True, func=False)
# from sympy import sympify
# f = coeff2.args[0].args[1]
# print(f)
# g = f.split_out_vars([y[5]], [z[2],z[3]])
# print(f"{f=}")
# print(f"{g=}")
# import sympy
# print(f"{sympy.expand(f-g, func=True)}")
# from sympy import sympify
# print(getattr(bob, '_eval_simplify', None))
# bob = sympify(bob, rational=True)
# print(getattr(bob, '_eval_simplify', None))
# print(bob._eval_simplify())