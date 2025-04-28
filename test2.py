from schubmult.schub_lib.schub_lib import compute_vpathdicts, elem_sym_perms_op, elem_sym_perms
from schubmult.perm_lib import *
from schubmult.utils.perm_utils import *
# from schubmult import QPDSx, QPSx
# from schubmult import uncode
# from symengine import S, expand
# # QPDSx = make_parabolic_quantum_basis([2, 3, 4])
# A = QPDSx(2,3,4)(uncode([1,2,0,2,3]), 0)
# #A2 = QPSx(2,3,4)(uncode([1,2,0,2,3]))
# #assert expand(A.as_polynomial() - A2.as_polynomial()) == S.Zero
# B = QPDSx(2,3,4)(uncode([1,3,0,1,2]), 0)
# #B2 = QPSx(2,3,4)(uncode([1,3,0,1,2]))
# #assert expand(B.as_polynomial() - B2.as_polynomial()) == S.Zero
# C = A*B
# #C2 = A2*B2
# print(f"{C=}")
# #print(f"{C2=}")

# residue = expand(C.as_polynomial()-A.as_polynomial()*B.as_polynomial())
# print(f"{residue=}")
# print(f"{expand(C.as_polynomial())=}")
# print(f"{expand(A.as_polynomial())=}")
# print(f"{expand(B.as_polynomial())=}")

# Single CEM iff avoids 1423 and 1432

from schubmult import *
from itertools import permutations
import sys
import sympy

n = int(sys.argv[1])
import symengine
perms = [Permutation(perm) for perm in list(permutations(list(range(1,n+1))))]
# ct = 0
# for perm in perms:
#     pl = Sx(perm)
#     xp = sympy.expand(pl.in_CEM_basis())
#     # cond1 = not perm.has_pattern([1,4,2,3]) and not perm.has_pattern([4,1,3,2]) and not perm.has_pattern([3,1,4,2]) and not perm.has_pattern([1,4,3,2])
#     cond2 = isinstance(xp, sympy.Add)
#     # if (cond1 and cond2) or (not cond1 and not cond2):
#     #     print(f"Fail {perm} {cond1} {cond2} {xp=}")
#     if not cond2:
#         ct += 1
# print(ct)

def is_dominant(p):
    for i in range(len(p)-2):
        if p.code[i] < p.code[i+1]:
            return False
    return True

def check(perm):
    return not perm.has_pattern([1,4,2,3]) and not perm.has_pattern([4,1,3,2]) and not perm.has_pattern([3,1,4,2]) and not perm.has_pattern([1,4,3,2])

def checkone(perm):
    return not perm.has_pattern([1,4,2,3]) or not perm.has_pattern([4,1,3,2]) or not perm.has_pattern([3,1,4,2]) or not perm.has_pattern([1,4,3,2])


def has_4132(perm,four,one,three,two):
    if perm[one] < perm[two] and perm[two] < perm[three] and perm[three] < perm[four]:
        return True
    return False

def has_3142(perm,three,one,four,two):
    if perm[one] < perm[two] and perm[two] < perm[three] and perm[three] < perm[four]:
        return True
    return False

ct_dict = {}
p_dict = {}

def check2(perm):
    # increase by at most 1
    decrease_by_3 = 0
    decrease_by_2 = 0
    min_value = None

    for i in range(len(perm.code) - 1):
        if not min_value or perm.code[i] < min_value:
            min_value = perm.code[i]
        if perm.code[i] - min_value > 1:
            return False
        if perm.code[i] + 1 < perm.code[i+1]:
            return False
        if decrease_by_2 > 0: 
            if perm.code[i] < perm.code[i + 1]:
                return False
            else:
                decrease_by_2 -= 1
        if decrease_by_3 > 0:
            if perm.code[i] < perm.code[i+1]:
                return False
            else:
                decrease_by_3 -= 1
        
        if perm.code[i] - perm.code[i + 1] >= 2:
            decrease_by_2 = 1
        if perm.code[i] - perm.code[i + 1] >= 3:
            decrease_by_3 = 2
    return True

bongo = 0
if __name__ == "__main__":
    # ct_fact = [0 for i in range(n)]
    # for perm in perms:
    #     # found_one = True
    #     # cond1 = check(perm)
    #     # cond2 = check2(perm)

    #     # if cond1 != cond2:
    #     #     print(f"Fail {perm=} {perm.code=} {cond1=} {cond2=}")
    #     # if perm.inv > 0 and check(perm) and not is_dominant(perm):
    #     if perm.inv == 0:
    #         continue
    #     th = theta(~perm)
    #     mu = uncode(th)
    #     th = trimcode(mu)
    #     #     th = theta(~perm)
    #     #     tmu = uncode(th)
    #     #     crud = [b - a for a,b in zip(trimcode(perm),trimcode(~tmu))]
    #     #     fonf = uncode(crud)
    #     #     print(f"{code(perm)} {code(~tmu)} {code(fonf)}")
    #         #print(f"{}")
    #     vmu = perm * mu
    #     vpathdicts = compute_vpathdicts(th, vmu)
        
    #     if not check(perm) and isinstance(Sx(perm).in_CEM_basis().simplify(),symengine.Mul):
    #         print(f"{perm} {trimcode(perm)}")
    #         print(Sx(perm).in_CEM_basis().simplify())
    #         for i in range(len(vpathdicts)):
    #             print(f"{i}: {vpathdicts[i]}")
    # print(f"{n=}")
    # print(f"{ct_fact=}")
    # for i in range(len(perms)):
    #     u = perms[i]
    #     if not u.has_pattern([3,1,4,2]) and not u.has_pattern([4,1,3,2]) and not u.has_pattern([1,4,2,3]) and u.inv > 0:
    #         print(f"{u} {u.code} {Sx(u).in_CEM_basis().simplify()}")
    #         assert not isinstance(Sx(u).in_CEM_basis().simplify(), symengine.Add)
        # # vmu = u*((~u).minimal_dominant_above())
        # # cd = trimcode(~(vmu))
        # # if check(u) and len(cd)>0:
        # #     #print(f"{cd} {trimcode(((~u).minimal_dominant_above()))}")
        # #     #print(f"{u.code}")
        # cd1 = [*u.code]
        # for i in range(1, len(cd1)):
        #     if cd1[i] > cd1[i-1]:
        #         cd1[i] -= 1
        # if check(u) and not is_dominant(u):
        #     print(f"{u.code} {trimcode(uncode(cd1))}")
        #     if not is_dominant(uncode(cd1)):
        #         print(f"Oh no")
        #         exit(1)
                
        
    #     if u.inv == 0:
    #         continue
    #     spock = Sx(u).in_CEM_basis()._sympy_().as_poly().factor_list() 
    #     if len(spock[1]) > 1:
    #         good = True
    #         flist = [Sx((v**e).as_expr().subs(Sx.elem_sym_subs(n))) for v,e in spock[1]]
    #         for flop in flist:
    #             if len(flop.coeff_dict.keys()) > 1:
    #                 good = False
    #                 break
    #         if not good:
    #             print(f"Blang {u.code} {flist}")
    #         else:
    #             bongo += 1
                
            
    #     # for j in range(i, len(perms)):
    #     #     v = perms[j]
    #     #     if v.inv == 0:
    #     #         continue
    #     #     P = Sx(u) * Sx(v)
    #     #     if len(P.coeff_dict.keys()) == 1 and (not is_dominant(u) or not is_dominant(v)):
    #     #         w = next(iter(P.coeff_dict.keys()))[0]
    #     #         print(f"{trimcode(u),trimcode(v)} = {trimcode(w)} {check(u),check(v)}")
    # print(bongo)
    # ct = 0
    # for perm in perms:
    #     # if not check(perm):
    #     #     for i in range(len(perm)-1,3,-1):
    #     #         if check(Permutation(perm[:i])) and not perm.has_pattern([3,1,4,2]) and not perm.has_pattern([4,1,3,2]) and not perm.has_pattern([1,4,3,2]):
    #     #             print(f"{i} {perm} {perm.pattern_at(*list(range(i)))} {Sx(perm).in_CEM_basis().simplify()}")
    #     #             break
    #     # if check(perm):
    #     #     print(f"{perm} {elem_sym_perms_op(perm,len(perm)-2,len(perm)-2)} {Sx(perm).in_CEM_basis()}")
    #     #if not perm.has_pattern([3,1,4,2]) and not perm.has_pattern([4,1,3,2]) and not perm.has_pattern([1,4,3,2])
    #     if not perm.has_pattern([3,1,4,2]):
    #         print(f"{perm} {Sx(perm).in_CEM_basis()}")
    total = 0
    for perm in perms:
        if perm.inv < 8 or not perm.has_pattern([1,3,2]):
            continue
        lastvar = max((~perm).descents()) + 1
        if lastvar == 1:
            continue
        spoly = DSx(perm)
        total += 1
        for i in range(lastvar-1,1,-1):
            try:
                if len(spoly.coproduct(*list(range(i,lastvar+1)),on_coeff_gens=True).keys()) == 1:
                    bongo+=1
                    print(f"{DSx(perm)} = {spoly.coproduct(*list(range(i,lastvar+1)),on_coeff_gens=True)}")
                    break
            except:
                pass
            
    # dct = [{} for i in range(n-1)]

    # for perm in perms:
    #     for i in range(1,n):
    #         esp = elem_sym_perms(perm, i, i)
    #         cts = {}
    #         for perm2, f in esp:
    #             if f > 0:
    #                 cts[f] = cts.get(f,[]) + [perm2]
    #         for f in cts:
    #             if len(cts[f]) == 1:
    #                 dct[i-1][perm] = dct[i-1].get(perm,[]) + [cts[f][0]]
    
    # dct2 = [{} for i in range(n-1)]
    # for i in range(len(dct)):
    #     for perm in dct[i]:
    #         for perm1 in dct[i][perm]:
    #             dct2[i][perm1] = perm

    # for perm in perms:
    #     perm0 = perm
    #     factorization = []
    #     while perm0 != Permutation([]):
    #         found = False
    #         for i in range(len(dct2)):
    #             if perm0 in dct2[i]:
    #                 perm1 = dct2[i][perm0]
    #                 factorization += [(perm0.inv-perm1.inv, i+1)]
    #                 found = True
    #                 perm0 = perm1
    #                 break
    #         if not found:
    #             factorization += [perm0]
    #             break
    #     print(f"{perm} {factorization} {Sx(perm).in_CEM_basis()}")
    print(f"{bongo} {total}")