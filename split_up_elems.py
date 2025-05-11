import symengine
import sympy
from sympy import Dict, div

from schubmult import *
from schubmult.abc import *
from schubmult.rings import DoubleSchubertRing, SingleSchubertRing
from schubmult.schub_lib.positivity import compute_positive_rep
from schubmult.symbolic import Add, Integer, Mul, Pow, S, expand, expand_func, prod, sympify
from schubmult.symmetric_polynomials import canonicalize_elem_syms_coeff, coeffvars, degree, genvars, is_of_func_type, numvars
from schubmult.utils.perm_utils import add_perm_dict

r = GeneratingSet("r")
znz = {z[i]: -z[i] for i in range(20)}
subs_dict = {y[1]: S.Zero}

for i in range(1, len(y[1:])):
    subs_dict[y[i + 1]] = subs_dict[y[i]] + r[i]

subs_dict2 = {r[i]: y[i + 1] - y[i] for i in range(1, 20)}
# symengine.var("y_(1:100)")
# symengine.var("z_(1:100)")
# bargain = E(1, 1, y_4, z_1)*E(1, 1, y_4, z_2)*E(1, 3, y_1, y_4, y_5, z_1, z_2, z_3) + E(1, 1, y_4, z_1)*E(1, 1, y_5, z_4)*E(1, 3, y_1, y_4, y_5, z_1, z_2, z_3) - E(1, 1, y_4, z_1)*E(2, 3, y_1, y_4, y_5, z_1, z_2) + E(1, 1, y_5, z_2)*E(1, 1, y_5, z_4)*E(1, 3, y_1, y_4, y_5, z_1, z_2, z_3) - E(1, 1, y_5, z_4)*E(2, 3, y_1, y_4, y_5, z_1, z_2) + E(3, 3, y_1, y_4, y_5, z_1)
sympy.init_printing(pretty_print=False)

z_ring = SingleSchubertRing(z)
zy_ring = DoubleSchubertRing(z, y)

# bagel1 = Permutation([4, 6, 1, 5, 7, 2, 3])
# porn = Permutation([5, 7, 1, 4, 6, 2, 3])
bagel1 = Permutation([4, 1, 5, 2, 3])
porn = Permutation([5, 1, 4, 2, 3])

ctgood = 0
ctbad = 0


def coeff_to_monom(monom, genset):
    dct = {}
    if genset.index(monom) != -1:
        dct[genset.index(monom)] = 1
    if isinstance(monom, Pow):
        dct[genset.index(monom.args[0])] = int(monom.args[1])
    else:
        for arg in monom.args:
            dct = add_perm_dict(dct, coeff_to_monom(arg, genset))
    return Dict(dct)


def genvars_monom(arg):
    if is_of_func_type(arg, FactorialElemSym) or is_of_func_type(arg, FactorialCompleteSym):
        monom = prod(genvars(arg))
    elif isinstance(arg, Mul):
        monom = S.One
        for arg2 in arg.args:
            if is_of_func_type(arg2, FactorialElemSym) or is_of_func_type(arg2, FactorialCompleteSym):
                monom *= prod(genvars(arg2))
            if isinstance(arg2, Pow) and (is_of_func_type(arg2.args[0], FactorialElemSym) or is_of_func_type(arg2.args[0], FactorialCompleteSym)):
                monom *= prod(genvars(arg2.args[0])) ** int(arg2.args[1])
    elif isinstance(arg, Pow):
        monom = prod(genvars(arg.args[0])) ** int(arg.args[1])
    else:
        monom = S.NegativeInfinity
    return monom


def coeffvars_monom(arg):
    if is_of_func_type(arg, FactorialElemSym) or is_of_func_type(arg, FactorialCompleteSym):
        monom = coeffvars(arg)[0] ** degree(arg)
    elif isinstance(arg, Mul):
        monom = S.One
        for arg2 in arg.args:
            if is_of_func_type(arg2, FactorialElemSym) or is_of_func_type(arg2, FactorialCompleteSym):
                monom *= coeffvars(arg2)[0] ** degree(arg2)
            if isinstance(arg2, Pow) and (is_of_func_type(arg2.args[0], FactorialElemSym) or is_of_func_type(arg2.args[0], FactorialCompleteSym)):
                monom *= coeffvars(arg2.args[0])[0] ** (degree(arg2.args[0]) * int(arg2.args[1]))
    elif isinstance(arg, Pow):
        monom = coeffvars(arg.args[0])[0] ** (degree(arg.args[0]) * int(arg.args[1]))
    else:
        monom = S.NegativeInfinity
    return monom


def splitupgenvars(pos_neg_part, comp=False):
    if isinstance(pos_neg_part, Add):
        bacon = pos_neg_part
        args = bacon.args
    else:
        args = [pos_neg_part]
    dct = {}
    for arg in args:
        if comp:
            monom = coeffvars_monom(arg)
        else:
            monom = genvars_monom(arg)
        dct[coeff_to_monom(monom, y)] = dct.get(coeff_to_monom(monom, y), S.Zero) + arg
    return dct


def splitupcoeffvars(pos_neg_part, comp=False):
    if isinstance(pos_neg_part, Add):
        bacon = pos_neg_part
        args = bacon.args
    else:
        args = [pos_neg_part]
    dct = {}
    for arg in args:
        if comp:
            monom = genvars_monom(arg)
        else:
            monom = coeffvars_monom(arg)
        dct[coeff_to_monom(monom, z)] = dct.get(coeff_to_monom(monom, z), S.Zero) + arg
    return dct


# def splitupcoeffvars(pos_neg_part):
#     if isinstance(pos_neg_part, Add):
#         bacon = pos_neg_part
#         args = bacon.args
#     else:
#         args = [pos_neg_part]
#     dct = {}
#     for arg in args:
#         if is_of_func_type(arg, FactorialElemSym):
#             monom = prod(genvars(arg))
#         elif isinstance(arg, Mul):
#             monom = S.One
#             for arg2 in arg.args:
#                 if is_of_func_type(arg2, FactorialElemSym):
#                     monom *= prod(genvars(arg2))
#                 if isinstance(arg2, Pow) and is_of_func_type(arg2.args[0], FactorialElemSym):
#                     monom *= prod(genvars(arg2.args[0]))  ** int(arg2.args[1])
#         elif isinstance(arg, Pow):
#             monom = prod(genvars(arg.args[0])) ** int(arg.args[1])
#         dct[coeff_to_monom(monom,y)] = dct.get(coeff_to_monom(monom,y), S.Zero) + arg
#     return dct
success = 0
fail = 0


def splitupallvars(pos_neg_part):
    if isinstance(pos_neg_part, Add):
        bacon = pos_neg_part
        args = bacon.args
    else:
        args = [pos_neg_part]
    dct = {}
    for arg in args:
        monom1 = Dict(coeff_to_monom(genvars_monom(arg), y))
        monom2 = Dict(coeff_to_monom(coeffvars_monom(arg), z))
        dct[(monom1, monom2)] = dct.get((monom1, monom2), S.Zero) + arg
    return dct


# simplify graph
for k, bargain in (DSx(bagel1, elem_sym=True) * DSx(porn, "z")).items():
    if expand(bargain, func=True) == S.Zero or isinstance(bargain, Integer):
        continue
    bargain = expand(bargain)
    bacon = bargain
    bacon = sympify(FactorialCompleteSym.from_expr_elem_sym(bacon))
    bacon_new = 1
    while bacon != bacon_new:
        bacon_new = bacon
        if isinstance(bacon, Add):
            for arg in bacon.args:
                if not (isinstance(arg, Mul) and isinstance(arg.args[0], Integer) and arg.args[0] < 0):
                    expr = sympy.sympify(arg)
                    if isinstance(arg, Mul):
                        args = arg.args
                    elif isinstance(arg, Pow):
                        args = int(arg.args[1]) * [arg.args[0]]
                    elif isinstance(arg, FactorialCompleteSym):
                        args[arg]
                    else:
                        continue
                    for fromp in args:
                        if isinstance(fromp, Pow):
                            fromp = fromp.args[0]
                        if not isinstance(fromp, FactorialCompleteSym):
                            continue
                        marfle = fromp
                        bacon = expand(split_out_vars(bacon, genvars(marfle)[:1], None))

    if isinstance(bacon, Integer):
        continue
        # break
        # bacon = expand(split_out_vars(bacon, None, coeffvars(marfle)[-1:]))
    # print(f"first {bacon=}")
    # #bacon = canonicalize_elem_syms(bacon)
    # bacon_new = 0
    # #bacon = canonicalize_elem_syms(bacon)
    # # while bacon != bacon_new:
    # #     bacon_new = bacon
    # #     if isinstance(bacon, Add):
    # #         for arg in bacon.args:
    # #             if isinstance(arg, Mul) and isinstance(arg.args[0],Integer) and arg.args[0] < 0:
    # #                 expr = sympy.sympify(arg)
    # #                 for fromp in arg.args[1:]:
    # #                     if isinstance(fromp, Pow):
    # #                         fromp = fromp.args[0]
    # #                     marfle = fromp
    # #                     bacon = expand(split_out_vars(bacon, genvars(marfle)[:len(genvars(marfle))//2], None))
    # #                     #break
    # #                 #bacon = expand(split_out_vars(bacon, None, coeffvars(marfle)[-1:]))
    # # print(f"second {bacon=}")
    # #bacon = canonicalize_elem_syms_coeff(bacon)
    # # if sympy.sstr(sympy.sympify(bacon)).find("-") != -1:
    # #     print("Trying to save it")
    # #     bacon = sympify(FactorialCompleteSym.from_expr_elem_sym(bacon))
    # #     assert (expand(bacon - bargain, func=True) == S.Zero)
    # #     print(f"iter {bacon=}")
    # #     bacon_new = 0
    # #     while bacon != bacon_new:
    # #         bacon_new = bacon
    # #         if isinstance(bacon, Add):
    # #             for arg in bacon.args:
    # #                 if isinstance(arg, Mul) and isinstance(arg.args[0],Integer) and arg.args[0] < 0:
    # #                     expr = arg
    # #                     for fromp in arg.args[1:]:
    # #                         if isinstance(fromp, Pow):
    # #                             fromp = fromp.args[0]
    # #                         marfle = fromp
    # #                         bacon = expand(split_out_vars(bacon, genvars(marfle)[:1], None))
    # #                         #break
    #                     #bacon = expand(split_out_vars(bacon, None, coeffvars(marfle)[-1:]))
    # print(f"{sympy.sympify(bacon)=}")
    # print(f"{expand(bacon - bargain,func=True)}")
    # assert (expand(bacon - bargain, func=True) == S.Zero)
    bacon = FactorialCompleteSym.to_expr_elem_sym(bacon)
    if str(sympy.sympify(bacon)).find("-") == -1:
        success += 1
    else:
        fail += 1
    continue
    dct = splitupcoeffvars(bacon, comp=True)
    if isinstance(bacon, Add):
        #     dct = {}
        #     for arg in bacon.args:
        #         if is_of_func_type(arg, FactorialElemSym):
        #             monom = coeffvars(arg)[0] ** degree(arg)
        #         elif isinstance(arg, Mul):
        #             monom = S.One
        #             for arg2 in arg.args:
        #                 if is_of_func_type(arg2, FactorialElemSym):
        #                     monom *= coeffvars(arg2)[0] ** degree(arg2)
        #                 if isinstance(arg2, Pow) and is_of_func_type(arg2.args[0], FactorialElemSym):
        #                     monom *= coeffvars(arg2.args[0])[0] ** (degree(arg2.args[0]) * int(arg2.args[1]))
        #         elif isinstance(arg, Pow):
        #             monom = coeffvars(arg.args[0])[0] ** (degree(arg.args[0]) * int(arg.args[1]))
        #         dct[monom] = dct.get(monom, S.Zero) + arg
        # dct = splitupallvars(bacon)
        # schuber = z_ring.from_expr(Add(*[expand(efficient_subs(vl,{y[i]: S.Zero for i in range(10)}),func=True,mul=False) for vl in dct.values()]))
        # schuber = z_ring.from_expr(expand(efficient_subs(bargain,znz),func=True,mul=False))
        # if schuber.expand() == S.Zero:
        #     print("Dodged a bullet")
        #     continue
        # if any(isinstance(val,Add) for val in dct.values()):
        #     raise Exception("bongdunket")
        # print(schuber.ring.from_dict({k: abs(v) for k, v in schuber.items()}))
        # dctyep = {}
        # dctnope = {}
        # dctbool = {}
        pos_part = S.Zero
        neg_part = S.Zero
        pos_neg_part = S.Zero
        anyn = False
        dct_boingus = dct
        fownp = dct
        # for key, voig in dct_boingus.items():
        #     dct = splitupgenvars(voig)
        #     for bagel in dct:
        #         fownp[bagel] = fownp.get(bagel,S.Zero) + dct[bagel]
        if True:
            dct = fownp
            for monom, bargle in dct.items():
                if expand(bargle, func=True, deep=True, mul=True) != 0:
                    # voib_bo = efficient_subs(expand(bargle, func=False, mul=True),znz)
                    # voib=expand(voib_bo,func=True,mul=False)
                    # plop = z_ring(voib)
                    # plop2 = z_ring.from_dict({k: expand(efficient_subs(v,subs_dict),mul=False) for k,v in plop.items()})
                    try:
                        compute_positive_rep(expand(bargle, func=True, mul=False), y, z, False, False)
                        # dctbool[monom] = True
                        pos_part += bargle
                        # dctyep[monom] = voib_bo
                    except Exception:
                        flip = True
                        if isinstance(bargle, Add):
                            for arg in bargle.args:
                                if isinstance(arg, Mul) and isinstance(arg.args[0], Integer) and int(arg.args[0]) < 0:
                                    pos_neg_part += arg
                                else:
                                    neg_part += arg
                        else:
                            arg = bargle
                            if isinstance(arg, Mul) and isinstance(arg.args[0], Integer) and int(arg.args[0]) < 0:
                                pos_neg_part += arg
                            else:
                                neg_part += arg
                        anyn = True
                        # dctnope[monom] = voib_bo
            if anyn:
                pos_part = pos_part.expand(deep=False)
                neg_part = neg_part.expand(deep=False)
                pos_neg_part = pos_neg_part.expand(deep=False)
            # print(f"{splitupcoeffvars(pos_part)=}")
            # print(f"{splitupallvars(neg_part)=}")
            # print(f"{splitupallvars(pos_neg_part)=}")
            # dctyep2 = {}
            # dctnope2 = {}

            # def forple(voib2):
            #     plop = zy_ring(efficient_subs(voib2,znz))
            #     plop2 = zy_ring.from_dict({k: expand(efficient_subs(v,subs_dict),mul=False) for k,v in plop.items()})
            #     return plop2
            # #forple = Add(*list(dctyep.values()))
            # for monom, blarp in dctyep.items():
            #     if isinstance(blarp, Add):
            #         dctyep2[monom] = [{k: efficient_subs(v,subs_dict) for k,v in z_ring.from_expr(expand_func(arg)).items() if k in pos_part} for arg in blarp.args]
            #     else:
            #         dctyep2[monom] = [{k: efficient_subs(v, subs_dict) for k,v in z_ring.from_expr(expand_func(blarp)).items() if k in pos_part}]
            # for monom, blarp in dctnope.items():
            #     if isinstance(blarp, Add):
            #         dctnope2[monom] = [{k: efficient_subs(v,subs_dict) for k,v in z_ring.from_expr(expand_func(arg)).items() } for arg in blarp.args]
            #     else:
            #         dctnope2[monom] = [{k: efficient_subs(v, subs_dict) for k,v in z_ring.from_expr(expand_func(blarp)).items()}]
            # print(f"{dctyep=}")
            # print(f"{dctyep2=}")
            # print(f"{dctnope=}")
            # print(f"{dctnope2=}")
            # dct_neg = splitupgenvars(neg_part)
            # dct_pos_neg = splitupgenvars(pos_neg_part)
            # dct_pos = splitupgenvars(pos_part)
            # print(f"{dct_neg=}")
            # print(f"{dct_pos_neg=}")
            # print(f"{dct_pos=}")
            #
            #  dct = {}
            # if isinstance(neg_part, Add):
            #     bacon = neg_part
            #     dct = {}
            #     for arg in bacon.args:
            #         if is_of_func_type(arg, FactorialElemSym):
            #             monom = prod(genvars(arg))
            #         elif isinstance(arg, Mul):
            #             monom = S.One
            #             for arg2 in arg.args:
            #                 if is_of_func_type(arg2, FactorialElemSym):
            #                     monom *= prod(genvars(arg2))
            #                 if isinstance(arg2, Pow) and is_of_func_type(arg2.args[0], FactorialElemSym):
            #                     monom *= prod(genvars(arg2.args[0]))  ** int(arg2.args[1])
            #         elif isinstance(arg, Pow):
            #             monom = prod(genvars(arg.args[0])) ** int(arg.args[1])
            #         dct[monom] = dct.get(monom,S.Zero)
            # else:
            #     dct[prod(genvars(neg_part))] = neg_part
            # pos_part += pos_neg_part
            # dct_remain = {}
            # dct_pos = {}
            # if isinstance(pos_part, Add):
            #     bacon = pos_part
            #     dct = {}
            #     for arg in bacon.args:
            #         if is_of_func_type(arg, FactorialElemSym):
            #             monom = prod(genvars(arg))
            #         elif isinstance(arg, Mul):
            #             monom = S.One
            #             for arg2 in arg.args:
            #                 if is_of_func_type(arg2, FactorialElemSym):
            #                     monom *= prod(genvars(arg2))
            #                 if isinstance(arg2, Pow) and is_of_func_type(arg2.args[0], FactorialElemSym):
            #                     monom *= prod(genvars(arg2.args[0]))  ** int(arg2.args[1])
            #         elif isinstance(arg, Pow):
            #             monom = prod(genvars(arg.args[0])) ** int(arg.args[1])
            #         if monom in dct_neg or monom in dct_pos_neg:
            #             dct_pos[monom] = dct_pos.get(monom, S.Zero) + arg
            #         else:
            #             dct_remain[monom] = dct_remain.get(monom, S.Zero) + arg
            # else:
            #     if isinstance(pos_part, Mul):
            #         args = pos_part.args
            #     elif isinstance(pos_part, Pow):
            #         args = int(pos_part.args[1])*[pos_part.args[0]]
            #     else:
            #         args = [pos_part]

            #     gv = []
            #     for arg in args:
            #         gv += genvars(arg)
            #     #[genvars(arg) for arg in args]
            #     monom = prod(gv)
            #     if monom in dct_neg or monom in dct_pos_neg:
            #         dct_pos[monom] = dct_pos.get(monom, S.Zero) + pos_part
            #     else:
            #         dct_remain[monom] = dct_remain.get(monom, S.Zero) + pos_part
            # print(f"{dct_pos=}")
            # print(f"{dct_remain=}")
            # #print(f"{expand(bargain - neg_part)=}")
            # continue
            # try:
            #     bagel = compute_positive_rep(expand(neg_part, func=True, mul=False), y, z, False, False)
            #     bagel = bagel+expand_func(pos_part)
            #     print(f"Yay: {bagel}")
            #     ctgood += 1
            # except Exception:
            #     ctbad+=1
            #     import traceback
            #     traceback.print_exc()

            # # def ctneg(expr):
            # #     return str(expr).count("-")
            # # yep_list = [(dctyep[monom],z_ring.from_dict(v)) for monom, v in dctyep2.items()]
            # # nope_list = [(dctnope[monom],z_ring.from_dict(v)) for monom, v in dctnope2.items()}]
            # # # total_list = [*yep_list, *nope_]
            # # # fun_list = []
            # # # ct = len(nope_list)
            # # # while ct:
            # # #     new_nope_list = []
            # # #     for spob, elem in nope_list:
            # # #         new_yep_list = []
            # # #         for spob2, elem2 in nope_list:
            # # #             if ctneg(elem + elem2) < ctneg(elem) + ctneg(elem2):
            # # #                 new_yep_list

    # Anything that isn't Schub will disappear
