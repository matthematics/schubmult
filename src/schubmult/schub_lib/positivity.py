from bisect import bisect_left
from functools import cache

import numpy as np
import psutil
import pulp as pu
from cachetools import cached
from cachetools.keys import hashkey
from sortedcontainers import SortedList
from sympy import Dict

import schubmult.rings as rings
from schubmult.perm_lib import (
    Permutation,
    code,
    cycle,
    dominates,
    inv,
    one_dominates,
    phi1,
    theta,
    uncode,
)
from schubmult.rings.poly_lib import _vars, efficient_subs, elem_sym_poly
from schubmult.rings.schub_poly import schubpoly
from schubmult.schub_lib.schub_lib import (
    divdiffable,
    is_coeff_irreducible,
    is_split_two,
    pull_out_var,
    reduce_coeff,
    reduce_descents,
    try_reduce_u,
    try_reduce_v,
    will_formula_work,
)
from schubmult.symbolic import Add, Integer, Mul, Pow, S, Symbol, expand, is_of_func_type, poly, prod, sympify, sympify_sympy
from schubmult.symmetric_polynomials import FactorialCompleteSym, FactorialElemSym, canonicalize_elem_syms, coeffvars, degree, genvars, numvars
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from .double import schubmult_double, schubmult_double_pair, schubmult_double_pair_generic

logger = get_logger(__name__)


def compute_positive_rep(val, var2=None, var3=None, msg=False, do_pos_neg=True):
    do_pos_neg = False  # noqa: F841
    notint = False
    try:
        val2 = int(expand(val))
        # val2 = expand(val)
    except Exception:
        notint = True
    if notint:
        z_ring = rings.SingleSchubertRing(var3)
        opt = Optimizer(z_ring, val)
        frees = val.free_symbols
        # logger.debug(f"{frees=}")
        # logger.debug(f"{[type(s) for s in frees]=}")
        varsimp2 = [m for m in frees if var2.index(m) != -1]
        varsimp3 = [m for m in frees if var3.index(m) != -1]
        varsimp2.sort(key=lambda k: var2.index(k))
        varsimp3.sort(key=lambda k: var3.index(k))
        # logger.debug(f"{varsimp2=}")
        # logger.debug(f"{varsimp3=}")
        var22 = [sympify_sympy(v) for v in varsimp2]
        var33 = [sympify_sympy(v) for v in varsimp3]
        # var22 = [sympify(m) for m in varsimp2]
        # var33 = [sympify(m) for m in varsimp3]
        n1 = len(varsimp2)

        # for i in range(len(varsimp2)):
        #     varsimp2[i] = var2[var2list.index(varsimp2[i])]
        # for i in range(len(varsimp3)):
        #     varsimp3[i] = var3[var3list.index(varsimp3[i])]

        base_vectors = []
        base_monoms = []
        vec = opt.vec0

        # if do_pos_neg:
        #     smp = val
        #     flat, found_one = flatten_factors(smp)
        #     while found_one:
        #         flat, found_one = flatten_factors(flat, varsimp2, varsimp3)
        #     pos_part = 0
        #     neg_part = 0
        #     if isinstance(flat, Add) and not is_flat_term(flat):
        #         for arg in flat.args:
        #             if expand(arg) == 0:
        #                 continue
        #             if not is_negative(arg):
        #                 pos_part += arg
        #             else:
        #                 neg_part -= arg
        #         if neg_part == 0:
        #             # print("no neg")
        #             return pos_part
        #     depth = 1

        #     mons = split_monoms(pos_part, varsimp2, varsimp3)
        #     mons = {tuple(mn) for mn in mons}
        #     mons2 = split_monoms(neg_part, varsimp2, varsimp3)
        #     mons2 = {tuple(mn2) for mn2 in mons2}

        #     # mons2 = split_monoms(neg_part)
        #     # for mn in mons2:
        #     # if mn not in mons:
        #     # mons.add(mn)
        #     # print(mons)
        #     status = 0
        #     size = len(mons)
        #     while status != 1:
        #         base_monoms, mons = find_base_vectors(mons, mons2, varsimp2, varsimp3, depth)
        #         if len(mons) == size:
        #             raise ValueError("Found counterexample")

        #         size = len(mons)
        #         base_vectors = []
        #         for i in range(len(base_monoms)):
        #             vec0 = opt.poly_to_vec(base_monoms[i])
        #             if vec0 is not None:
        #                 base_vectors += [vec0]

        #         vrs = [pu.LpVariable(name=f"a{i}", lowBound=0, cat="Integer") for i in range(len(base_vectors))]
        #         lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
        #         lp_prob += 0
        #         eqs = opt.base_vec
        #         for j in range(len(base_vectors)):
        #             for i in base_vectors[j]:
        #                 bvi = base_vectors[j][i]
        #                 if bvi == 1:
        #                     eqs[i] += vrs[j]
        #                 else:
        #                     eqs[i] += bvi * vrs[j]
        #         for i in eqs.keys():
        #             lp_prob += eqs[i] == vec[i]
        #         try:
        #             solver = pu.PULP_CBC_CMD(msg=msg)
        #             status = lp_prob.solve(solver)
        #         except KeyboardInterrupt:
        #             current_process = psutil.Process()
        #             children = current_process.children(recursive=True)
        #             for child in children:
        #                 child_process = psutil.Process(child.pid)
        #                 child_process.terminate()
        #                 child_process.kill()
        #             raise KeyboardInterrupt()
        #         status = lp_prob.status
        # else:
        # logger.debug("this")
        val_poly = poly(expand(val), *var22, *var33)
        vec = opt.poly_to_vec(val)
        mn = val_poly.monoms()
        L1 = tuple([0 for i in range(n1)])
        mn1L = []
        lookup = {}
        # logger.debug("this")
        for mm0 in mn:
            key = mm0[n1:]
            if key not in lookup:
                lookup[key] = []
            mm0n1 = mm0[:n1]
            st = set(mm0n1)
            if len(st.intersection({0, 1})) == len(st) and 1 in st:
                lookup[key] += [mm0]
            if mm0n1 == L1:
                mn1L += [mm0]
        # logger.debug("this")
        for mn1 in mn1L:
            comblistmn1 = [1]
            for i in range(n1, len(mn1)):
                if mn1[i] != 0:
                    arr = np.array(comblistmn1)
                    comblistmn12 = []
                    mn1_2 = (*mn1[n1:i], 0, *mn1[i + 1 :])
                    for mm0 in lookup[mn1_2]:
                        comblistmn12 += (
                            arr
                            * np.prod(
                                [varsimp2[k] - varsimp3[i - n1] for k in range(n1) if mm0[k] == 1],
                            )
                        ).tolist()
                    comblistmn1 = comblistmn12
            for i in range(len(comblistmn1)):
                b1 = comblistmn1[i]
                vec0 = opt.poly_to_vec(b1)
                if vec0:
                    base_vectors += [vec0]
                    base_monoms += [b1]
        vrs = [pu.LpVariable(name=f"a{i}", lowBound=0, cat="Integer") for i in range(len(base_vectors))]
        lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
        lp_prob += 0
        eqs = {}
        for j in range(len(base_vectors)):
            for i in base_vectors[j]:
                bvi = int(base_vectors[j][i])
                if bvi == 1:
                    if i not in eqs:
                        eqs[i] = vrs[j]
                    else:
                        eqs[i] += vrs[j]
                elif bvi != 0:
                    if i not in eqs:
                        eqs[i] = bvi * vrs[j]
                    else:
                        eqs[i] += bvi * vrs[j]
        for i in eqs:
            try:
                lp_prob += eqs[i] == vec.get(i, 0)
            except KeyError:
                # print(f"{vec=} {val=}")
                raise
        # print(f"{vec=}")
        # print(lp_prob.constraints)
        try:
            # logger.debug("I IS SOLVING BOLVING")
            solver = pu.PULP_CBC_CMD(msg=msg)
            status = lp_prob.solve(solver)  # noqa: F841
        except KeyboardInterrupt:
            current_process = psutil.Process()
            children = current_process.children(recursive=True)
            for child in children:
                child_process = psutil.Process(child.pid)
                child_process.terminate()
                child_process.kill()
            raise KeyboardInterrupt()
        # print(f"{pos_part=}")
        # print(f"{neg_part=}")
        # else:
        # print(f"No dice {flat=}")
        # exit(1)
        # #val = pos_part - neg_part

        # depth+=1
        val2 = 0
        for k in vrs:
            x = vrs[k].value()
            if x != 0 and x is not None:
                val2 += int(x) * k
        if expand(val - val2) != 0:
            # print(f"{val=}")
            # print(f"{val2=}")
            # print(f"{vec=}")
            raise Exception
    # print(f"{val2=}")
    return val2


def coeffvars_monom(arg, genset):
    # actually put in the vars
    monom = {}
    if is_of_func_type(arg, FactorialElemSym) or is_of_func_type(arg, FactorialCompleteSym):
        monom[genset.index(coeffvars(arg)[0])] = [*genvars(arg)]
    elif isinstance(arg, Mul):
        for arg2 in arg.args:
            if is_of_func_type(arg2, FactorialElemSym) or is_of_func_type(arg2, FactorialCompleteSym):
                i = genset.index(coeffvars(arg2)[0])
                monom[i] = monom.get(i,[]) + [*genvars(arg2)]
            if isinstance(arg2, Pow) and (is_of_func_type(arg2.args[0], FactorialElemSym) or is_of_func_type(arg2.args[0], FactorialCompleteSym)):
                i = genset.index(coeffvars(arg2.args[0])[0])
                monom[i] = monom.get(i,[])  + [*genvars(arg2.args[0])] * int(arg2.args[1])
    elif isinstance(arg, Pow):
        i = genset.index(coeffvars(arg.args[0])[0])
        monom[i] = monom.get(i,[])  + [*genvars(arg.args[0])] * int(arg.args[1])
    return Dict({k: tuple(v) for k, v in monom.items()})


def coeff_to_monom(monom, genset):
    dct = {}
    if monom == S.One:
        return Dict(dct)
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


def splitupcoeffvars(pos_neg_part, genset):
    if isinstance(pos_neg_part, Add):
        bacon = pos_neg_part
        args = bacon.args
    else:
        args = [pos_neg_part]
    dct = {}
    for arg in args:
        dct[coeffvars_monom(arg, genset)] = arg
    return dct


def splitupgenvars(pos_neg_part, genset):
    if isinstance(pos_neg_part, Add):
        bacon = pos_neg_part
        args = bacon.args
    else:
        args = [pos_neg_part]
    dct = {}
    for arg in args:
        monom = genvars_monom(arg)
        dct[coeff_to_monom(monom, genset)] = dct.get(coeff_to_monom(monom, genset), S.Zero) + arg
    return dct

def splitupallvars(pos_neg_part, gs1, gs2):
    if isinstance(pos_neg_part, Add):
        bacon = pos_neg_part
        args = bacon.args
    else:
        args = [pos_neg_part]
    dct = {}
    for arg in args:
        monom1 = Dict(coeff_to_monom(genvars_monom(arg, gs1),gs1))
        monom2 = Dict(coeff_to_monom(coeffvars_monom(arg, gs2),gs2))
        dct[(monom1, monom2)] = dct.get((monom1, monom2), S.Zero) + arg
    return dct


def compute_positive_rep_new(val, var2=None, var3=None, msg=False, do_pos_neg=True):
    val_expr = expand(val, func=True)
    z_ring = rings.SingleSchubertRing(var3)
    opt = Optimizer(z_ring,val_expr)
    vec_base = opt.vec0
    val = canonicalize_elem_syms(expand(val))
    mnset = splitupcoeffvars(val, var3)
    #mndct = {k: set(splitupgenvars(v,var2).keys()) for k, v in dctcv.items()}

    # mndct = {}
    # for k, st in dctmn.items():
    #     if k not in mndct:
    #         mndct[k] = set()
    #     mndct[k].update(
    #combcache={}
    base_vectors = {}
    base_monoms = []
    #base_monoms += [b1]s
    # lookup = {}
    frees = val_expr.free_symbols
        # logger.debug(f"{frees=}")
        # logger.debug(f"{[type(s) for s in frees]=}")
    varsimp2 = [m for m in frees if var2.index(m) != -1]
    varsimp3 = [m for m in frees if var3.index(m) != -1]
    varsimp2.sort(key=lambda k: var2.index(k))
    varsimp3.sort(key=lambda k: var3.index(k))
    # logger.debug(f"{varsimp2=}")
    # logger.debug(f"{varsimp3=}")
    var22 = [sympify_sympy(v) for v in varsimp2]
    var33 = [sympify_sympy(v) for v in varsimp3]
    # var22 = [sympify(m) for m in varsimp2]
    # var33 = [sympify(m) for m in varsimp3]
    n1 = len(varsimp2)
    val_poly = S.Zero
    for k, expr in mnset.items():
        val_poly += poly(expand(expr, func=True),*var22, *var33)
    mn = val_poly.monoms()
    L1 = tuple([0 for i in range(n1)])
    mn1L = []
    lookup = {}
    # logger.debug("this")
    for mm0 in mn:
        key = mm0[n1:]
        if key not in lookup:
            lookup[key] = []
        mm0n1 = mm0[:n1]
        st = set(mm0n1)
        if len(st.intersection({0, 1})) == len(st) and 1 in st:
            lookup[key] += [mm0]
        if mm0n1 == L1:
            mn1L += [mm0]
    # logger.debug("this")
    for mn1 in mn1L:
        comblistmn1 = [1]
        for i in range(n1, len(mn1)):
            if mn1[i] != 0:
                arr = np.array(comblistmn1)
                comblistmn12 = []
                mn1_2 = (*mn1[n1:i], 0, *mn1[i + 1 :])
                for mm0 in lookup[mn1_2]:
                    comblistmn12 += (
                        arr
                        * np.prod(
                            [varsimp2[k] - varsimp3[i - n1] for k in range(n1) if mm0[k] == 1],
                        )
                    ).tolist()
                comblistmn1 = comblistmn12
        for i in range(len(comblistmn1)):
            b1 = comblistmn1[i]
            vec0 = opt.poly_to_vec(b1)
            if vec0:
                base_vectors[b1] = vec0
        # logger.debug("this")
    # mnsetset = {vv: {k: set(voing) for k, voing in vv.items()} for vv in mnset}
    # for vv in mnset:
    #     key = vv
    #     if key not in lookup:
    #         lookup[key] = {}
    #     vvset = mnsetset[vv]
    #     bad = False
    #     poinkset = set()
    #     for z_index, pw in vvset.items():
    #         lookup[key][z_index] = []
    #         for vv2 in mnset:
    #             vv2set = mnsetset[vv2]
    #             bad = True
    #             if vv2set.get(z_index,set()).issubset(pw):
    #                 for pork, bingo in vv2set.items():
    #                     if pork == z_index:
    #                         continue
    #                     if pork in vv and not vvset[pork].issubset(bingo):
    #                         bad = True
    #                         break
    #                     bad = False
    #                     poinkset.update(bingo)
    #             if not bad:
    #                 lookup[key][z_index] += [frozenset(poinkset)]
    # for mn1 in mnset:
    #     comblistmn1 = [S.One]
    #     for z_index, vs in mn1.items():
    #         arr = np.array(comblistmn1)
    #         comblistmn12 = []
    #         #n1 = mn1[z_index]
    #         #combinations degree in the lists
    #         # fslower but we can do better
            
    #         # print(lst)
    #         #cached combs
    #         from itertools import combinations
    #         for lst in lookup[mn1][z_index]:
    #             if (lst, len(vs)) not in combcache:
    #                 combcache[(lst, len(vs))] = list(combinations(lst, len(vs)))
    #             combs = combcache[(lst,len(vs))] 
    #             for comb in combs:
    #                 comblistmn12 += (
    #                     arr
    #                     * np.prod(
    #                         [k - var3[int(z_index)] for k in comb],
    #                     )
    #                 ).tolist()
    #         comblistmn1 = comblistmn12
    #     # print(comblistmn12)
    #     for i in range(len(comblistmn1)):
    #         b1 = comblistmn1[i]
    #         vec0 = opt.poly_to_vec(b1)
            
    #         # if b1 in base_monoms:
    #         #     continue
    #         if vec0 is not None:
    #             base_vectors[b1] = vec0
    vrs = {bv: pu.LpVariable(name=f"a{bv}", lowBound=0, cat="Integer") for bv in base_vectors}
    lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
    lp_prob += 0
    eqs = {}
    for bv, vec in base_vectors.items():
        for i in vec:
            bvi = int(vec[i])
            if bvi == 1:
                if i not in eqs:
                    eqs[i] = vrs[bv]
                else:
                    eqs[i] += vrs[bv]
            elif bvi != 0:
                if i not in eqs:
                    eqs[i] = bvi * vrs[bv]
                else:
                    eqs[i] += bvi * vrs[bv]
    for i in eqs:
        try:
            lp_prob += eqs[i] == vec_base[i]
        except KeyError:
            raise
    # print(f"{vec=}")
    # print(lp_prob.constraints)
    try:
        # logger.debug("I IS SOLVING BOLVING")
        solver = pu.PULP_CBC_CMD(msg=msg)
        status = lp_prob.solve(solver)  # noqa: F841
    except KeyboardInterrupt:
        current_process = psutil.Process()
        children = current_process.children(recursive=True)
        for child in children:
            child_process = psutil.Process(child.pid)
            child_process.terminate()
            child_process.kill()
        raise KeyboardInterrupt()
    # print(f"{pos_part=}")
    # print(f"{neg_part=}")
    # else:
    # print(f"No dice {flat=}")
    # exit(1)
    # #val = pos_part - neg_part

    # depth+=1
    val2 = 0
    for k in base_vectors:
        x = vrs[k].value()
        if x != 0 and x is not None:
            val2 += int(x) * k
    if expand(val - val2, func=True) != 0:
        # print(f"{vec=}")
        raise Exception
    # print(f"{val2=}")
    return val2


@cached(
    cache={},
    key=lambda val, u2, v2, w2, var2=None, var3=None, msg=False, do_pos_neg=False, sign_only=False, optimize=True, elem_sym=None: hashkey(val, u2, v2, w2, var2, var3, msg, do_pos_neg, sign_only, optimize),
)
def posify(
    val,
    u2,
    v2,
    w2,
    var2=None,
    var3=None,
    msg=False,
    do_pos_neg=False,
    sign_only=False,
    optimize=True,
    n=_vars.n,
    elem_sym=None,
):
    if not v2.has_pattern([1, 4, 2, 3]) and not v2.has_pattern([4, 1, 3, 2]) and not v2.has_pattern([3, 1, 4, 2]) and not v2.has_pattern([1, 4, 3, 2]):
        logger.debug("Recording new characterization was used")
        return schubmult_double({u2: 1}, v2, var2, var3).get(w2, 0)
    # logger.debug(f"NEW {val=} {u2=} {v2=} {w2=}")
    oldval = val
    if inv(u2) + inv(v2) - inv(w2) == 0:
        # logger.debug(f"Hmm this is probably not or val inty true {val=}")
        return val

    if not sign_only and expand(val) == 0:
        # logger.debug(f"Hmm this is probably not true {u2=} {v2=} {w2=} {val=}")
        return 0
    # logger.debug("proceeding")
    u, v, w = u2, v2, w2
    # u, v, w = try_reduce_v(u2, v2, w2)
    if is_coeff_irreducible(u2, v2, w2):
        u, v, w = try_reduce_u(u2, v2, w2)
        if is_coeff_irreducible(u, v, w):
            u, v, w = u2, v2, w2
            if is_coeff_irreducible(u, v, w):
                w0 = w
                u, v, w = reduce_descents(u, v, w)
                if is_coeff_irreducible(u, v, w):
                    u, v, w = reduce_coeff(u, v, w)
                    if is_coeff_irreducible(u, v, w):
                        while is_coeff_irreducible(u, v, w) and w0 != w:
                            w0 = w
                            u, v, w = reduce_descents(u, v, w)
                            if is_coeff_irreducible(u, v, w):
                                u, v, w = reduce_coeff(u, v, w)

    if w != w2 and sign_only:
        # logger.debug(f"Return 0 ")
        return 0
    # logger.debug(f"Reduced to {u2=} {v2=} {w2=} {val=}")
    if is_coeff_irreducible(u, v, w):
        u3, v3, w3 = try_reduce_v(u, v, w)
        if not is_coeff_irreducible(u3, v3, w3):
            u, v, w = u3, v3, w3
        else:
            u3, v3, w3 = try_reduce_u(u, v, w)
            if not is_coeff_irreducible(u3, v3, w3):
                u, v, w = u3, v3, w3
    split_two_b, split_two = is_split_two(u, v, w)
    # logger.debug("Recording line number")
    if len([i for i in code(v) if i != 0]) == 1:
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        cv = code(v)
        for i in range(len(cv)):
            if cv[i] != 0:
                k = i + 1
                p = cv[i]
                break
        inv_u = inv(u)
        r = inv(w) - inv_u
        val = 0
        w2 = w
        hvarset = [w2[i] for i in range(min(len(w2), k))] + [i + 1 for i in range(len(w2), k)] + [w2[b] for b in range(k, len(u)) if u[b] != w2[b]] + [w2[b] for b in range(len(u), len(w2))]
        # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        return elem_sym_poly(
            p - r,
            k + p - 1,
            [-var3[i] for i in range(1, n)],
            [-var2[i] for i in hvarset],
        )
        # if expand(val - oldval) != 0:
        #     # logger.debug("This is bad")
        #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
    if will_formula_work(v, u) or dominates(u, w):
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        return dualcoeff(u, v, w, var2, var3)
        # if expand(val - oldval) != 0:
        # logger.debug("This is bad")
        # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=} {will_formula_work(v,u)=} {dominates(u,w)=}")
        # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
    if not v.has_pattern([1, 4, 2, 3]) and not v.has_pattern([4, 1, 3, 2]) and not v.has_pattern([3, 1, 4, 2]) and not v.has_pattern([1, 4, 3, 2]):
        logger.debug("Recording new characterization was used")
        return schubmult_double({u: 1}, v, var2, var3).get(w, 0)

    if inv(w) - inv(u) == 1:
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        a, b = -1, -1
        for i in range(len(w)):
            if a == -1 and u[i] != w[i]:
                a = i
            elif (i >= len(u) and w[i] != i + 1) or (b == -1 and u[i] != w[i]):
                b = i
        arr = [[[], v]]
        d = -1
        for i in range(len(v) - 1):
            if v[i] > v[i + 1]:
                d = i + 1
        for i in range(d):
            arr2 = []
            if i in [a, b]:
                continue
            i2 = 1
            if i > b:
                i2 += 2
            elif i > a:
                i2 += 1
            for vr, v2 in arr:
                dpret = pull_out_var(i2, v2)
                for v3r, v3 in dpret:
                    arr2 += [[[*vr, v3r], v3]]
            arr = arr2
        val = 0
        for L in arr:
            v3 = L[-1]
            if v3[0] < v3[1]:
                continue
            v3 = v3.swap(0, 1)
            toadd = 1
            for i in range(d):
                if i in [a, b]:
                    continue
                i2 = i
                if i > b:
                    i2 = i - 2
                elif i > a:
                    i2 = i - 1
                oaf = L[0][i2]
                if i >= len(w):
                    yv = i + 1
                else:
                    yv = w[i]
                for j in range(len(oaf)):
                    toadd *= var2[yv] - var3[oaf[j]]
            toadd *= schubpoly(v3, [0, var2[w[a]], var2[w[b]]], var3)
            val += toadd
        # if expand(val - oldval) != 0:
        # logger.debug("This is bad")
        # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
        # logger.debug(f"good to go {u2=} {v2=} {w2=}")
        return val
    # if split_two_b:
    #     # logger.debug("Recording line number")
    #     if sign_only:
    #         return 0
    #     cycles = split_two
    #     a1, b1 = cycles[0]
    #     a2, b2 = cycles[1]
    #     a1 -= 1
    #     b1 -= 1
    #     a2 -= 1
    #     b2 -= 1
    #     spo = sorted([a1, b1, a2, b2])
    #     real_a1 = min(spo.index(a1), spo.index(b1))
    #     real_a2 = min(spo.index(a2), spo.index(b2))
    #     real_b1 = max(spo.index(a1), spo.index(b1))
    #     real_b2 = max(spo.index(a2), spo.index(b2))

    #     good1 = False
    #     good2 = False
    #     if real_b1 - real_a1 == 1:
    #         good1 = True
    #     if real_b2 - real_a2 == 1:
    #         good2 = True
    #     a, b = -1, -1
    #     if good1 and not good2:
    #         a, b = min(a2, b2), max(a2, b2)
    #     if good2 and not good1:
    #         a, b = min(a1, b1), max(a1, b1)
    #     arr = [[[], v]]
    #     d = -1
    #     for i in range(len(v) - 1):
    #         if v[i] > v[i + 1]:
    #             d = i + 1
    #     for i in range(d):
    #         arr2 = []

    #         if i in [a1, b1, a2, b2]:
    #             continue
    #         i2 = 1
    #         i2 += len([aa for aa in [a1, b1, a2, b2] if i > aa])
    #         for vr, v2 in arr:
    #             dpret = pull_out_var(i2, v2)
    #             for v3r, v3 in dpret:
    #                 arr2 += [[[*vr, (v3r, i + 1)], v3]]
    #         arr = arr2
    #     val = 0

    #     if good1:
    #         arr2 = []
    #         for L in arr:
    #             v3 = L[-1]
    #             if v3[real_a1] < v3[real_b1]:
    #                 continue
    #             v3 = v3.swap(real_a1, real_b1)
    #             arr2 += [[L[0], v3]]
    #         arr = arr2
    #         if not good2:
    #             for i in range(4):
    #                 arr2 = []

    #                 if i in [real_a2, real_b2]:
    #                     continue
    #                 if i == real_a1:
    #                     var_index = min(a1, b1) + 1
    #                 elif i == real_b1:
    #                     var_index = max(a1, b1) + 1
    #                 i2 = 1
    #                 i2 += len([aa for aa in [real_a2, real_b2] if i > aa])
    #                 for vr, v2 in arr:
    #                     dpret = pull_out_var(i2, v2)
    #                     for v3r, v3 in dpret:
    #                         arr2 += [[[*vr, (v3r, var_index)], v3]]
    #                 arr = arr2
    #     if good2:
    #         arr2 = []
    #         for L in arr:
    #             v3 = L[-1]
    #             try:
    #                 if v3[real_a2] < v3[real_b2]:
    #                     continue
    #                 v3 = v3.swap(real_a2, real_b2)
    #             except IndexError:
    #                 continue
    #             arr2 += [[L[0], v3]]
    #         arr = arr2
    #         if not good1:
    #             for i in range(4):
    #                 arr2 = []

    #                 if i in [real_a1, real_b1]:
    #                     continue
    #                 i2 = 1
    #                 i2 += len([aa for aa in [real_a1, real_b1] if i > aa])
    #                 if i == real_a2:
    #                     var_index = min(a2, b2) + 1
    #                 elif i == real_b2:
    #                     var_index = max(a2, b2) + 1
    #                 for vr, v2 in arr:
    #                     dpret = pull_out_var(i2, v2)
    #                     for v3r, v3 in dpret:
    #                         arr2 += [[[*vr, (v3r, var_index)], v3]]
    #                 arr = arr2

    #         for L in arr:
    #             v3 = L[-1]
    #             tomul = 1
    #             doschubpoly = True
    #             if (not good1 or not good2) and v3[0] < v3[1] and (good1 or good2):
    #                 continue
    #             if (good1 or good2) and (not good1 or not good2):
    #                 v3 = v3.swap(0, 1)
    #             elif not good1 and not good2:
    #                 doschubpoly = False
    #                 if v3[0] < v3[1]:
    #                     dual_u = uncode([2, 0])
    #                     dual_w = Permutation([4, 2, 1, 3])
    #                     coeff = perm_act(dualcoeff(dual_u, v3, dual_w, var2, var3), 2, var2)

    #                 elif len(v3) < 3 or v3[1] < v3[2]:
    #                     if len(v3) <= 3 or v3[2] < v3[3]:
    #                         coeff = 0
    #                         continue
    #                     v3 = v3.swap(0, 1).swap(2, 3)
    #                     coeff = perm_act(schubpoly(v3, var2, var3), 2, var2)
    #                 elif len(v3) <= 3 or v3[2] < v3[3]:
    #                     if len(v3) <= 3:
    #                         v3 += [4]
    #                     v3 = v3.swap(2, 3)
    #                     coeff = perm_act(
    #                         posify(
    #                             schubmult_one(Permutation([1, 3, 2]), v3, var2, var3).get(
    #                                 Permutation([2, 4, 3, 1]),
    #                                 0,
    #                             ),
    #                             Permutation([1, 3, 2]),
    #                             v3,
    #                             Permutation([2, 4, 3, 1]),
    #                             var2,
    #                             var3,
    #                             msg,
    #                             do_pos_neg,
    #                             optimize=optimize,
    #                         ),
    #                         2,
    #                         var2,
    #                     )
    #                     # logger.debug(f"{coeff=}")
    #                 else:
    #                     coeff = perm_act(
    #                         schubmult_one(Permutation([1, 3, 2]), v3, var2, var3).get(
    #                             Permutation([2, 4, 1, 3]),
    #                             0,
    #                         ),
    #                         2,
    #                         var2,
    #                     )
    #                 # logger.debug(f"{coeff=}")
    #                 # if expand(coeff) == 0:
    #                 #     # logger.debug("coeff 0 oh no")
    #                 tomul = sympify(coeff)
    #             toadd = 1
    #             for i in range(len(L[0])):
    #                 var_index = L[0][i][1]
    #                 oaf = L[0][i][0]
    #                 if var_index - 1 >= len(w):
    #                     yv = var_index
    #                 else:
    #                     yv = w[var_index - 1]
    #                 for j in range(len(oaf)):
    #                     toadd *= var2[yv] - var3[oaf[j]]
    #             if (not good1 or not good2) and (good1 or good2):
    #                 varo = [0, var2[w[a]], var2[w[b]]]
    #             else:
    #                 varo = [0, *[var2[w[spo[k]]] for k in range(4)]]
    #             if doschubpoly:
    #                 toadd *= schubpoly(v3, varo, var3)
    #             else:
    #                 subs_dict3 = {var2[i]: varo[i] for i in range(len(varo))}
    #                 toadd *= efficient_subs(tomul, subs_dict3)
    #             val += toadd
    #             # logger.debug(f"accum {val=}")
    #         #logger.debug(f"{expand(val-oldval)=}")
    #         # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
    #         return val
    if will_formula_work(u, v):
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        return forwardcoeff(u, v, w, var2, var3)
        # if expand(val - oldval) != 0:
        #     # logger.debug("This is bad")
        #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
    # logger.debug("Recording line number")
    # c01 = code(u)
    # c02 = code(w)
    # c03 = code(v)

    c1 = code(~u)
    c2 = code(~w)

    if one_dominates(u, w):
        if sign_only:
            return 0
        while c1[0] != c2[0]:
            w = w.swap(c2[0] - 1, c2[0])
            v = v.swap(c2[0] - 1, c2[0])
            # w[c2[0] - 1], w[c2[0]] = w[c2[0]], w[c2[0] - 1]
            # v[c2[0] - 1], v[c2[0]] = v[c2[0]], v[c2[0] - 1]
            # w = tuple(w)
            # v = tuple(v)
            c2 = code(~w)
            # c03 = code(v)
            # c01 = code(u)
            # c02 = code(w)
        # if is_reducible(v):
        #     # logger.debug("Recording line number")
        #     if sign_only:
        #         return 0
        #     newc = []
        #     elemc = []
        #     for i in range(len(c03)):
        #         if c03[i] > 0:
        #             newc += [c03[i] - 1]
        #             elemc += [1]
        #         else:
        #             break
        #     v3 = uncode(newc)
        #     coeff_dict = schubmult_one(
        #         u,
        #         uncode(elemc),
        #         var2,
        #         var3,
        #     )
        #     val = 0
        #     for new_w in coeff_dict:
        #         tomul = coeff_dict[new_w]
        #         newval = schubmult_one(new_w, uncode(newc), var2, var3).get(
        #             w,
        #             0,
        #         )
        #         # logger.debug(f"Calling posify on {newval=} {new_w=} {uncode(newc)=} {w=}")
        #         newval = posify(newval, new_w, uncode(newc), w, var2, var3, msg, do_pos_neg, optimize=optimize,elem_sym=elem_sym)
        #         val += tomul * shiftsubz(newval)
        #     # if expand(val - oldval) != 0:
        #     #     # logger.debug("This is bad")
        #     #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
        #         # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        #     return val
        # removed, iffy (hard to implement)
        # if c01[0] == c02[0] and c01[0] != 0:
        #     # logger.debug("Recording line number")
        #     if sign_only:
        #         return 0
        #     varl = c01[0]
        #     u3 = uncode([0] + c01[1:])
        #     w3 = uncode([0] + c02[1:])
        #     val = 0
        #     val = schubmult_one(u3, v, var2, var3).get(
        #         w3,
        #         0,
        #     )
        #     # logger.debug(f"Calling posify on {val=} {u3=} {v=} {w3=}")
        #     val = posify(val, u3, v, w3, var2, var3, msg, do_pos_neg, optimize=optimize,elem_sym=elem_sym)
        #     for i in range(varl):
        #         val = perm_act(val, i + 1, var2)
        #     # if expand(val - oldval) != 0:
        #     #     # logger.debug("This is bad")
        #     #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
        #     # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        #     return val
        if c1[0] == c2[0]:
            # logger.debug("Recording line number")
            if sign_only:
                return 0
            vp = pull_out_var(c1[0] + 1, v)
            u3 = phi1(u)
            w3 = phi1(w)
            val = 0
            for arr, v3 in vp:
                tomul = 1
                for i in range(len(arr)):
                    tomul *= var2[1] - var3[arr[i]]

                val2 = schubmult_double_pair(u3, v3, var2, var3).get(
                    w3,
                    0,
                )
                val2 = posify(val2, u3, v3, w3, var2, var3, msg, do_pos_neg, optimize=optimize,elem_sym=elem_sym)
                val += tomul * shiftsub(val2, var2)
            # if expand(val - oldval) != 0:
            #     # logger.debug("This is bad")
            #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=")
            # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
            return val
    # logger.debug("Fell all the way through. Cleverness did not save us")
    if not sign_only:
        # logger.debug("Recording line number")
        if optimize:
            if elem_sym:
                #print(f"{elem_sym=}")
                val2 = compute_positive_rep_new(elem_sym, var2, var3, msg, False)
            elif inv(u) + inv(v) - inv(w) == 1:
                val2 = compute_positive_rep(val, var2, var3, msg, False)
            else:
                val2 = compute_positive_rep(val, var2, var3, msg, do_pos_neg)
            if val2 is not None:
                val = val2
            return val
        if optimize is None:
            raise Exception
        # logger.debug("RETURNINGOLDVAL")
        return oldval
    # logger.debug("Recording line number")
    d = expand(val).as_coefficients_dict()
    for v in d.values():
        if v < 0:
            return -1
    return 1


def shiftsub(pol, var2=None):
    subs_dict = {var2[i]: var2[i + 1] for i in range(99)}
    return efficient_subs(sympify(pol), subs_dict)


def shiftsubz(pol, var3=None):
    subs_dict = {var3[i]: var3[i + 1] for i in range(99)}
    return efficient_subs(sympify(pol), subs_dict)


def split_flat_term(arg, genset):
    arg = expand(arg)
    ys = []
    zs = []
    for arg2 in arg.args:
        if any(genset.index(x) != -1 for x in arg2.free_symbols):
            if isinstance(arg2, Mul):
                for i in range(int(arg2.args[0])):
                    ys += [arg2.args[1]]
            else:
                ys += [arg2]
        elif isinstance(arg2, Mul):
            for i in range(abs(int(arg2.args[0]))):
                zs += [-arg2.args[1]]
        else:
            zs += [arg2]
    return ys, zs


def is_flat_term(term):
    if isinstance(term, Integer) or isinstance(term, int):
        return True
    dc = expand(term).as_coefficients_dict()
    return all(isinstance(x, Symbol) for x in dc.keys())


def flatten_factors(term):
    found_one = False
    if is_flat_term(term):
        return term, False
    if isinstance(term, Pow):
        if is_flat_term(term.args[0]) and len(term.args[0].args) > 2:
            ys, zs = split_flat_term(term.args[0])
            terms = [1]
            for i in range(len(ys)):
                terms2 = []
                for j in range(len(term.args[1])):
                    for t in terms:
                        terms2 += [t * (ys[i] + zs[i])]
                terms = terms2
            return Add(*terms)
        if is_flat_term(term.args[0]):
            return term, False
        return flatten_factors(term.args[0]) ** term.args[1], True
    if isinstance(term, Mul):
        terms = [1]
        for arg in term.args:
            terms2 = []
            if isinstance(arg, Add) and not is_flat_term(expand(arg)):
                found_one = True
                for term3 in terms:
                    for arg2 in arg.args:
                        flat, found = flatten_factors(arg2)
                        terms2 += [term3 * flat]
            elif isinstance(arg, Add) and is_flat_term(arg) and len(arg.args) > 2:
                found_one = True
                ys, zs = split_flat_term(arg)
                for term3 in terms:
                    for i in range(len(ys)):
                        terms2 += [term3 * (ys[i] + zs[i])]
            else:
                flat, found = flatten_factors(arg)
                if found:
                    found_one = True
                for term3 in terms:
                    terms2 += [term3 * flat]
            terms = terms2
        if len(terms) == 1:
            term = terms[0]
        else:
            term = Add(*terms)
        return term, found_one
    if isinstance(term, Add):
        res = 0
        for arg in term.args:
            flat, found = flatten_factors(arg)
            if found:
                found_one = True
            res += flat
        return res, found_one
    return None


def fres(v):
    for s in v.free_symbols:
        return s
    return None


def split_mul(arg0, var2=None, var3=None):
    monoms = SortedList()

    var2s = {fres(var2[i]): i for i in range(len(var2))}
    var3s = {fres(var3[i]): i for i in range(len(var3))}
    # print(f"{type(arg0)=} {arg0=}")
    if isinstance(arg0, Pow):
        arg = arg0
        arg2 = expand(arg.args[0])
        yval = arg2.args[0]
        zval = arg2.args[1]
        if str(yval).find("z") != -1:
            yval, zval = zval, yval
        if str(zval).find("-") != -1:
            zval = -zval
        if str(yval).find("-") != -1:
            yval = -yval
        tup = (var2s[fres(yval)], var3s[fres(zval)])
        for i in range(int(arg0.args[1])):
            monoms += [tup]
    else:
        for arg in arg0.args:
            if is_flat_term(arg):
                if isinstance(arg, Integer) or isinstance(arg, int):
                    continue
                arg = expand(arg)
                if arg == 0:
                    break
                yval = arg.args[0]
                zval = arg.args[1]
                if str(yval).find("z") != -1:
                    yval, zval = zval, yval
                if str(zval).find("-") != -1:
                    zval = -zval
                if str(yval).find("-") != -1:
                    yval = -yval
                monoms += [(var2s[fres(yval)], var3s[fres(zval)])]
            elif isinstance(arg, Pow):
                arg2 = arg.args[0]
                yval = arg2.args[0]
                zval = arg2.args[1]
                if str(yval).find("z") != -1:
                    yval, zval = zval, yval
                if str(zval).find("-") != -1:
                    zval = -zval
                if str(yval).find("-") != -1:
                    yval = -yval
                tup = (var2s[fres(yval)], var3s[fres(zval)])
                for i in range(int(arg.args[1])):
                    monoms += [tup]
    return monoms


def split_monoms(pos_part, var2, var3):
    arrs = SortedList()
    if isinstance(pos_part, Add):
        for arg0 in pos_part.args:
            monoms = split_mul(arg0, var2, var3)
            arrs += [monoms]
    elif isinstance(pos_part, Mul) or isinstance(pos_part, Pow):
        monoms = split_mul(pos_part, var2, var3)
        arrs += [monoms]
    else:
        return [pos_part]
    return arrs


def is_negative(term):
    sign = 1
    if isinstance(term, Integer) or isinstance(term, int):
        return term < 0
    if isinstance(term, Mul):
        for arg in term.args:
            if isinstance(arg, Integer):
                sign *= arg
            elif isinstance(arg, Add):
                if str(arg).find("-y") != -1:
                    sign *= -1
            elif isinstance(arg, Pow):
                mulsign = 1
                if str(arg.args[0]).find("-y") != -1:
                    mulsign = -1
                sign *= mulsign**term.index
    elif isinstance(term, Pow):
        mulsign = 1
        if str(term.args[0]).find("-y") != -1:
            mulsign = -1
        sign *= mulsign**term.index
    return sign < 0


def find_base_vectors(monom_list, var2, var3, depth):
    size = 0
    mn_fullcount = {}
    # pairs_checked = set()
    monom_list = {tuple(mn) for mn in monom_list}
    ct = 0
    while ct < depth and size != len(monom_list):
        size = len(monom_list)
        monom_list2 = set(monom_list)
        additional_set2 = set()
        for mn in monom_list:
            mncount = mn_fullcount.get(mn, {})
            if mncount == {}:
                for tp in mn:
                    mncount[tp] = mncount.get(tp, 0) + 1
                mn_fullcount[mn] = mncount
            for mn2 in monom_list:
                mn2count = mn_fullcount.get(mn2, {})
                if mn2count == {}:
                    for tp in mn2:
                        mn2count[tp] = mn2count.get(tp, 0) + 1
                    mn_fullcount[mn2] = mn2count
                num_diff = 0
                for tp in mncount:
                    pt = mn2count.get(tp, 0) - mncount[tp]
                    num_diff += abs(pt)
                    if num_diff > 1:
                        break
                if num_diff == 1:
                    diff_term1 = None
                    diff_term2 = None
                    for tp in mn2count:
                        if mn2count[tp] > mncount.get(tp, 0):
                            diff_term2 = tp
                            break
                    for tp2 in mncount:
                        if mncount[tp2] > mn2count.get(tp2, 0):
                            diff_term1 = tp2
                            break
                    # print(f"{mn,mn2}")
                    if diff_term1 is None or diff_term2 is None:
                        raise Exception(f"{mn=} {mn2=}")
                    if diff_term2[1] == diff_term1[1]:
                        continue
                    new_term1 = (diff_term1[0], diff_term2[1])
                    new_term2 = (diff_term2[0], diff_term1[1])
                    # mn3 = [*mn]
                    # mn4 = list(mn2)
                    index = bisect_left(mn, diff_term1)
                    mn3 = list(mn[:index]) + list(mn[index + 1 :])
                    index = bisect_left(mn3, new_term1)
                    mn3_t = tuple(mn3[:index] + [new_term1] + mn3[index:])
                    index2 = bisect_left(mn2, diff_term2)
                    mn4 = list(mn2[:index2]) + list(mn2[index2 + 1 :])
                    index2 = bisect_left(mn4, new_term2)
                    mn4_t = tuple(mn4[:index2] + [new_term2] + mn4[index2:])
                    if mn3_t not in monom_list2:
                        additional_set2.add(mn3_t)
                    monom_list2.add(mn3_t)
                    if mn4_t not in monom_list2:
                        additional_set2.add(mn4_t)
                    monom_list2.add(mn4_t)
        monom_list = monom_list2
        ct += 1
    ret = []
    for mn in monom_list:
        if len(mn) != len(set(mn)):
            continue
        res = 1
        for tp in mn:
            res *= var2[tp[0]] - var3[tp[1]]
        ret += [res]
    return ret, monom_list


def posify_generic_partial(val, u2, v2, w2):
    val2 = val
    val = posify(val, u2, v2, w2, var2=None, var3=None, msg=True, do_pos_neg=False, sign_only=False, optimize=False)
    if expand(val - val2) != 0:
        # logger.debug("Warning, failed on a case")
        raise Exception(f"{val=} {val2=} {u2=} {v2=} {w2=}")
    # print("FROFL")
    return val


@cache
def schubmult_generic_partial_posify(u2, v2):
    return {w2: posify_generic_partial(val, u2, v2, w2) for w2, val in schubmult_double_pair_generic(u2, v2).items()}


def forwardcoeff(u, v, perm, var2=None, var3=None):
    th = theta(v)
    muv = uncode(th)
    vmun1 = (~v) * muv

    w = perm * vmun1
    if inv(w) == inv(vmun1) + inv(perm):
        coeff_dict = schubmult_double_pair(u, muv, var2, var3)
        # logger.debug(f"{coeff_dict.get(w,0)=} {w=} {perm=} {vmun1=} {v=} {muv=}")
        return coeff_dict.get(w, 0)
    return 0


def dualcoeff(u, v, perm, var2=None, var3=None):
    if inv(u) == 0:
        # logger.debug("Recording line number")
        vp = v * (~perm)
        if inv(vp) == inv(v) - inv(perm):
            return schubpoly(vp, var2, var3)
    dpret = []
    ret = 0
    if dominates(u, perm):
        dpret = dualpieri(u, v, perm)
    else:
        # logger.debug("Recording line number")
        dpret = []
        # logger.debug("Recording line number")
        th = theta(u)
        muu = uncode(th)
        umun1 = (~u) * muu
        w = perm * umun1
        # logger.debug("spiggle")
        # logger.debug(f"{u=} {muu=} {v=} {w=} {perm=}")
        # logger.debug(f"{w=} {perm=}")
        if inv(w) == inv(umun1) + inv(perm):
            dpret = dualpieri(muu, v, w)
            # logger.debug(f"{muu=} {v=} {w=}")
            # logger.debug(f"{dpret=}")
    for vlist, vp in dpret:
        # logger.debug("Recording line number")
        toadd = 1
        for i in range(len(vlist)):
            for j in range(len(vlist[i])):
                toadd *= var2[i + 1] - var3[vlist[i][j]]
        toadd *= schubpoly(vp, var2, var3, len(vlist) + 1)
        ret += toadd
    return ret
    # logger.debug("Recording line number")
    # schub_val = schubmult_one(u, v, var2, var3)
    # val_ret = schub_val.get(perm, 0)
    # if expand(val - val_ret) != 0:
    #     # logger.debug(f"{schub_val=}")
    #     # logger.debug(f"{val=} {u=} {v=} {var2[1]=} {var3[1]=}  {perm=} {schub_val.get(perm,0)=}")
    # logger.debug(f"good to go {ret=}")


def dualpieri(mu, v, w):
    # logger.debug(f"dualpieri {mu=} {v=} {w=}")
    lm = code(~mu)
    cn1w = code(~w)
    while len(lm) > 0 and lm[-1] == 0:
        lm.pop()
    while len(cn1w) > 0 and cn1w[-1] == 0:
        cn1w.pop()
    if len(cn1w) < len(lm):
        return []
    for i in range(len(lm)):
        if lm[i] > cn1w[i]:
            return []
    c = Permutation([1, 2])
    for i in range(len(lm), len(cn1w)):
        c = cycle(i - len(lm) + 1, cn1w[i]) * c
    # c = permtrim(c)
    # logger.debug("Recording line number")
    res = [[[], v]]
    # logger.debug(f"{v=} {type(v)=}")
    for i in range(len(lm)):
        # logger.debug(f"{res=}")
        res2 = []
        for vlist, vplist in res:
            vp = vplist
            vpl = divdiffable(vp, cycle(lm[i] + 1, cn1w[i] - lm[i]))
            # logger.debug(f"{vpl=} {type(vpl)=}")
            if len(vpl) == 0:
                continue
            vl = pull_out_var(lm[i] + 1, vpl)
            # logger.debug(f"{vl=}")
            for pw, vpl2 in vl:
                res2 += [[[*vlist, pw], vpl2]]
        res = res2
    if len(lm) == len(cn1w):
        return res
    res2 = []
    for vlist, vplist in res:
        vp = vplist
        vpl = divdiffable(vp, c)
        if len(vpl) == 0:
            continue
        res2 += [[vlist, vpl]]
    # logger.debug(f"{res2=}")
    return res2


class Optimizer:
    def __init__(self, z_ring, poly):
        self.z_ring=z_ring
        self.monom_to_vec = {}
        self.vec0 = None
        self.vec0 = self.poly_to_vec(poly)

    def _init_basevec(self, dc):
        self.monom_to_vec = {}
        index = 0
        for mn in dc:
            if dc[mn] == 0:
                continue
            self.monom_to_vec[mn] = index
            index += 1

    def poly_to_vec(self, poly):
        poly = expand(sympify(poly).xreplace({self.z_ring.genset[1]: 0}))

        dc = poly.as_coefficients_dict()

        if self.vec0 is None:
            self._init_basevec(dc)

        vec = {}
        for mn in dc:
            cf = dc[mn]
            if cf == 0:
                continue
            cf = abs(int(cf))
            try:
                index = self.monom_to_vec[mn]
            except KeyError:
                return None
            if self.vec0 is not None and self.vec0[index] < cf:
                return None
            vec[index] = cf
        return vec


class SchubOptimizer(Optimizer):
    # spake = z_ring(sympify(expand_func(v)).xreplace(vnv))
    # spake = zs_ring.from_dict({k:expand(spa) for k,spa in spake.items()})
    # #print(f"{v}={spake}")
    # fartspangle = {}
    # for k,pork in spake.items():
    #     spink = pork.as_coefficients_dict()
    #     for sperg, po in spink.items():
    #         if po == 0:
    #             continue
    #         fartspangle[sperg] = fartspangle.get(sperg,0) + 1
    def __init__(self, z_ring, poly):
        self.z_ring = z_ring
        self.pos_dict = {z_ring.genset[i]: -z_ring.genset[i] for i in range(100)}
        self.vec0 = self.poly_to_vec(poly)

    # def _init_basevec(self, dc):
    #     self.monom_to_vec = {}
    #     index = 0
    #     for mn in dc:
    #         if dc[mn] == 0:
    #             continue
    #         self.monom_to_vec[mn] = index
    #         index += 1

    def poly_to_vec(self, poly):
        spake = self.z_ring(poly.xreplace(self.pos_dict))
        spake = self.z_ring.from_dict({k: expand(spa) for k, spa in spake.items()})
        # #print(f"{v}={spake}")
        monom_dct = {}
        for k, pork in spake.items():
            spink = pork.as_coefficients_dict()
            for sperg, po in spink.items():
                if po == 0:
                    continue
                monom_dct[sperg] = monom_dct.get(sperg, {})
                monom_dct[sperg][k] = monom_dct[sperg].get(k, 0) + int(po)
        return monom_dct
        # poly = expand(sympify(poly).xreplace({self.z_ring.genset[1]: 0}))

        # dc = poly.as_coefficients_dict()

        # if self.vec0 is None:
        #     self._init_basevec(dc)

        # vec = {}
        # for mn in dc:
        #     cf = dc[mn]
        #     if cf == 0:
        #         continue
        #     cf = abs(int(cf))
        #     try:
        #         index = self.monom_to_vec[mn]
        #     except KeyError:
        #         return None
        #     if self.vec0 is not None and self.vec0[index] < cf:
        #         return None
        #     vec[index] = cf
        # return vec
