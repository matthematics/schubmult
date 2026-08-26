from schubmult import *
from schubmult.symbolic.common_polys import grothendieck_poly
from schubmult.symbolic import S, sympify, prod, expand_func, expand, Integer, Mul, Add, Pow, Symbol
from schubmult.abc import E, y, z, H
from schubmult.symbolic.common_polys import efficient_subs
from functools import cache
from schubmult.symbolic.common_polys import _vars, efficient_subs, elem_func_func_mul, elem_sym_func, elem_sym_poly
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict, add_perm_dict_with_coeff
from schubmult.utils.schub_lib import (
    compute_vpathdicts,
    elem_sym_perms,
    elem_sym_perms_op,
    elem_sym_positional_perms,
    pull_out_var,
)
from schubmult.mult.separated_descents import grothmult_double_plus
#from schubmult.utils._mul_utils import add_perm_dict_with_coeff
beta = Gx._beta

#negy_dict = {y[i]: -y[i] for i in range(20)}

@cache
def _elem_sym_perms(perm, k):
    return elem_sym_perms(perm, k, k)

def groth_w0(n, ring, beta):
    ring = DSx([]).ring
    genset = ring.genset
    coeff_genset = ring.coeff_genset
    start_elem = ring.one
    for i in range(n - 1, 0, -1):
        new_start_elem = 0
        yvar = coeff_genset[n - i]
        bvar = (1 + beta * yvar)
        var2 = [coeff_genset[j] * bvar for j in range(1, 25)]
        
        for permo, coeff in start_elem.items():
            #new_elem = ringa(permo) * E(i, i, ring.genset[1:], ringb.coeff_genset[1:])
            for elem_perm, diff in _elem_sym_perms(permo, i):
                new_start_elem += coeff * (bvar**diff) * prod([var2[permo[p] - 1] + yvar for p in range(i) if permo[p] == elem_perm[p]]) * ring(elem_perm)
            #new_start_elem += coeff * ring.from_dict({k: (1 + beta*ring.coeff_genset[n - i])**((k.inv - permo.inv)) * v for k, v in new_elem.items()})
        # a_j to y_j(1 + y_{n - i})
        start_elem = new_start_elem
        #start_elem = new_start_elem
    return start_elem


def dom_groth(dom_perm, ring, beta):
    ring = DSx([]).ring
    genset = ring.genset
    coeff_genset = ring.coeff_genset
    start_elem = ring.one    
    lengths = (~dom_perm).trimcode
    n = len(lengths) + 1
    for i in range(n - 1, 0, -1):
        new_start_elem = 0
        yvar = coeff_genset[n - i]
        bvar = (1 + beta * yvar)
        var2 = [coeff_genset[j] * bvar for j in range(1, 25)]
        length = lengths[n - i - 1]
        for permo, coeff in start_elem.items():
            for elem_perm, diff in _elem_sym_perms(permo, length):
                new_start_elem += coeff * (bvar**diff) * prod([var2[permo[p] - 1] + yvar for p in range(length) if permo[p] == elem_perm[p]]) * ring(elem_perm)
        start_elem = new_start_elem
    return start_elem

def groth_poly(perm, beta, return_dict=False):
    
    ring = DSx([]).ring
    ringt = DSx([], "t").ring
    n = len(perm)
    swp = (~perm) * Permutation.w0(n)
    start_dict = {Permutation([]): S.One}
    n = len(perm)
    for i in range(1, n):
        new_start_dict = {}
        for permo, coeff in start_dict.items():
            #new_start_dict = add_perm_dict_with_coeff(grothmult_double_plus({uncode([1] * (n - i)): coeff}, permo, ringt.genset, ring.coeff_genset[i - 1:], beta=beta), new_start_dict, coeff=1)
            new_start_dict = add_perm_dict_with_coeff(grothmult_double_plus({uncode([1] * (n - i)): coeff}, permo, ringt.genset, ring.coeff_genset[i - 1:], beta=beta), new_start_dict, coeff=1)
        start_dict = new_start_dict
    start_dict = {k: v.subs({ringt.genset[i]: ring.genset[(~k)[i - 1]] for i in range(1, 20)}).expand() for k, v in start_dict.items()}
    if return_dict:
        return {k: v.expand() for k, v in start_dict.items()}
    return start_dict[swp]

def groth_poly_pull(perm, beta, varnum):    
    ring = DSx([]).ring
    # start_elem = ring.one
    start_dict = {Permutation([]): S.One}
    n = len(perm)
    # if perm == Permutation.w0(n):
    #     return groth_w0(n, beta)
    met = False
    t = DSx([], "t").ring.coeff_genset
    dom_elem1 = ~uncode(list(range(n - 1, n - 1 - varnum, -1)))
    dct = grothmult_double_plus({uncode([1] * (n - varnum))}, dom_elem1, ring.coeff_genset,  ring.genset[varnum-1:], beta=beta)
    low_dom_elem = ~uncode(list(range(n - 1 - varnum, 0, -1)))
    new_dct = {k * low_dom_elem: v for k, v in dct.items() if (k * low_dom_elem).inv == k.inv + low_dom_elem.inv}
    new_new_dict = {~(k * (~perm)): v for k, v in new_dct.items() if (k * (~perm)).inv == k.inv - (~perm).inv}
    return new_new_dict


def isobaric_strip_on_dschub(start, length, schub_perm, beta):
    start_schub = DSx(schub_perm)
    
    ring = start_schub.ring
    #var_strip = ring.genset[start + 1:start + length + 1]
    #one_genset = [-(beta**(-1)) for _ in range(50)]
    bigger_schub = ring.zero
    positions = list(range(start + 1, start + length + 1))
    for perm, coeff in start_schub.items():
        for elem_perm, diff, sign in elem_sym_positional_perms(perm, length, *positions):
            bigger_schub += sign * coeff * (beta**diff) * prod([ring.coeff_genset[perm[positions[p] - 1]]*beta + 1 for p in range(length) if perm[positions[p] - 1] == elem_perm[positions[p] - 1]]) *ring(elem_perm)
    stripness = list(range(start, start + length))
    ret_schub = ring.from_dict({k: v for k, v in bigger_schub.items()})
    for desc in stripness:
        ret_schub = ring.from_dict({perm2.swap(desc - 1, desc): v for perm2, v in ret_schub.items() if perm2[desc - 1] > perm2[desc]})
    return ret_schub

@cache
def apply_isobaric_to_schub(diff_perm, schub_perm, beta):
    elem = DSx(schub_perm)
    
    strips = [[i, (diff_perm).trimcode[i - 1]] for i in range(1, (diff_perm).max_descent + 1)]
    for strip in reversed(strips):
        if strip[1] == 0:
            continue
        new_elem = elem.ring.zero
        for perm, coeff in elem.items():
            new_elem += coeff * isobaric_strip_on_dschub(strip[0], strip[1], perm, beta=beta)
        elem = new_elem
    return elem

def alt_grothendieck_poly(perm, beta):
    from schubmult.abc import z
    #dom_perm = Permutation.w0(len(perm))
    dom_perm = perm.minimal_dominant_above()
    #Permutation.w0(len(perm))
    diff_perm = (~perm) * dom_perm
    ring = DSx([]).ring
    #first_potato = groth_w0(len(perm), ring, beta=beta)
    first_potato = dom_groth(dom_perm, ring, beta=beta)
    schub_elem = ring.zero
    for perm2, coeff in first_potato.items():
        schub_elem += coeff * apply_isobaric_to_schub(diff_perm, perm2, beta=beta)
    return ring.from_dict({k: v.expand() for k, v in schub_elem.items()})

def pull_out_groth_var2(perm, beta, varnum):
    from schubmult.abc import z
    #dom_perm = Permutation.w0(len(perm))
    ring = DSx([]).ring
    if varnum > perm.max_descent:
        return {perm: S.One}
    work_perm = (~perm)
    
    dom_perm = work_perm.minimal_dominant_above()

    cd = (~dom_perm).trimcode[:varnum - 1]
    if varnum  < len((~dom_perm).trimcode):
        dom_chop_cd = (~dom_perm).trimcode[varnum:]
    else:
        dom_chop_cd = []
    new_dom_perm = ~uncode(cd)
    chopdom = ~uncode(dom_chop_cd)
    start_dict = grothmult_double_plus({uncode([1]* (~dom_perm).trimcode[varnum - 1]): 1}, new_dom_perm, ring.coeff_genset, ring.genset[varnum - 1:], beta=beta)
    #Permutation.w0(len(perm))
    diff_perm = (~work_perm) * dom_perm
    new_dict = {~((k*chopdom)*(~diff_perm)): v for k, v in start_dict.items() if ((k*chopdom)*(~diff_perm)).inv == k.inv + chopdom.inv - diff_perm.inv}
    # ring = DSx([]).ring
    # #first_potato = groth_w0(len(perm), ring, beta=beta)
    # first_potato = dom_groth(dom_perm, ring, beta=beta)
    # schub_elem = ring.zero
    # for perm2, coeff in first_potato.items():
    #     schub_elem += coeff * apply_isobaric_to_schub(diff_perm, perm2, beta=beta)
    # return ring.from_dict({k: v.expand() for k, v in schub_elem.items()})
    print(new_dict)
    return new_dict

def pull_out_groth_var(perm, beta, varnum):
    from schubmult.abc import z
    #dom_perm = Permutation.w0(len(perm))
    ring = DSx([]).ring
    if varnum > perm.max_descent:
        return {perm: S.One}
    work_perm = (~perm)
    n = len(perm)
    dom_perm = Permutation.w0(n)

    cd = (~dom_perm).trimcode[:varnum - 1]
    if varnum  < len((~dom_perm).trimcode):
        dom_chop_cd = (~dom_perm).trimcode[varnum:]
    else:
        dom_chop_cd = []
    new_dom_perm = ~uncode(cd)
    chopdom = ~uncode(dom_chop_cd)
    start_dict = grothmult_double_plus({uncode([1]* (~dom_perm).trimcode[varnum - 1]): 1}, new_dom_perm, ring.coeff_genset, ring.genset[varnum - 1:], beta=beta)
    #Permutation.w0(len(perm))
    diff_perm = (~work_perm) * dom_perm
    print(f"{start_dict=}")
    new_dict = {~((k*chopdom)*(~diff_perm)): v for k, v in start_dict.items() if ((k*chopdom)*(~diff_perm)).inv == k.inv + chopdom.inv - diff_perm.inv}
    
    print(f"{new_dict=}")
    return new_dict


def _check_perm(perm_list):
    import time
    start = time.monotonic()
    perm = Permutation(perm_list)
    ring = DSx([]).ring
    _beta = Gx._beta
    groth1 = alt_grothendieck_poly(perm, beta=_beta)
    # groth2 = grothendieck_poly(perm, ring.genset, ring.coeff_genset, beta=_beta).expand()
    # diff = (groth1.as_polynomial().expand() - groth2).expand()
    
    # if diff != 0:
    #     return (perm_list, f"Grothendieck polynomial for {perm} does not match: {diff=}, \n{groth1=}\n {groth2=}", elapsed)
    # check positivity
    for coeff in groth1.values():
        #coeff = expand(coeff)
        if isinstance(coeff, int):
            if coeff < 0:
                return (perm_list, f"Grothendieck polynomial for {perm} has negative coefficient: {coeff=}, \n{groth1=}", elapsed)
        else:
            result = True
            def walk_expr(args):
                nonlocal result
                for arg in args:
                    if isinstance(arg, int | Integer):
                        if arg < 0:
                            result = False
                            return
                    elif isinstance(arg, Mul | Add | Pow):
                        walk_expr(arg.args)
                    elif isinstance(arg, Symbol):
                        return
                    else:
                        raise ValueError(f"Unexpected type {type(arg)} in expression {arg}")
            walk_expr(coeff.args)
            if not result:
                return (perm_list, f"Grothendieck polynomial for {perm} has negative coefficient", elapsed)
    elapsed = time.monotonic() - start
    return (perm_list, None, elapsed)


if __name__ == "__main__":
    import sys
    from multiprocessing import Pool, cpu_count

    n = int(sys.argv[1])
    num_procs = int(sys.argv[2]) if len(sys.argv) > 2 else cpu_count()
    perms = [perm for perm in Permutation.all_permutations(n) if perm.inv != 0]
    # hardest (highest-inversion) permutations first so a few slow stragglers
    # don't get scheduled last and stall the whole pool at the end
    perms.sort(key=lambda p: p.inv, reverse=True)
    perm_lists = [list(perm) for perm in perms]

    failures = []
    with Pool(num_procs) as pool:
        for i, (perm_list, error, elapsed) in enumerate(pool.imap_unordered(_check_perm, perm_lists), 1):
            if error is not None:
                failures.append(error)
                print(f"FAIL {perm_list}", file=sys.stderr, flush=True)
            print(f"[{i}/{len(perm_lists)}] {perm_list} {'ok' if error is None else 'FAILED'} ({elapsed:.2f}s)", flush=True)

    if failures:
        raise AssertionError("\n".join(failures))
    print(f"All {len(perm_lists)} permutations verified.", flush=True)
    