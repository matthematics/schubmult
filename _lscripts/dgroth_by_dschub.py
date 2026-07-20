from schubmult import *
from schubmult.symbolic.common_polys import grothendieck_poly
from schubmult.symbolic import S, sympify, prod, expand_func
from schubmult.abc import E, y, z, H
from schubmult.symbolic.common_polys import efficient_subs
from functools import cache
from schubmult.symbolic.common_polys import _vars, efficient_subs, elem_func_func_mul, elem_sym_func, elem_sym_poly
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base
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
one_genset = [S.NegativeOne for _ in range(50)]
#negy_dict = {y[i]: -y[i] for i in range(20)}

@cache
def _elem_sym_perms(perm, k):
    return elem_sym_perms(perm, k, k)

def groth_w0(n, beta):
    ring = DSx([]).ring
    start_elem = ring.one
    for i in range(n - 1, 0, -1):
        new_start_elem = 0
        yvar = ring.coeff_genset[n - i]
        bvar = (1 + beta * yvar)
        var2 = [ring.coeff_genset[j] * bvar for j in range(1, 25)]
        
        for permo, coeff in start_elem.items():
            #new_elem = ringa(permo) * E(i, i, ring.genset[1:], ringb.coeff_genset[1:])
            for elem_perm, diff in _elem_sym_perms(permo, i):
                new_start_elem += coeff * (bvar**diff) * prod([var2[permo[p] - 1] + yvar for p in range(i) if permo[p] == elem_perm[p]]) * ring(elem_perm)
            #new_start_elem += coeff * ring.from_dict({k: (1 + beta*ring.coeff_genset[n - i])**((k.inv - permo.inv)) * v for k, v in new_elem.items()})
        # a_j to y_j(1 + y_{n - i})
        start_elem = new_start_elem
        #start_elem = new_start_elem
    return start_elem

# def groth_poly(perm, beta):
    
#     ring = DSx([]).ring
#     # start_elem = ring.one
#     start_dict = {Permutation([]): S.One}
#     n = len(perm)
#     # if perm == Permutation.w0(n):
#     #     return groth_w0(n, beta)
#     for i in range(1, n):
#         new_start_dict = {}
#         for permo, coeff in start_dict.items():
#             new_start_dict = add_perm_dict_with_coeff(grothmult_double_plus({uncode([1] * (n - i)): 1}, permo, ring.genset, ring.coeff_genset[i:], beta=beta), new_start_dict, coeff=coeff)
#         start_dict = new_start_dict
#     return start_dict[(~perm) * Permutation.w0(n)]

def groth_poly(perm, beta, return_dict=False):
    
    ring = DSx([]).ring
    print("Pantacka fart ")
    t = DSx([], "t").ring.coeff_genset
    subs_dict = {t[i]: -ring.genset[i]/(1 + beta * ring.genset[i]) for i in range(20)}
    #return grothmult_double_plus({Permutation.w0(len(perm)): 1}, Permutation([]), [-yy / (1 + beta * yy) for yy in t[:20]], ring.coeff_genset, beta=beta)[perm].subs(subs_dict).expand().simplify().expand()
    if perm.inv == 0:
        return S.One
    n = len(perm)
    #d = min([i for i in range(1, n) if i not in perm.descents()])
    # perm2 = perm.swap(d, d + 1)
    # swp = Permutation.ref_product(d + 1)
    perm2 = Permutation.w0(n)
    swp = (~perm) * Permutation.w0(n)
    # ret = grothmult_double_plus({perm: 1}, Permutation([]), ring.genset, ring.coeff_genset, beta=beta)
    ret = grothmult_double_plus({perm2: 1}, Permutation([]), [0, *[ring.genset[swp[i]] for i in range(20)]], ring.coeff_genset, beta=beta)
    #print({k: v.expand() for k, v in ret.items()})
    if return_dict:
        return {k: v.expand() for k, v in ret.items()}
    return ret[swp]
    # # start_elem = ring.one
    # start_dict = {Permutation([]): S.One}
    # n = len(perm)
    # # if perm == Permutation.w0(n):
    # #     return groth_w0(n, beta)
    # for i in range(1, n):
    #     new_start_dict = {}
    #     for permo, coeff in start_dict.items():
    #         new_start_dict = add_perm_dict_with_coeff(grothmult_double_plus({uncode([1] * (n - i)): 1}, permo, ring.genset, ring.coeff_genset[i:], beta=beta), new_start_dict, coeff=coeff)
    #     start_dict = new_start_dict
    # #return start_dict[perm]
    # return start_dict[(~perm) * Permutation.w0(n)]

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
    dct = grothmult_double({uncode([1] * (n - varnum))}, dom_elem1, ring.coeff_genset,  ring.genset[varnum-1:], beta=beta)
    low_dom_elem = ~uncode(list(range(n - 1 - varnum, 0, -1)))
    new_dct = {k * low_dom_elem: v for k, v in dct.items() if (k * low_dom_elem).inv == k.inv + low_dom_elem.inv}
    new_new_dict = {~(k * (~perm)): v for k, v in new_dct.items() if (k * (~perm)).inv == k.inv - (~perm).inv}
    return new_new_dict
    # for i in range(1, n):
    #     the_var = ring.coeff_genset[i:]
    #     if met:
    #          the_var = ring.coeff_genset[i - 1:]
    #     if i == varnum:
    #         met = True
    #         the_var = z
    #     new_start_dict = {}
    #     for permo, coeff in start_dict.items():
    #         new_start_dict = add_perm_dict_with_coeff(grothmult_double_plus({uncode([1] * (n - i)): 1}, permo, ring.genset, , beta=beta), new_start_dict, coeff=coeff)
    #     start_dict = new_start_dict
    # return start_dict[(~perm) * Permutation.w0(n)]

    #     new_start_elem = 0
    #     yvar = ring.coeff_genset[n - i]
    #     bvar = (1 + beta * yvar)
    #     var2 = [ring.coeff_genset[j] * bvar for j in range(1, 25)]
        
    #     for permo, coeff in start_elem.items():
    #         #new_elem = ringa(permo) * E(i, i, ring.genset[1:], ringb.coeff_genset[1:])
    #         for elem_perm, diff in _elem_sym_perms(permo, i):
    #             new_start_elem += coeff * (bvar**diff) * prod([var2[permo[p] - 1] + yvar for p in range(i) if permo[p] == elem_perm[p]]) * ring(elem_perm)
    #         #new_start_elem += coeff * ring.from_dict({k: (1 + beta*ring.coeff_genset[n - i])**((k.inv - permo.inv)) * v for k, v in new_elem.items()})
    #     # a_j to y_j(1 + y_{n - i})
    #     start_elem = new_start_elem
    #     #start_elem = new_start_elem
    # return start_elem


def isobaric_strip_on_dschub(start, length, schub_perm):
    start_schub = DSx(schub_perm)
    
    ring = start_schub.ring
    var_strip = ring.genset[start + 1:start + length + 1]
    bigger_schub = E(length, length, var_strip, one_genset) * start_schub
    stripness = list(range(start, start + length))
    ret_schub = ring.from_dict({k: v for k, v in bigger_schub.items()})
    for desc in stripness:
        ret_schub = ring.from_dict({perm2.swap(desc - 1, desc): v for perm2, v in ret_schub.items() if perm2[desc - 1] > perm2[desc]})
    return ret_schub

@cache
def apply_isobaric_to_schub(diff_perm, schub_perm):
    elem = DSx(schub_perm)
    
    strips = [[i, diff_perm.trimcode[i - 1]] for i in range(1, diff_perm.max_descent + 1)]
    for strip in reversed(strips):
        if strip[1] == 0:
            continue
        new_elem = elem.ring.zero
        for perm, coeff in elem.items():
            new_elem += coeff * isobaric_strip_on_dschub(strip[0], strip[1], perm)
        elem = new_elem
    return elem

# def dom_groth_factor(dom, ring):
#     from schubmult.rings.schubert.double_schubert_ring import DoubleSchubertRing
#     from schubmult.symbolic.poly.variables import CustomGeneratingSet
#     one_set = CustomGeneratingSet([1 for _ in range(50)])
#     ybeta_set = CustomGeneratingSet([-beta*ring.coeff_genset[i] for i in range(50)])
#     pants_ring = DoubleSchubertRing(one_set, ybeta_set)
    
#     return pants_ring(dom).as_polynomial()

def alt_grothendieck_poly(perm):
    from schubmult.abc import z
    dom = perm.minimal_dominant_above()
    #Permutation.w0(len(perm))
    diff_perm = (~perm) * dom
    ring = DSx([]).ring
    first_potato = ring.from_expr(dom_groth(dom, ring.genset, ring.coeff_genset, beta))
    schub_elem = ring.zero
    for perm2, coeff in first_potato.items():
        schub_elem += coeff * apply_isobaric_to_schub(diff_perm, perm2)
    fat_scubs = {y[i]: z[i] - 1 for i in range(1, 50)}
    return ring.from_dict({k: v.expand() for k, v in schub_elem.items()})


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    ring = DSx([]).ring
    # for perm in perms:
    #     poly1 = alt_grothendieck_poly(perm)
            
    #     print(f"Grothendieck polynomial for {perm} verified.")
    #     print(poly1)
    for perm in perms:
        mx = 2
        if perm.inv == 0:
            continue
        for varnum in  range(1, mx):
            print(f"Trying {perm} {varnum}")
            #groth1 = groth_poly(perm, 1)
            #genset = MaskedGeneratingSet(ring.coeff_genset, index_mask = [varnum])
            groth1 = groth_poly(perm, beta=Gx._beta, return_dict=False)
            print(f"Fat computer")
            good = True
            #for perm2, coeff in groth1_dct.items():
            #    groth1 = coeff
            groth2 = grothendieck_poly(perm, ring.genset, ring.coeff_genset, beta=Gx._beta).expand()
            diff = (groth1 - groth2).expand()
            try:
                assert diff == 0, f"Grothendieck polynomial for {perm} does not match: {diff=}, \n{groth1=}\n {groth2=}"
            except AssertionError as e:
                good = False
            if good:
                print("Happy pants")
            else:
                print("Sad pants")