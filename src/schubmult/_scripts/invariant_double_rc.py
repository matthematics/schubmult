from schubmult import *
from schubmult.symbolic.symmetric_polynomials.elem_sym import FactorialElemSym
from schubmult.rings.polynomial_algebra import *
#from schubmult.symbolic import expand
from sympy import pretty_print, Add, Mul, Pow, sympify, expand, S
import numpy as np
x = GeneratingSet("x")
t = GeneratingSet("t")

r = BoundedRCFactorAlgebra()
rr = RCGraphRing()

def as_elem_sym_rc(weight, size):
    if not set(weight).issubset({0,1}):
        raise ValueError("Weight must be binary")
    indexes = [i + 1 for i in range(len(weight)) if weight[i] == 1]
    k = len(weight)
    p = sum(weight)
    result = 1
    for spot, ind in enumerate(indexes, start=1):
        result *= (r(r.make_key((RCGraph.monk_rc(ind, k),),size)))
    return result

def as_double_elem_sym_rc(weight, genset, size):
    if not set(weight).issubset({0,1}):
        raise ValueError("Weight must be binary")
    indexes = [i + 1 for i in range(len(weight)) if weight[i] == 1]
    k = len(weight)
    p = sum(weight)
    result = 1
    for spot, ind in enumerate(indexes, start=1):
        result *= (r(r.make_key((RCGraph.monk_rc(ind, k),),size)) - genset[k - p - (ind - spot)])
    return result

def elem_sym_as_rc_sum(elem_sym_elem, size, retvars=False):
    p = elem_sym_elem.degree
    k = elem_sym_elem.numvars
    result = 0
    vardict = {}
    for subb in itertools.combinations(range(1, k+1), p):
        weight = [1 if i in subb else 0 for i in range(1, k + 1)]
        result += as_double_elem_sym_rc(weight, elem_sym_elem.coeffvars, size)
        vardict[((k,tuple(weight)),)] = elem_sym_elem.coeffvars
    if retvars:
        return result, vardict
    return result

def symbol_expr_as_rc_sum(symbol_expr, size):
    result = 0
    symbol_expr = sympify(symbol_expr).expand(deep=True)
    if isinstance(symbol_expr, Add):
        for term in symbol_expr.args:
            result += symbol_expr_as_rc_sum(term, size)
        return result
    if isinstance(symbol_expr, Mul):
        result = 1
        for term in symbol_expr.args:
            result *= symbol_expr_as_rc_sum(term, size)
        return result
    if isinstance(symbol_expr, Pow):
        
        base_sum = symbol_expr_as_rc_sum(symbol_expr.base, size)
        result = 1
        for _ in range(term.exp):
            result *= base_sum
        return result
    if isinstance(symbol_expr, FactorialElemSym):
        return elem_sym_as_rc_sum(symbol_expr, size)
    return symbol_expr

def symbol_expr_as_weight_coeff_dict(symbol_expr, genset, size):
    symbol_expr = sympify(symbol_expr).expand(deep=True)
    if isinstance(symbol_expr, Add):
        result = {}
        for term in symbol_expr.args:
            result.update({k: v + result.get(k, 0) for k, v in symbol_expr_as_weight_coeff_dict(term, genset, size).items()})
        return result
    if isinstance(symbol_expr, Mul):
        new_result = {(): 1}
        for term in symbol_expr.args:
            new_new_result = {}
            old_result = symbol_expr_as_weight_coeff_dict(term, genset, size)
            for k, v in old_result.items():
                for k2, v2 in new_result.items():
                    key = tuple(sorted(k + k2))
                    new_new_result[key] =  new_new_result.get(key, 0) + v * v2
            new_result = new_new_result
        return new_result
    if isinstance(symbol_expr, Pow):
        
        base_sum = symbol_expr_as_weight_coeff_dict(symbol_expr.base, genset, size)
        result = {(): 1}
        for _ in range(symbol_expr.exp):
            new_result = {}
            for k, v in result.items():
                for k2, v2 in base_sum.items():
                    key = tuple(sorted(k + k2))
                    new_result[key] = new_result.get(key, 0) + v * v2
            result = new_result
        return result
    if isinstance(symbol_expr, FactorialElemSym):
        _, dct = elem_sym_as_rc_sum(symbol_expr, size, retvars=True)
        result = {}
        for ((k, weight),), vars in dct.items():
            monom_weight = [0] * size
            var_tuple = tuple(sorted([genset.index(yy) for yy in vars]))
            indexes = [i + 1 for i in range(len(weight)) if weight[i] == 1]
            for spot, index in enumerate(indexes, start=1):
                monom_weight[var_tuple[(index - spot)] - 1] += 1
            result[((k,weight,tuple(monom_weight)),)] = result.get(((k,weight,tuple(monom_weight)),), 0) + 1
        return result
    return {(): symbol_expr}


def double_monk_squash(rc, monk_rc, genset1, genset2 = None):
    if genset2 is None:
        genset2 = genset1
    if isinstance(rc, RCGraph):
        rc2 = rc.squash_product(monk_rc)
        monk_row, monk_col = monk_rc.left_to_right_inversion_coords(0)
        values = sorted([rc.perm[i] for i in range(max(monk_rc.perm.descents()) + 1)])
        result = rr(rc2)
        for mrc in RCGraph.all_rc_graphs(monk_rc.perm):
            result += (mrc.polyvalue(genset1, genset2)) * rr(rc)
        return result
    result = rr.zero
    for rcc, coeff in rc.items():
        result += coeff * double_monk_squash(rcc, monk_rc, genset1, genset2)
    return result



if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    #length = int(sys.argv[2])
    y = GeneratingSet("y")
    t = GeneratingSet("t")
    #comps = [comp for comp in itertools.combinations_with_replacement(range(k), length)]
    perms = [perm for perm in Permutation.all_permutations(n) if perm.inv ==0]#!= 0 and perm.descents() == {k - 1}]

    for perm in perms:
        actual = DSx(perm) * DSx(Permutation.ref_product(*tuple(range(1,k + 1))), "t")
        base = rr.schub(perm, k).resize(k)
        #result = rr.zero
        result = base
        for i in range(1, k + 1):
            monk_rc = RCGraph.monk_rc(i, k)
            #result += double_monk_squash(base, monk_rc, y, t)
            result = double_monk_squash(result, monk_rc, y, [t[1]] * 50)
        #new_result = double_monk_squash(result, )
        failed = False
        for rcc, coeff in result.items():
            coeff2 = actual.get(rcc.perm, 0)
            if (coeff - coeff2).expand() != 0:
                print(f"Mismatch for {perm} term {rcc}: got {coeff}, expected {coeff2}")
                print(f"Difference: {(coeff - coeff2).expand()}")
                print(f"Actual: {actual.get(rcc.perm, 0)}")
                print(f"Result: {result[rcc]}")
                input("Press Enter to continue...")
                failed = True
        if not failed:
            print(f"Match for {perm}")

    # DForest = PolynomialAlgebra(DoubleForestPolyBasis(x, t))
    # dring = DSx([], "t").ring

    # for p in range(1, k + 1):
    #     for subb in itertools.combinations(range(1, k+1), p):
    #         weight = [1 if i in subb else 0 for i in range(1, k + 1)]
    #         pretty_print(as_double_elem_sym_rc(weight, t))
    # size = 2 * n
    # #for perm1, perm2 in itertools.product(perms, repeat=2):
    # for perm in perms:
    #     poly1 = dring(perm)
    #     # poly2 = dring(perm2)
    #     base1 = poly1.in_SEM_basis().expand(deep=True)
    #     # base2 = poly2.in_SEM_basis().expand(deep=True)
    #     rc_rep = r.from_dict({key: v for key, v in symbol_expr_as_rc_sum(base1, size).items() if r.key_to_rc_graph(key).resize(n-1).forest_weight == r.key_to_rc_graph(key).resize(n-1).perm.pad_code(n-1)})
    #     bacon_test = DForest.from_expr(rc_rep.to_rc_graph_ring_element().polyvalue(x), length=n-1)
    #     assert bacon_test.almosteq(DForest(*perm.pad_code(n-1))), f"Mismatch for {perm}: got {bacon_test}, expected {DForest(*perm.pad_code(n-1))}"
    #     if False:
    #         factor_elem1 = symbol_expr_as_weight_coeff_dict(base1, t, n-1)
    #         print(factor_elem1)
    #         for tup, coeff in factor_elem1.items():
    #             rc = 1
    #             total_monom_weight = [0] * (n - 1)
    #             path_weight = []
    #             check_k = 1
    #             vanishes = False
    #             for k, weight, monom_weight in tup:
    #                 while check_k < k:
    #                     path_weight = [check_k, *path_weight]
    #                     check_k += 1
    #                 path_weight = [k - sum(weight), *path_weight]
    #                 #for index in var_tuple:
    #                 # indexes = [i + 1 for i in range(len(weight)) if weight[i] == 1]
    #                 # for spot, index in enumerate(indexes, start=1):
    #                 #     monom_weight[var_tuple[k - sum(weight) - (index - spot)] - 1] += 1
    #                 total_monom_weight = (np.array(total_monom_weight) + np.array(monom_weight)).tolist()
    #                 if any((monom_weight[i] > 0 and weight[i] == 1) for i in range(len(weight))):
    #                     vanishes = True
    #                 rc *= as_elem_sym_rc(weight, size)
    #                 check_k += 1
    #             print(f"{tup=}")
    #             print(f"Monomial with weight {total_monom_weight} has coefficient {coeff} and\n{rc.to_rc_graph_ring_element().resize(n-1)=}")
    #             monom_coeff = Schub(~perm, n - 1).change_basis(MonomialBasis).get(tuple(total_monom_weight), 0)
    #             path_coeff = Schub(perm * Permutation.w0(n), n - 1).change_basis(MonomialBasis).get(tuple(path_weight), 0)
    #             print(f"coeff in v^{{-1}}: {monom_coeff}")
    #             print(f"coeff in vw_0: {path_coeff}")
    #             rc_graph = next(iter(rc.to_rc_graph_ring_element()))
    #             rc_is = rc_graph.perm == perm
    #             print(f"{rc_is=}")
    #             print(f"{path_weight=}")
    #             print(f"column weight {rc_graph.transpose(n-1).length_vector}")
    #             print(f"{vanishes=}")
    #             if not vanishes:
    #                 assert not rc_is
    #                 print("Paint")
    #                 input()
    #         # if rc_is:
    #         #     assert monom_coeff != 0, f"Mismatch for {perm} monomial weight {total_monom_weight}: got {rc_is} but expected nonzero coeff in v^-1"
    #         #     print("Gotthatright")
    #         # if rc_is:
    #         #     assert monom_coeff * coeff > 0, f"Mismatch for {perm} monomial weight {monom_weight}: got {rc_is} but expected {monom_coeff*coeff > 0}"
    #         #assert rc_is == (monom_coeff*coeff > 0), f"Mismatch for {perm} monomial weight {monom_weight}: got {rc_is} but expected {monom_coeff*coeff > 0}"
    #     # factor_elem2 = symbol_expr_as_rc_sum(base2, size)

    #     # dothis = (factor_elem1 * factor_elem2).to_rc_graph_ring_element()
    #     # dothis = dothis.ring.from_dict({k.resize(n - 1): v for k, v in dothis.items() if expand(v, deep=True) != 0})
        
    #     # actual_result = poly1 * poly2
        
    #     # # for perm, coeff2 in actual_result.items():
    #     # #     # tump = sum([expand(v, deep=True) * Sx.from_expr(rc.polyvalue(t)) for  rc, v in dothis.divdiff_perm(perm).items()])
    #     # #     # coeff = tump.get(Permutation([]), S.Zero)
    #     # #     # diffdiff = (coeff2 - coeff).expand(deep=True)

    #     # while any(expand(v, deep=True) !=0 for v in dothis.values()):
    #     #     dothis2 = dothis.ring.from_dict(dothis)
    #     #     top_term = RCGraph.principal_rc(*max([(k.perm, len(k)) for k,v in dothis2.items() if expand(v, deep=True) != 0], key=lambda rc: (rc[0].inv, rc[0])))

    #     #     coeff = dothis2[top_term].expand(deep=True)
            
    #     #     diffdiff = (coeff - actual_result.get(top_term.perm, 0)).expand(deep=True)
    #     #     if diffdiff != 0:
    #     #         print(f"Mismatch for {perm1} * {perm2} at {top_term.perm}: got {coeff}, expected {actual_result.get(top_term.perm, 0)}")
    #     #         break
    #     #     else:
    #     #         print(f"Match for {perm1} * {perm2} at {top_term.perm}: got {coeff}, expected {actual_result.get(top_term.perm, 0)}")
    #     #     dothis3 = dothis2 - coeff * symbol_expr_as_rc_sum(dring(top_term.perm).in_SEM_basis().expand(deep=True), size).to_rc_graph_ring_element().resize(n - 1)
    #     #     dothis = dothis3
            
    #     # for rc, coeff in dothis.items():
    #     #     diffdiff = (coeff - actual_result.get(rc.perm, 0)).expand()
    #     #     if diffdiff != 0:
    #     #         print(f"Mismatch for {perm1} * {perm2} at {rc}: got {coeff}, expected {actual_result.get(rc.perm, 0)}")
    #     #     else:
    #     #         print(f"Match for {perm1} * {perm2} at {rc}: got {coeff}, expected {actual_result.get(rc.perm, 0)}")


        
    

    # for comp in comps:
    #     perm = uncode(comp)
    #     if not perm.is_dominant:
    #         continue
    #     #tst = DForest(*comp).expand()
    #     tst = DSx(perm, "t").expand()
    #     dfor = DForest.from_expr(tst, length=length)
    #     print(f"{comp}: {dfor}")
    #grass_perms = [perm for perm in perms if perm.descents() == {k - 1}]
    #check_dict = {}
    # t = GeneratingSet("t")
    # r = RCGraphRing()
    # for seq in itertools.product(range(1, k + 1), repeat=length):
    #     print(f"{seq=}")
    #     the_rc = RCGraph([()]).resize(k)
    #     the_rc_test = r(RCGraph([()]).resize(k))
    #     for a in list(seq):
    #         old_the_rc = the_rc
    #         the_rc = the_rc.squash_product(RCGraph.monk_rc(a, k))
    #         coldiff = set([b for row in range(len(the_rc)) for b in the_rc[row]]).difference([b for row in range(len(old_the_rc)) for b in old_the_rc[row]])
    #         if len(coldiff) > 1:
    #             raise RuntimeError(f"Unexpectedly got more than one new crossing: {coldiff} for seq={seq} at a={a}")
    #         #the_rc_test = the_rc_test.broadcast._double_squash_monk(a, k, t, t[next(iter(coldiff)) + 1 - a])
    #         the_rc_test = the_rc_test.broadcast.squash_product(RCGraph.monk_rc(a, k)) - t[next(iter(coldiff))]*the_rc_test
    #     print(repr(the_rc))
    #     if the_rc in check_dict:
    #         assert the_rc_test.almosteq(check_dict[the_rc]), f"mismatch for seq={seq} the_rc={the_rc} test={the_rc_test} expected={check_dict[the_rc]}"
    #         print(f"Good: ", seq, f"{dict(the_rc_test)=}")
    #     else:
    #         check_dict[the_rc] = the_rc_test
    
        
