from schubmult import *
from schubmult.symbolic.symmetric_polynomials.elem_sym import FactorialElemSym
from schubmult.rings.polynomial_algebra import *
#from schubmult.symbolic import expand
from sympy import pretty_print, Add, Mul, Pow, sympify, expand, S

x = GeneratingSet("x")
t = GeneratingSet("t")

r = BoundedRCFactorAlgebra()

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

def elem_sym_as_rc_sum(elem_sym_elem, size):
    p = elem_sym_elem.degree
    k = elem_sym_elem.numvars
    result = 0
    for subb in itertools.combinations(range(1, k+1), p):
        weight = [1 if i in subb else 0 for i in range(1, k + 1)]
        result += as_double_elem_sym_rc(weight, elem_sym_elem.coeffvars, size)
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


if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    #length = int(sys.argv[2])

    #comps = [comp for comp in itertools.combinations_with_replacement(range(k), length)]
    perms = [perm for perm in Permutation.all_permutations(n) if perm.inv != 0]

    # DForest = PolynomialAlgebra(DoubleForestPolyBasis(x, t))
    dring = DSx([], "t").ring

    # for p in range(1, k + 1):
    #     for subb in itertools.combinations(range(1, k+1), p):
    #         weight = [1 if i in subb else 0 for i in range(1, k + 1)]
    #         pretty_print(as_double_elem_sym_rc(weight, t))
    size = 2 * n
    for perm1, perm2 in itertools.product(perms, repeat=2):
        poly1 = dring(perm1)
        poly2 = dring(perm2)
        base1 = poly1.in_SEM_basis().expand(deep=True)
        base2 = poly2.in_SEM_basis().expand(deep=True)

        factor_elem1 = symbol_expr_as_rc_sum(base1, size)
        factor_elem2 = symbol_expr_as_rc_sum(base2, size)

        dothis = (factor_elem1 * factor_elem2).to_rc_graph_ring_element()
        dothis = dothis.ring.from_dict({k.resize(n - 1): v for k, v in dothis.items() if expand(v, deep=True) != 0})
        
        actual_result = poly1 * poly2
        
        # for perm, coeff2 in actual_result.items():
        #     # tump = sum([expand(v, deep=True) * Sx.from_expr(rc.polyvalue(t)) for  rc, v in dothis.divdiff_perm(perm).items()])
        #     # coeff = tump.get(Permutation([]), S.Zero)
        #     # diffdiff = (coeff2 - coeff).expand(deep=True)

        while any(expand(v, deep=True) !=0 for v in dothis.values()):
            dothis2 = dothis.ring.from_dict(dothis)
            top_term = RCGraph.principal_rc(*max([(k.perm, len(k)) for k,v in dothis2.items() if expand(v, deep=True) != 0], key=lambda rc: (rc[0].inv, rc[0])))

            coeff = dothis2[top_term].expand(deep=True)
            
            diffdiff = (coeff - actual_result.get(top_term.perm, 0)).expand(deep=True)
            if diffdiff != 0:
                print(f"Mismatch for {perm1} * {perm2} at {top_term.perm}: got {coeff}, expected {actual_result.get(top_term.perm, 0)}")
                break
            else:
                print(f"Match for {perm1} * {perm2} at {top_term.perm}: got {coeff}, expected {actual_result.get(top_term.perm, 0)}")
            dothis3 = dothis2 - coeff * symbol_expr_as_rc_sum(dring(top_term.perm).in_SEM_basis().expand(deep=True), size).to_rc_graph_ring_element().resize(n - 1)
            dothis = dothis3
            
        # for rc, coeff in dothis.items():
        #     diffdiff = (coeff - actual_result.get(rc.perm, 0)).expand()
        #     if diffdiff != 0:
        #         print(f"Mismatch for {perm1} * {perm2} at {rc}: got {coeff}, expected {actual_result.get(rc.perm, 0)}")
        #     else:
        #         print(f"Match for {perm1} * {perm2} at {rc}: got {coeff}, expected {actual_result.get(rc.perm, 0)}")


        
    

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
    
        
