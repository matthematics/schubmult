from schubmult import *
from schubmult.rings import *
from schubmult.symbolic import S
from itertools import permutations
import sys
from functools import cache
import numpy as np
import numpy.linalg as nla

perm = uncode([2, 2, 1])


@cache
def generate(perm):
    val = ASx(perm)

    expr = []

    while val != val.ring.zero:
        mx = [k[0].code for k in val.keys() if val[k] != S.Zero ]
        mx.sort(reverse=True)
        cd = mx[0]
        print(cd)
        mx_key = next(iter([k for k in val.keys() if k[0].code == cd]))
        cd = [*cd]
        fv = cd.pop(0)
        while len(cd) > 1 and cd[-1] == 0:
            cd.pop()
        expr2 = [[1,cd]]
        cf = val[mx_key]
        val -= cf * ASx(uncode([fv]),1)*ASx(uncode(cd))
        if len(cd) > 1:
            expr2 = generate(uncode(cd))

        for dingbats in expr2:
            cf2 = cf * dingbats[0]
            toappend = [fv, *dingbats[1]]
            expr += [[cf2, toappend]]

        print(val)
        
    return expr

@cache
def single_coprod(p, n, T):
    res = T.zero
    for i in range(p+1):
        res += T.from_dict({((uncode([i]), n), (uncode([p-i]), n)): S.One})
    return res

def coprod(val):
    R = (ASx([]).ring @ ASx([]).ring)
    res = R.zero

    for (perm, n), v in val.items():
        res += v * FA([]).ring.tensor_schub_expand(schub_elem(perm, n).coproduct())
    return res


    # if val == val.ring.zero:
    #     return T.zero
    
    # cprd_val = T.zero

    # while val != val.ring.zero:
    #     mx = [k[0].code for k in val.keys() if val[k] != S.Zero ]
    #     mx.sort(reverse=True)
    #     cd = mx[0]
        

    #     mx_key = next(iter([k for k in val.keys() if k[0].code == cd]))
    #     if len(cd) == 0:
    #         return cprd_val + T.from_dict({((Permutation([]),mx_key[1]),(Permutation([]),mx_key[1])): val[mx_key] * S.One})
    #     cd = [*cd]
    #     fv = cd.pop(0)
    #     while len(cd) > 1 and cd[-1] == 0:
    #         cd.pop()
    #     cf = val[mx_key]
    #     cprd_val += (T.from_dict({((Permutation([]),0),(Permutation([]),0)): cf*S.One}))* single_coprod(fv, 1, T) * coprod(ASx(uncode(cd), mx_key[1] - 1), T)
    #     val -= cf * ASx(uncode([fv]),1)*ASx(uncode(cd), mx_key[1] - 1)
    #     print(val)
    
    return cprd_val

def split_hom(a, b, p, m):
    res = a.ring.zero
    for k, v in a.items():
        for k2, v2 in b.items():
            try:
                assert (len(k) <= p and len(k2) == 0) or (len(k) == p and len(k2) <= m) or (len(k) >= p and len(k2) == m)
            except AssertionError:
                print(f"{k=}, {k2=}, {p=}, {m=}")
                raise
            if len(k) < p:
                res += v2* v * FA(k)
            elif len(k) == p:
                res += v * v2 * FA(k) * FA(k2)
            else:
                res += v * v2 * FA(k[:p])*FA(k2)*FA(k[p:])
    return res

def schub_split_hom(a, b, p, m):
    res = ASx([]).ring.zero
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            res += v1*v2*split_hom(FA([]).ring.schub_elem(*k1), FA([]).ring.schub_elem(*k2), p, m).schub_expand()
    return res


from schubmult.symbolic import *
    

def schub_elem(perm, numvars):
    res = FA([]).ring.zero
    expr = Sx(perm*~uncode(list(range(perm.inv + numvars, perm.inv, -1)))).in_SEM_basis().expand()
    #print(f"{expr=}")   
    args = expr.args
    if not isinstance(expr, Add):
        args = [expr]
    for arg in args:
        tup = list(range(perm.inv + numvars, perm.inv, -1))
        coeff = S.One
        #print(f"{arg=}")
        if is_of_func_type(sympify(arg), FactorialElemSym):
            arg = sympify_sympy(arg)
            tup[perm.inv + numvars - arg.numvars] = arg.numvars - arg.degree
        elif isinstance(arg, Mul):
            for arg0 in arg.args:
                if is_of_func_type(sympify(arg0), FactorialElemSym):
                    arg0 = sympify_sympy(arg0)
                    tup[perm.inv + numvars - arg0.numvars] = arg0.numvars - arg0.degree
                else:
                    coeff = Integer(arg0)
        else:
            coeff = Integer(arg)
        res += coeff * FA(tup)
    return res

def varlen(perm):
    if perm.inv == 0:
        return 0
    return max(perm.descents()) + 1

if __name__ == "__main__":
    n0 = 5
    perms = Permutation.all_permutations(n0)

    
    #perms2 = Permutation.all_permutations(6)
    

    from schubmult.abc import *
    from schubmult.perm_lib import *
    Permutation.print_as_code = True
    numvars = 3
    for perm in perms:
        print(f"{perm}: {coprod(ASx(perm))}")    


        # fa = FA([]).ring.schub_elem(perm, varlen(perm))
        # res = Sx([]).ring.zero
        # for burger, v in fa.items():
        #     to_add = Sx([])
        #     for letter in burger:
        #         if letter > 0:
        #             to_add *= h(letter, n0, x[1:])
        #     res += v*to_add
        # print(f"{perm.code}: {fa}")
        # print(f"{perm.code}: {res}")

    sys.exit(0)

    numvars0 = 2
    numvars1 = 2
    p = 1
    m = 2
    pad = 2
    # tups = [trimcode(perm2) for perm2 in perms2 if len(trimcode(perm2)) <= numvars]
    # tups2 = []
    # for tup in tups:
    #     tups2 += [(*tup,*([0] * (numvars - len(tup))))]
    # tups = tups2
    results = {}
    results2 = {}
    bingo = set(perms)
    Permutation.print_as_code = True
    for perm3 in perms:
        if varlen(perm3) != 5:
            continue
        cprd = Sx(perm3).coproduct(1,2, 5)
        for (perm, perm2), v in cprd.items():
            if varlen(perm) != 2:
                continue
            if varlen(perm2) > 2:
                continue
            boing = (ASx(perm, 2) * ASx(perm2, 2) * ASx([], 1))
            frc = boing[(perm3,5)]/((ASx(perm, 2) * ASx([],1))[(perm, 3)])
            #print(frc)
            boing0 = schub_split_hom(ASx(perm, 3), ASx(perm2,2), 2,2)
            # print(f"{boing0[(perm3,5)]=}")
            # print(f"{boing[(perm3,5)]=}")
            try:
                assert boing0[(perm3,5)] == boing[(perm3,5)]
            except AssertionError:
                print("Nope")
            #     print(boing0)
            #     print(boing)
                print(f"{perm.code=} {perm2.code=} {perm3=}")
            #print(f"{boing0[(perm3,5)]} != {frc}")
            #    sys.exit(1)

    #     twotwo = ASx([],1) * ASx(perm2)
    #     a1 = ASx(perm1) * twotwo
    #     a2 = schub_split_hom(ASx(perm1) * ASx(perm2), ASx([], 1), 2, 1)
    # try:
    #     assert a1 == a2
    # except AssertionError:
    #     print(f"Nope")
    #     print(f"{a1}")
    #     print(f"{a2}")    

        # split_at_m = Sx(perm).coproduct(1,2,3)
        # for (perm0, perm1), v in split_at_m.items():
        #     for perm2 in perm2:
        #         lower_perm = perm0 * (~perm2)
        #         if lower_perm.inv == perm0.inv - perm2.inv:
        #             split_at_p = Sx(lower_perm).coproduct(1)
        #             for (perm00, perm11), v2 in split_at_p.items():
        #                 if perm11.inv == 0
        # if varlen(perm) != numvars0 + numvars1:
        #     continue
        # for perm2 in perms:
        #     if varlen(perm2) != 1:
        #         continue
        #     lower_perm = (perm * (~perm2))
        #     if lower_perm.inv != perm.inv - perm2.inv:
        #         continue
        #     splitter = Sx(lower_perm).coproduct(1)
        #     for (perm0, perm00), v in splitter.items():
        #         if perm00.inv == 0:
        #             res1 = (v * split_hom(FA([]).ring.schub_elem(lower_perm,numvars0), FA([]).ring.schub_elem(ID_PERM,numvars1),p,m)* ).schub_expand()
        #             print(f"{perm} @ {perm2}: {res}")
        #             for perm3, v in res.items():
        #                 assert v == Sx(perm3[0]).coproduct(1,4)[(perm, perm2)]
        # if perm.inv>0 and max(perm.descents()) + 1 > numvars:
        #     continue
        # delem = FA([]).ring.zero
        # for perm2 in perms:
        #     if (perm2*(~perm)).inv == perm2.inv - perm.inv and (perm2.inv == 0 or max(perm2.descents()) < numvars):
        #         coeff = efficient_subs(Sx(~((perm2)*(~perm))).as_polynomial(),dict([(xx, yy) for (xx,yy) in zip(Sx([]).ring.genset,[-yy for yy in DSx([]).ring.coeff_genset])]))
        #         #print(f"{coeff=}")
        #         delem += coeff*delem.ring.schub_elem(perm2, numvars)
        #         delem = delem.ring.from_dict({k: v for k,v in delem.items() if expand(v) != S.Zero})
        # print(f"{perm}: {delem}")
        # results[perm] = delem
        #print(f"{perm}: {split_hom(FA([]).ring.schub_elem)} ")

    # for perm1 in perms:
    #     if perm1 not in results:
    #         continue
    #     for perm2 in perms:
    #         if perm2 not in results:
    #             continue
    #         frof = results[perm1] * results[perm2]
    #         kk = list(frof.keys())
    #         for key in kk:
    #             if key[0] not in results:
    #                 del frof[key]
    #         binf = ASx([]).ring.zero
    #         boing0 = ADSx(perm1,numvars) * ADSx(perm2,numvars)
    #         for k, v in boing0.items():
    #             if k[0] in results:
    #                 binf += v * ASx([]).ring.from_dict({(kkk[0], numvars*2): vvv for kkk, vvv in results[k[0]].items()})
    #         #print(f"{binf=}")
    #         #print(f"{frof=}")
    #         flerf = binf - frof
    #         for k, v in flerf.items():
    #             try:
    #                 assert expand(v) == S.Zero
    #                 print(f"success {1}")
    #             except AssertionError:
    #                 print(f"fail {1}")
                    # print(f"{k=}")
                    # print(f"{binf=}")
                    # print(f"{frof=}")


    sys.exit(0)
                
    #     res = ASx([]).ring.zero
    #     if perm.inv == 0 or max(perm.descents()) < numvars:
    #         spoingle = Sx(DSx(perm).expand())
    #         for k, v in spoingle.items():
    #             if k.inv == 0 or max(k.descents()) < numvars:
    #                 results[k] = results.get(k, res) + v*ASx(perm,numvars)

    # for perm in perms:                
    #     res = ADSx([]).ring.zero
    #     if perm.inv == 0 or max(perm.descents()) < numvars:
    #         spoingle = DSx(Sx(perm).expand())
    #         for k, v in spoingle.items():
    #             if k.inv == 0 or max(k.descents()) < numvars:
    #                 results2[k] = results2.get(k, res) + v*ADSx(perm,numvars)

        # for tup in tups:
        #     #print(tup)
        #     mulval = DSx([])    
        #     #print(f"{tup=}")
        #     for i, deg in enumerate(tup):                
        #         if deg > 0:
        #             mulval *= H(deg, 1, [x[i+1]], [0,0,0,0,0,0,0,0,0,0,0] )
        #     #print(mulval)

        #     res += mulval.get(perm, S.Zero) * FA(tup)
        #     #print(res)      
        # print(f"{perm}: {res}")
        # voible = DSx([]).ring.zero
        # for k, v in res.items():
        #     fiff = ASx([])
        #     for iff in reversed(k):
        #         fiff = (ASx(uncode([iff]), 1) * fiff)
        #     voible += v * Sx([]).ring.from_dict({k2[0]: v2 for k2, v2 in fiff.items()})
    #     print(f"result: {voible}")
        # assert res == FA([]).ring.schub_elem(perm, numvars)
    Permutation.print_as_code = True
    a = results[uncode([1,1])]
    b = results[uncode([2,1])]
    res =  a * b 
    print(f"{res=}")
    rois = ADSx([])
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            rois += v1*v2*results2[k1[0]]*results2[k2[0]]
    print(f"{rois=}")

    sys.exit(0)

    perm_index = {}
    for i, perm in enumerate(perms):
        perm_index[perm] = i

    def conv(f, g, elem):
        blorg = FA(elem).coproduct()
        res = FA([]).ring.zero
        for (key1, key2), val in blorg.items():
            res += val * f(key1)*g(key2,len(elem))
        return res

    def idfunc(b):
        return FA(b)
    
    elem_index = {}

    gensets = {}
    def gfunc(key, n):
        global perms, perm_index, elem_index, gensets
        
        res = FA([]).ring.zero
        for perm in perms:
            xp = GeneratingSet(f"x{perm_index[perm]}")
            gensets[perm] = xp
            the_key = perm.code
            while len(the_key)>0 and the_key[-1] == 0:
                the_key.pop()
            m = len(the_key)
            
            
            for i in range(m,n):
                
                if i == m:
                    key0 = (*the_key,)
                else:
                    key0 = (*the_key,*([0]*(n-m)))
                if key0 not in elem_index:
                    elem_index[key0] = len(elem_index.keys())
                res += xp[elem_index[key0]] * FA(key0)
        return res
    
    elem = (2, 0, 1, 3)

    
    for perm in perms:
        res_elem = FA([]).ring.zero
        S = conv(idfunc, gfunc, (*perm.code,))
        print(f"S({perm.code}) = {S}")
        # for key, val in S.items():
        #     if key == ():
        #         print(f"{key} {val}=1")
        #         for perm in gensets:
        #             L = [elem0 for elem0 in elem_index.keys() if elem_index[elem0] == gensets[perm].index(val)]
        #             if len(L) > 0:
        #                 res_elem += sum([FA(a) for a in L])
        #     else:
        #         print(f"{key} {val}=0")
    
    





    sys.exit(0)
    matrix_rows = {}
    perm_index = {}
    for i, perm in enumerate(perms):
        perm_index[perm] = i
        matrix_rows[perm] = Sx(perm).to_genset_dict(trim=True)
    
    poly_index = {}
    num_poly = 0
    matrix = np.ndarray((len(perms), len(perms)), dtype = np.int64)
    matrix.fill(0)
    for perm, row in matrix_rows.items():
        for giblet, val in row.items():
            if giblet not in poly_index:
                poly_index[giblet] = num_poly
                num_poly += 1
            print(poly_index)
            matrix[perm_index[perm], poly_index[giblet]] = val
    print(matrix)
    invmatrix = nla.inv(matrix).astype(np.int64)
    print(invmatrix)

    print(perm_index)
    print(poly_index)


    import sys
    sys.exit(0)
    expr = generate(perm)

    res = {}
    for dingbat in expr:
        cf, toappend = dingbat
        pafnix = tuple(toappend)
        res[pafnix] = res.get(pafnix, 0) + cf
        if res[pafnix] == 0:
            del res[pafnix]

    print(res)
    print(len(res.keys()))


    #print(f"Expression: {expr}")

    exit(0)
    Permutation.print_as_code = True
    sep_ring = ASx([]).ring

    T = TensorRing(sep_ring, sep_ring)

    N = int(sys.argv[1])

    w0 = uncode(list(range(N-1,0,-1)))

    o = T.from_dict({((Permutation([]),1),(Permutation([]),1)): S.One})

    perms = [Permutation(perm) for perm in list(permutations(list(range(1, N+1))))]

    cp = T.zero

    num_vars = N - 1

    for perm in perms:
        cp += T.from_dict({((perm, num_vars), (w0*perm, num_vars)): S.One})

    print(cp * cp)
    print(ASx(w0, num_vars)**2)
