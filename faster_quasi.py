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
from schubmult.abc import *

def antipode(val):
    res = val.ring.zero
    for k, v in val.items():
        return None

def reverse_elem(val):
    res = val.ring.zero
    for k, v in val.items():
        k0 = tuple(reversed(k))
        res += (S.NegativeOne ** sum(k))*v * val.ring(k0)
    return res

def schub_elem(perm, numvars):
    res = FA([]).ring.zero
    expr = Sx(perm*~uncode(list(range(perm.inv + numvars, perm.inv, -1)))).in_SEM_basis().expand()
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

def double_element(perm, gen_y, max_deg, n):
    SR = SingleSchubertRing(gen_y)
    elem = ASx([]).ring.zero
    stack = [perm]
    while len(stack) > 0:
        new_perm = stack.pop()
        if new_perm.inv > max_deg:
            continue
        elem += ASx(new_perm) * SR(new_perm *(~perm)).as_polynomial()
        l_new_perm = ~new_perm
        for i in range(n):
            if l_new_perm[i] < l_new_perm[i+1]:
                to_add_perm = ~(l_new_perm.swap(i, i+1))
                stack.append(to_add_perm)
    return elem

def noncomsymhom(elem):
    ret = elem.ring.zero
    for bongo, v in elem.items():
        new_tup = tuple([b for b in bongo if b !=0])
        ret += v * elem.ring(new_tup)
    return ret

def compute_vpathdicts_order(th, vmu):
    #th may not be in order
    vpathdicts = [[{},th[i]] for index in range(len(th))]
    vpathdicts[-1][0][vmu] = None
    thL = len(th)

    top = code(~Permutation(uncode(list(sorted(th)))))
    for i in range(thL - 1, -1, -1):
        top2 = code(~Permutation(uncode(top)))
        while top2[-1] == 0:
            top2.pop()
        top2.pop()
        top = code(~Permutation(uncode(top2)))
        # print(f"{top=}")
        monoperm = Permutation(uncode(top))
        # print(f"{monoperm=}")
        k = i + 1
        for last_perm in vpathdicts[i]:
            newperms = kdown_perms(last_perm, monoperm, th[i], k)
            vpathdicts[i][last_perm] = newperms
            if i > 0:
                for trip in newperms:
                    vpathdicts[i - 1][trip[0]] = None
    vpathdicts2 = [{} for i in range(len(th))]
    for i in range(len(th)):
        for key, valueset in vpathdicts[i].items():
            for value in valueset:
                key2 = value[0]
                if key2 not in vpathdicts2[i]:
                    vpathdicts2[i][key2] = set()
                v2 = value[2]
                vpathdicts2[i][key2].add((key, value[1], v2))
    return vpathdicts2


def shuff_hom(key, index):
    frem = ASx(*key).free_element()
    new_free = frem.ring.zero
    for k, v in frem.items():
        new_free += v * frem.ring((*k[:index],0,*k[index:]))
    return new_free.schub_expand()

def monom_sym_rec(partition, numvars, genset):
    if numvars == 0:
        return S.One
    if numvars < 0:
        return S.Zero
    if len(partition) < numvars:
        partition = [*partition, *([0]*(numvars-len(partition)))]
    pm1 = -1
    res = S.Zero
    for i, p in enumerate(partition):
        if pm1 != p:
            pm1 = p
            res += (genset[1] ** p)*monom_sym_rec(partition[:i]+partition[i+1:], numvars-1, genset[1:])
    return res

def monom_sym(tup, numvars):
        from schubmult.abc import e, x, h
        from schubmult.symbolic import Mul, Pow, prod, sympy_Mul
        from schubmult.symmetric_polynomials import degree
        from schubmult.utils.perm_utils import p_trans

        # mu = list(range(numvars, 0, -1))
        # if len(mu) < len(tup):
        #toopunk = p_trans(tup)
        #tup = toopunk
        #mu = [*([numvars] * (len(tup)))]

        flat_part = [*tup]
        painted_bagel = Sx([]).ring.zero
        # flippant_partition = sorted([mu[i] - flat_part[i] for i in range(len(flat_part))], reverse=True)
        flippant_partition = sorted([flat_part[i] for i in range(len(flat_part))], reverse=True)
        # we're multiplying by monomial fat donkey
        # just the inverse kost
        # first command that mu flippant good
        # ptrans e, just schubert

        # penguin = list(sorted(flat_part, reverse=True))
        # #print(f"{penguin=}")
        # flippant_partition = [mu[i] - penguin[i] for i in range(len(flat_part))]
        print(f"{flippant_partition=}")
        lambda_prime = flippant_partition  # p_trans?
        print(f"{lambda_prime=}")
        # lambda_prime.reverse()
        # lambda_prime = p_trans(lambda_prime)
        # #lambda_prime = p_trans(penguin)
        # #lambda_prime = penguin
        total_len = len(tup) - numvars + 1
        perm_construct = ElementaryBasis.grassman_permutation(lambda_prime, numvars+5)
        # results = Sx(perm_construct).in_CEM_basis().expand()

        # e coeffs make Schubert
        def elem_func(p, k, *args):
            from schubmult.rings.free_algebra import FA

            return FA(p, k)

        # get eleme subs Schubert
        froiple = Sx(perm_construct).cem_rep(elem_func)
        print(froiple)
        # print(f"{results=}")
        # # we want a way to convert this to a partition
        # fartbonkkafoof = Sx([]).zero
        park_ranger = Sx([]).ring.zero
        for tupbonk, coeff in froiple.items():
            partition = sorted(tupbonk[2::2], reverse=True)
            #fipple = [0]*numvars
            # for pickles in partition
            # tootbucket = Sx([])
            # tootbucket *= prod([partition])
            print(f"{partition=}")
            park_ranger += coeff * prod([h(a,numvars,[x[0]] + x[1:]) for a in partition])
            # painted_bale
            #         partition1 = mu[:total_len]
            #         for i in range(len(partition)):
            #             partition1[i] = partition1[i] - partition[i]
            #         partition = partition1
            #         partition.reverse()
            painted_bagel[ElementaryBasis.grassman_permutation(partition,numvars)] = coeff
        # print(f"{painted_bagel=} bucketfat")
        print(f"{painted_bagel=}")
        print(f"{painted_bagel.expand()=}")
        # painted_bagel *= prod([x[i + 1] ** (mu[i] - boink_part[i - len(flat_part)]) for i in range(len(flat_part), len(mu))])
        # w0 = ~uncode(mu)
        # monom = {}
        # for k, v in painted_bagel.items():
        #     if (k * w0).inv != w0.inv - k.inv:
        #         raise Exception
        #     monom[(k * w0, numvars)] = v
        # return dict(monom)
        print(f"{park_ranger=}")
        return painted_bagel

def lyndon_words():
    n = 4
    S = [1, 2, 3]
    k = len(S)
    S.sort()

    # To store the indices
    # of the characters
    w = [-1]
    bubbles = []
    # Loop till w is not empty
    while w:

        # Incrementing the last character
        w[-1] += 1
        m = len(w)
        if m == n:
            #print(''.join(S[i] for i in w))
            bubbles.append([S[i] for i in w])
    
        # Repeating w to get a
        # n-length string
        while len(w) < n:
            w.append(w[-m])
    
        # Removing the last character
        # as long it is equal to
        # the largest character in S
        while w and w[-1] == k - 1:
            w.pop()
    return bubbles

def lyndon_factorization(full_word):
    m = 1
    k = 0
    factorization = [[]]
    #working_word = [full_word[0]]
    word = [*full_word]
    while len(word) > 0:
        factorization[-1] += [word[0]]
        while k < len(word) and m < len(word):
            if word[k] == word[m]:
                factorization[-1] += [word[m]]
                k += 1
                m += 1
            elif word[k] < word[m]:
                m += 1            
                k = 0
            else:
                factorization[-1] += [*word[:m-k]]
                word = word[m-k:]
                m = 1
                k = 0
        factorization[-1] += [*word[:m-k]]
        factorization += [[]]
        word = word[m-k:]
        m = 1
        k = 0
    factorization[-1] = word
    return factorization

def duval_factorization(s: str) -> list[str]:
    """
    Factorizes a string into its unique sequence of Lyndon words using Duval's algorithm.

    Args:
        s: The input string.

    Returns:
        A list of strings representing the Lyndon factorization of the input string.
    """
    n = len(s)
    i = 0
    factorization = []

    while i < n:
        j = i + 1
        k = i
        while j < n and s[k] <= s[j]:
            if s[k] < s[j]:
                k = i
            else:
                k += 1
            j += 1

        while i <= k:
            factorization.append(s[i:j - k])
            i += (j - k)
    return factorization

def chen_fox_lyndon_factorization(s):
    n = len(s)
    i = 0
    factorization = []
    while i < n:
        j, k = i + 1, i
        while j < n and s[k] <= s[j]:
            if s[k] < s[j]:
                k = i
            else:
                k += 1
            j += 1
        while i <= k:
            factorization.append(s[i:i + j - k])
            i += j - k
    #assert ''.join(factorization) == s
    return factorization

def factor_lyndon(lyndon_word):
    if len(lyndon_word) <= 1:
        return lyndon_word
    for pants in range(1, len(lyndon_word)):
        word1 = lyndon_word[:pants]
        word2 = lyndon_word[pants:]
        if len(chen_fox_lyndon_factorization(word1)) == 1 and len(chen_fox_lyndon_factorization(word2)) == 1:
            return FA(*word1) * FA(*word2) - FA(*word2)*FA(*word1)

def bacon(perm):
    cd = [*perm.code]
    while len(cd) > 0 and cd[-1] == 0:
        cd = cd[:-1]
    while len(cd) > 0 and cd[0] == 0:
        cd = cd[1:]
    if 0 in cd:
        return None
    return uncode(cd)

def trunc(af):
    zob = af.ring.zero
    print(af)
    for (k, n), v in af.items():
        key = bacon(k)
        if key:
            zob += v * ASx(key)
    return zob

def insert_in_middle(a, b, spot):
    a_word = a.change_basis(WordBasis)
    b_word = b.change_basis(WordBasis)
    woven_fat = FA().ring.zero
    for k1, v1 in a_word.items():
        for k2, v2 in b_word.items():
            woven_fat += FA(*k1[:spot], *k2, *k1[spot:]) * v1 * v2
    return woven_fat.change_basis(a.ring._basis)

sc_dict = {}

@cache
def precompute_schub_quasisym(n):
    from sage.all import Compositions, QuasiSymmetricFunctions, ZZ    
    
    QSym = QuasiSymmetricFunctions(ZZ)
    M = QSym.M()
    dct = {(): QSym.one()}
    for i in range(1, n + 1):
        for c in Compositions(i):
            cc = tuple([int(a) for a in c])
            j = FA(*cc).change_basis(JBasis)
            for k, v in j.items():
                dct[k] = dct.get(k, QSym.zero()) + v * M[*cc]
    return dct

@cache
def schub_quasisym(comp):
    dct = precompute_schub_quasisym(int(sum([c for c in comp])))
    return dct[tuple([int(cc) for cc in comp])]
    # from sage.all import QuasiSymmetricFunctions, ZZ, Compositions
    # from schubmult.rings.variables import genset_dict_from_expr
    # QSym = QuasiSymmetricFunctions(ZZ)
    # M = QSym.M()
    # perm = uncode(([0]*5)+[int(cc) for cc in comp])
    # # if numvars < varlen(perm= perm):
    # #     return S.Zero
    # # if numvars == varlen(perm):
    # #     if 0 in perm.code[:numvars]:
    # #         return S.Zero
    # # cd = perm.code[:varlen(perm)]
    # # cd = [*([0] * (perm.inv-varlen(perm))), *cd]
    # # polyfuck = Sx(uncode(cd)).expand()
    # # dct = genset_dict_from_expr(polyfuck, Sx([]).ring.genset)
    # ret = QSym.zero()
    # compint = tuple([int(cc) for cc in comp])
    # for comp2 in Compositions(perm.inv):
    #     comp2int = [int(cc) for cc in comp2]
    #     stib = FA(*comp2int).change_basis(JBasis)
    #     ret += stib.get(compint, QSym.zero()) * M[*comp2int]
    # for k, v in dct.items():
    #     try:
    #         viv = k.index(0)
    #     except ValueError:
    #         viv = -1
    #     if viv != -1:
    #         if k[viv:].count(0) != len(k[viv:]):
    #             continue
    #         ret += v * M[*k[:viv]]
    #     else:
    #         ret += v * M[*k]
    
    #print(f"{comp} {ret=}")
    return ret
    
def unsage(stinkfuck):   
    val = sympify(str(stinkfuck))
    subs_dict = {Symbol(f"x{i}"): Sx([]).ring.genset[i+1] for i in range(50)}
    val = efficient_subs(val, subs_dict)
    return val
    
              
if __name__ == "__main__":
    #bubbles = lyndon_words()

    Permutation.print_as_code = True

    #lyndon words positive
    

    n0 = 6
    n1 = 5  
    numvars = 3
    w0 = uncode(list(range(n1-1,0,-1)))
    perms = Permutation.all_permutations(n0)
    #perms1 = Permutation.all_permutations(n1)

    from schubmult.rings import *
    from schubmult.rings.free_algebra_basis import *
    from sage.all import Compositions
    dct = {}
    from sage.all import QuasiSymmetricFunctions, ZZ, Compositions, PolynomialRing
    tt = PolynomialRing(ZZ, 't')
    dct = {}
    QSym = QuasiSymmetricFunctions(tt)
    M = QSym.M()
    
    alphagod = (2,1,2)
    ret = QSym.zero()
    stack = [[alphagod,[]]]
    while len(stack) > 0:
        this_alpha = stack.pop()
        if len(this_alpha[0]) == 0:
            ret += M[*this_alpha[1]]
        else:
            # stinkbag = Sx(uncode(this_alpha[0]))
            # pilfer = stinkbag.pull_out_gen(stinkbag.ring.genset[1])
            asum = sum(this_alpha[0])
            # stinkset = set()
            # for k, _ in pilfer.items():
            #     fingbat = trimcode(k)
            #     fsum = sum(fingbat)
            #     new_alpha = [*fingbat]
            #     new_data = [*this_alpha[1], asum - fsum]
            #     stack.append([new_alpha, new_data])

            stinkbag2 = Sx(uncode([0] + [*this_alpha[0]]))
            pilfer2 = stinkbag2.pull_out_gen(stinkbag2.ring.genset[1])
            for k, _ in pilfer2.items():
                fingbat = trimcode(k)
                if 0 in fingbat[1:]:
                    continue
                fsum = sum(fingbat)
                if asum != fsum and (len(fingbat)!=0 and fingbat[0] != 0):
                    new_alpha = [*fingbat]
                    new_data = [*this_alpha[1], asum - fsum]
                    stack.append([new_alpha, new_data])
                elif len(fingbat) == 0 or fingbat[0] == 0:
                    new_alpha = [*fingbat[1:]]
                    new_data = [*this_alpha[1], asum - fsum]
                    stack.append([new_alpha, new_data])
        # Process this_alpha
        # If you need to add new elements to the stack, do so here
        # For example:
        # stack.append([new_alpha, new_data])
    # for i in range(1, n0):
    #     for c in Compositions(i):
    #         for pw in range(1):
    #             cint0 = tuple([int(a) for a in c])
    #             cint = tuple([0]*pw + [int(a) for a in c])
    #             #print(f"{cint} {schub_quasisym(cint)}")
    #             bobjones = FA(*list(cint)).change_basis(JTBasis)
    #             for k2, v2 in bobjones.items():
    #                 dct[k2[0]] = dct.get(k2[0], QSym.zero()) + int(v2) * (tt.gens()[0]** (k2[1]-pw)) * M[*cint0]
                #print(f"{cint}: {bobjones}")
    # QS = QSym.QS()
    print(f"J({alphagod}) = {ret}")
    for k, v in dict(ret).items():
        kk = tuple(k)
        fa_elem = FA(*kk).change_basis(JTBasis)
        for k2, v2 in fa_elem.items():
            if k2[0] == alphagod:
                print(f"{kk}: {v2}")
    sys.exit(0)
    def comp(p, perm):
        if perm.inv > 0 and max(perm.descents()) >= p:
            return None
        mx = max([perm.code[i] + i if i<len(perm.code) else i for i in range(p)]) + 1
        return uncode([mx - 1 - i  - perm.code[i] if i<len(perm.code) else mx - 1 - i for i in range(p)])

    EE = FreeAlgebra(ElementaryBasis)
    bigness = EE(tuple(list(range(1,n0))),n0-1)
    a = 4
    b = 4
    mat = []
    for i in range(a):
        
        for j in range(b):
            pol = DSx([]) * x[1]**i * x[2] ** j
            mat += [pol]
    keylist = set()
    for dingbat in mat:
        keylist = keylist.union(set((k for k in dingbat.keys() if dingbat[k] != S.Zero)))
    rowbucket = list(sorted(keylist))
    from sympy import maple_code
    for matrow in mat:
        print(f"[{', '.join([maple_code(sympify_sympy(matrow.get(row,S.Zero))) for row in rowbucket])}]")
        

        # for p in range(1, n0):
        #     res = comp(p, perm)
        #     if res:
        #         print(f"{trimcode(perm),p} {trimcode(res)}")
        # word = ASx(perm).change_basis(WordBasis).change_basis(JTBasis)
        # print(word,file=sys.stderr)
        # print(f"{trimcode(perm)}: {word}")
        # loin = max(perm.descents()) + 1
        # for compkey in sorted(dct.keys()):
        #     if perm.inv != sum(compkey):
        #         continue
        # # print(f"{compkey}: {dct[compkey]}")
        # # for j in range(2, sum(compkey)+1):
        # #     print(f"{j}: {Sx(unsage(dct[compkey].expand(j)).subs(Symbol('t'),S.One))}")
        # # print("")
        #     pancreas = {tuple([int(a) for a in k]): v for k, v in dict(dct[compkey]).items()}
            
            
        #     smep = 0
        #     for tup, v in word.items():
        #         if tup[0] in pancreas:# and (loin - len(tup) == 0 or str(pancreas[tup]).find(f"**{len(tup) - loin}") != -1):
        #             smep += v * pancreas[tup[0]] * (tt.gens()[0]**tup[1])
        #     smep = unsage(smep)
        #     #smep = sympify(smep).subs(Symbol('t'), S.One)
        #     if smep != 0:
        #         print(f"{trimcode(perm),compkey}: {smep}")
    #     pigpoly = unsage(dct[compkey].expand(len(compkey)+1))
    #     print(f"{Sx(pigpoly)}")
        #print(f"SS{compkey}: {QS(dct[compkey])}")
    # for perm in perms:
    #     # if 0 in trimcode(perm) or perm.inv == 0:
    #     #     continue
    #     # for i in range(1, n0):
    #         # if 0 in trimcode(perm2):
    #         #     continue
    #     # if 0 in perm.code[:varlen(perm)]:
    #     #     continue
    #     #jstick = J(*perm.code[:varlen(perm)]).change_basis(WordBasis)
    #     #print(f"{jstick=}")
    #     bobjones = ASx(perm).change_basis(WordBasis).change_basis(JTBasis)
    #     # for i in range(1, n0):
    #     #     for c in Compositions(i):
    #     #         cint = tuple([int(a) for a in c])
    #     #         #print(f"{cint} {schub_quasisym(cint)}")
    #     #         bobjones = FA(*cint).change_basis(SchubertBasis).change_basis(WordBasis).change_basis(JTBasis)
    #     print(f"{trimcode(perm)}: {bobjones}")
        #print(f"{trimcode(perm)}: {bobjones}")
        # jstick = ASx(perm).change_basis(WordBasis).kill_zero()
        # for perm2 in perms:
        #     if perm.inv != perm2.inv or 0 in perm2.code[:varlen(perm2)] or varlen(perm2) > varlen(perm):
        #         continue
        #     # gong = J(*perm.code[:varlen(perm)]).poly_inner_product(Sx(perm2).expand(),Sx([]).ring.genset)
        #     # gong = J(*perm.code[:varlen(perm)]).poly_inner_product(schub_quasisym(perm2, varlen(perm)+3), Sx([]).ring.genset, None)
        #     # gong = ASx(perm).poly_inner_product(unsage(schub_quasisym(perm2).expand(varlen(perm))), Sx([]).ring.genset, varlen(perm))
        #     gong = 0
        #     sq = schub_quasisym(tuple(perm2.code[:varlen(perm2)]))
        #     # print(f"{sq=}")
            
            
        #     for k, v in dict(sq).items():
        #         k = tuple(k)
        #         # print(f"{k=}")
        #         gong += int(v) * int(jstick.get(k, 0))

        #     # gong = J(*perm.code[:varlen(perm)]).poly_inner_product(Sx(perm2).expand(), Sx([]).ring.genset, varlen(perm2))
        #     if gong != 0:
        #         print(f"{trimcode(perm),trimcode(perm2)}: {gong}")
        # vl = varlen(perm)
        # fat_pig = Sx(perm)
        # for i in range(1, vl + 1):
        #     dicky = fat_pig.coporudct(*list(range(1,i+1)))
        #     for (k1, k2), v in dicky.items():
        #         if 0 not in k1.code[:varlen(k1)] and 0 not in k2.code[:varlen(k2)]:
        #             dct[(k1,k2)] = dct.get((k1,k2), ASx([]).ring.zero) + v * fat_pig.ring(k1,k2)
    # for perm in perms:
    #     spiff = perm.code[:varlen(perm)]
    #     piffnol = FA(*spiff).change_basis(SchubertBasis)
    #     boing = next(iter(piffnol.kill_zero().keys()))
    #     dct[boing] = dct.get(boing, ASx([]).ring.zero) + piffnol
    #     # fif = ASx(perm).change_basis(WordBasis)
    #     # fif.kill_zero()
    #     # for k, v in fif.items():
    #     #     dct[k] = dct.get(k, ASx([]).ring.zero) + v * ())
    
    # for dingbet, v in dct.items():
    #     print(f"{dingbet}: {v}")
    # for perm1 in perms:
    #     # print(f"{trimcode(perm1)}: {ASx(perm1).kill_zero(True,Symbol('t')).change_basis(WordBasis).change_basis(JBasis)}")
    #     print(f"{trimcode(perm1)}: {ASx(perm1).kill_zero(True,Symbol('t')).change_basis(WordBasis).change_basis(ZBasis)}")
    # for n0 in range(1,6):
    # #     # for n1 in range(1,6):
    #     for c in Compositions(n0):
    # #         zoingle = ASx(uncode([int(cc) for cc in c])).kill_zero(True, Symbol('t')).change_basis(WordBasis)
    # #         #zingle = J(*[int(cc) for cc in c])
    #         zatbong = FA(*[int(cc) for cc in c]).change_basis(JBasis)
    #         print(f"{c}: {zatbong}")
    #         print(f"{c}: {zoingle.change_basis(WordBasis).change_basis(ZBasis)}")
            #print(f"{c}: {zingle.change_basis(WordBasis).change_basis(ZBasis)}")
                # for c2 in Compositions(n1):
                #     piss = Z(*[int(cc)+1 for cc in c])
                #     piss1 = Z(*[int(cc)+1 for cc in c2])
                #     pisspiss = piss.change_basis(WordBasis) * piss1.change_basis(WordBasis)
                #     piss2 = FA.from_dict({tuple([a-1 for a in k]): v for k,v in pisspiss.items()})
                #     jj = piss2.kill_zero(True)
                #     pinto = jj.change_basis(JBasis)
                #     pants = J(*[int(cc) for cc in c])
                #     pants2 = J(*[int(cc) for cc in c2])
                #     fisbut = (pinto - (pants*pants2))
                #     assert len(fisbut.items()) == 0 or all([v == S.Zero for v in (pinto - (pants * pants2)).values()]), print(f"{pants*pants2=} {pinto=} {fisbut=}")
                #     print(f"{c}: Yay")
                    
        #print(f"{c}: {piss2.change_basis(SchubertBasis)}")
        # if perm1.inv == 0:
        #     continue
        # for perm2 in perms:
        #     print(insert_in_middle(ASx(perm1), ASx(perm2), 1))
        # for pin in range(varlen(perm1), n0):
        #     spingo = ASx(perm1,pin).change_basis(WordBasis)
        #     newbucket = FA().ring.zero
        #     for tup, v in spingo.items():
        #         newbucket += v * FA(*[a for a in tup if a != 0])
        #     print(f"{perm1.code,pin}: {newbucket.change_basis(JBasis)}")
        # for perm2 in perms:
        #     # for a1 in range(varlen(perm1),n0):
        #     #     for a2 in range(varlen(perm2),n0):
        #     perm1 = bacon(perm1)
        #     perm2 = bacon(perm2)
        #     if perm1 is None or perm2 is None:
        #         continue
        #     a1 = varlen(perm1)
        #     a2 = varlen(perm2)
        #     spoing = ASx([]).ring.from_dict((ASx(perm1) @ ASx(perm2)))
        #     try:
        #         assert all(a >= 0 for a in spoing.values()) or len(perm2.descents()) > 1
        #     except AssertionError:
        #         print(f"Failed {perm1.code, a1} {perm2.code, a2}")
        #         print(f"{spoing=}")
        #         sys.exit(1)
        #     print(f"Success {perm1.code, a1} {perm2.code, a2}")
        

    # can just add 1
    sys.exit(0)
    SX = PolynomialAlgebra(basis=SchubertPolyBasis(numvars, Sx([]).ring))

    # for perm in perms:
    #     if varlen(perm) > numvars:
    #         continue
    #     if 0 in perm.code[:numvars]:
    #         continue
    #     wordbungle = ASx(perm, numvars).change_basis(WordBasis)
    #     if len(perm.descents()) == 1:
    #         print(f"{perm.code=}")
    #         for element, coeff in wordbungle.items():
    #             if 0 not in element:
    #                 print(f"{chen_fox_lyndon_factorization(element)}: {coeff}")

    print("PINGUS")
    for perm in perms:
        Permutation.print_as_code = False
        pednool = ASx(perm).change_basis(WordBasis)
        # born = FreeAlgebraBasis.change_tensor_basis(pednool, WordBasis, WordBasis)
        reborn = pednool.ring.zero
        for k, v in pednool.items():
            if 0 not in k:
                reborn += v * reborn.ring(*k)
        reborn = reborn.coproduct()
        fatbucket = reborn.ring.zero
        for (k1, k2), v in reborn.items():
            fatbucket += v * reborn.ring((tuple([zob for zob in k1 if zob != 0]), tuple([zob for zob in k2 if zob != 0])))
        pigbasin = FreeAlgebraBasis.change_tensor_basis(fatbucket, SchubertBasis, SchubertBasis)
        repigbasin = pigbasin.ring.zero
        for (k1, k2), v in pigbasin.items():
            if 0 not in k1[0].code[:k1[1]] and 0 not in k2[0].code[:k2[1]]:
                repigbasin += v * repigbasin.ring((k1, k2))
        print(f"{perm.code}: {repigbasin}")
        # if varlen(perm) > numvars:
        #     continue
        # if 0 in perm.code[:numvars]:
        #     continue
        # wordbungle = ASx(perm, numvars).change_basis(WordBasis)
        # if len(perm.descents()) == 1:
        #     if len(perm.descents()) > 1:
        #         print("MOOGLEFONK")
        #     else:
        #         print("GRASSGRASS")
        #     print(f"{perm.code=}")
        #     for element, coeff in wordbungle.items():
        #         if 0 not in element:
        #             print("BANGTHEDONKEY")
        #             fat_pig = chen_fox_lyndon_factorization(element)
        #             for pingus in fat_pig:
        #                 print(f"{factor_lyndon(pingus)}: {coeff}")
        #         else:
        #             print("FATSNICKERNO")
        # if varlen(perm) > numvars:
        #     continue
        # print(SX(perm).change_basis(SX._basis.monomial_basis).change_basis(SX._basis))
        # #print(PA.from_expr(x[5]**4).change_basis(SX._basis))
        # pants = SX(perm)
        # print(f"{pants=}")
        # print(f"{pants.change_basis(ElemSymPolyBasis(SX._basis.numvars)).coproduct()=}")
        #raise Exception("boafns")
    
    #perms2 = Permutation.all_permutations(6)
    
    sys.exit(0)
    from schubmult.abc import *
    from schubmult.perm_lib import *
    from schubmult.rings import WordBasis, SchubertBasis
    #from schubmult.schub_lib import kdown_perms
    Permutation.print_as_code = True
    
    #@cache

    good = 0
    bad = 0
    # basis class
    # SSR = FreeAlgebra(basis=SchubertSchurBasis)
    # ring = FreeAlgebra(basis=SchubertSchurBasis)
    def elem_func(p, k,*args):
        return FA(p,k)
    
    #for perm in perms:
        # print(Sx(perm).cem_rep(elem_func))
    Elem = FreeAlgebra(basis=ElementaryBasis)
    #print(Elem((3,3,3,0,0,0),4).change_basis(SchubertBasis))
    pants = Elem((1,0,2,2,2,3,4),4).change_basis(SchubertBasis)
    print(f"{pants=}")
    print(f"{pants.change_basis(ElementaryBasis)=}")
    #print(monom_sym_rec((3,1,1),4,x))
    sys.exit(0)
    # for perm0 in perms:
    #     spangledong0 = SSR((0,)*(n0), perm0).change_basis(SchubertBasis)
    # perm_set = set()
    # for k in range(3,5):
    #     for perm1 in perms1:
    #         #perm0 = Permutation([])
    #         numvars = varlen(perm1)
    #         # if numvars < k:
    #         #     continue
    #         # if varlen(perm1) > 4:
    #         #     continue
    #         # if numvars < 3:
    #         #     continue
    #         spangledong1 = ASx(perm1,numvars)#SSR((1,)*(n0), perm1).change_basis(SchubertBasis)
    #         spangledong2 = spangledong1.change_basis(SeparatedDescentsBasis(k))
    #         if len(spangledong2.items()) == 1:
    #             perm_set.add(perm1*w0)
    #             print(f"{perm1.code}, {numvars}: {spangledong2}, {(perm1*w0).code=}")
            
        # spangledong = SSR((*((0,)*1),*((1,)*(n0))), perm0).change_basis(SchubertBasis)
        # res = (ASx(perm1)*spangledong).change_basis(SchubertSchurBasis)
        # print(f"{perm.code}, {perm0.code}, fattybagel: {res}")
            
        # if perm.inv == 0 or varlen(perm) >  numvars:
        #     continue
        # print(f"{perm0=}")
        # elem1 = ASx(perm0,3)
        # cprd = elem1.coproduct()
        # result = FreeAlgebra.change_tensor_basis(cprd, SchubertSchurBasis, SchubertSchurBasis)
        # #print(f"{prd=}")
        # print(f"{cprd=}")
        # print(f"{result=}")
            # print(f"{elem=}")
            # word_elem = elem.change_basis(SchubertBasis)
            # print(f"{word_elem=}")
            # elem_elem = (word_elem*word_elem).change_basis(SchubertSchurBasis)
            # print(elem_elem)
        
        # print(f"{perm.code}: {Sx([]).ring.in_schubert_schur_basis(perm, numvars)}")
        # print(f"{perm.code}: {Sx([]).ring.in_descending_schur_basis(perm, numvars)}")
                              #(perm, numvars)}")
        # spots = [2,3]
        # paco = Sx(perm).coproduct(*spots)
        # # off the end
        # for (p2, p1), v in paco.items():
        #     # 4 5
        #     pringles = Permutation.sorting_perm([-500]+list(p1[1:]))
        #     pringles2 = Permutation.sorting_perm(list(range(-500,-500+3))+list(p1[1:]))
        #     print(f"{pringles=}")
        #     print(f"{perm=}")
        #     print(f"{p1=}")
        #     new_perm = perm * (~pringles2)
        #     new_p1 = p1*(~pringles)
        #     print(f"{new_perm=}")
        #     print(f"{new_p1=}")
            
        #     assert perm.inv - new_perm.inv == p1.inv - new_p1.inv
        #     # numvars 1
        #     # numvars n0 - 2
        #     print(f"{p1.code=}")
        #     print(f"{p2.code=}")
        #     bango = ASx(new_p1,1)*ASx(p2,n0-2)
        #     print(f"{bango=}")
        #     assert bango.get((new_perm,n0-1), S.Zero) == v
        # #print(f"{perm.code}: {shuff_hom((perm, numvars), 2)}")
        # print(f"Success {perm.code}")
    # print(f"{len(perm_set)=}")
    # for perm in perm_set:
    #     print(perm.code)
    sys.exit(0)

    MAX_DEGREE = 6
    #numvars = 4
    bacon = set()
    flonf = set()
    result = []
    # from schubmult.rings.split_free_algebra_module import SplitFreeAlgebraModule
    #mod = SplitFreeAlgebraModule(FreeAlgebra, 1)
    #mod_elem=mod.from_dict({tuple(): S.One})
    #mod_elem=FA([]).ring.schub_elem(uncode([1,2]),numvars)*mod_elem
    for perm in perms:
        
        # if varlen(perm) > splitter:
        #     #print(f"{perm.code}: {}")
        if len(perm.descents()) > 0:
            continue
        minvars = varlen(perm)
        for numvars in range(minvars,n0):
            elem0 = ASx(perm,numvars).free_element()
            for splitter in range(1,numvars):
                elem = FA([]).ring.tensor_schub_expand(elem0.split(splitter))
                #if not any(v < 0 for v in elem.values()):
                print(f"({perm.code},{numvars}): {elem}")
                impute = sum([v*ASx(*a)*ASx(*b) for (a,b), v in elem.items()])
                assert impute == ASx(perm, numvars)
        #for numvars in range(varlen(perm), n0):
        #print(f"{perm}: {coprod(ASx(perm))}")    
        # if varlen(perm) > numvars:
        #     continue
        # print(f"{perm.code}: {reverse_elem(schub_elem(perm, numvars)).schub_expand()}")
        #numvars = varlen(perm)
        # 3 step solves schubert schur!!!!
        #print(f"{perm}: {FA([]).ring.from_dict(FA([]).ring.schub_elem(perm,numvars)*mod_elem).schub_expand()}")
        # # floff = FA([]).ring.schub_elem(perm, numvars).nsymexpand()
        # # print(f"{perm.code}: {floff}")
        # # floff = FA([]).ring.schub_elem(perm, numvars).remove_zeros(inserter = S.Zero).schub_expand()
        # # print(f"{perm.code} bacon: {floff}")
        # print(f"{perm.code}: {FA([]).ring.tensor_schub_expand(FA([]).ring.schub_elem(perm, numvars).coproduct())}")
    #     if tuple(floff.items()) not in bacon and any(v != S.Zero for v in floff.values()):
    #     # new_floff = ASx([]).ring.zero
    #     # old_new_floff = floff
    #     # while any(floff[k] - new_floff.get(k, S.Zero) != S.Zero for k in floff):
    #     #     # print(f"    {floff=}")
    #     #     # print(f"{new_floff=}")
    #     #     if new_floff != ASx([]).ring.zero:
    #     #         old_new_floff = new_floff
    #     #         floff = new_floff
    #     #     floff = floff.free_element().remove_zeros().schub_expand()
    #     #     new_floff = ASx([]).ring.zero
    #     #     for k, v in floff.items():
    #     #         b0 = list(tuple([*k[0].code]))
    #     #         while len(b0) > 0 and b0[-1] == 0:
    #     #             b0.pop()

    #     #         while len(b0) > 0 and b0[0] == 0:
    #     #             b0.pop(0)
    #     #         new_floff += v * ASx(uncode(b0), len(b0))
    #     #     if old_new_floff == new_floff:
    #     #         break
    #         flonf.add(perm)
    #         bacon.add(tuple(floff.items()))
    #         result += [(perm, floff)]
    
    # for perm, bonk in result:
    #     if perm not in flonf:
    #         continue
    #     print(f"{trimcode(perm)}: {bonk}")
        # print(f"{perm.code}, {numvars}: {new_floff.free_element()}")
        # borb = Sx(perm).in_SEM_basis()
        # fipple = sympify_sympy(borb)
        # fipple = fipple.replace(lambda x: is_of_func_type(x, FactorialElemSym), lambda xx: h(xx.degree, max([i + c + 1 for i, c in enumerate(perm.code)]) - xx.numvars, x[1:]))
        # print(f"{perm.code}: {fipple*Sx([])}")
        
        #print(f"{perm.code}, {numvars}: {FA([]).ring.schub_elem(perm, numvars).nsymexpand()}")

        # fa = FA([]).ring.schub_elem(perm, numvars)
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
