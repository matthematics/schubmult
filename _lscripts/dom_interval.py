from schubmult import *
import math
from sympy import prod
from functools import cache

def dominant_interval_size(dom_perm):
    if dom_perm.inv == 0:
        return 1
    # if not dom_perm.is_dominant:
    #     raise ValueError(f"{dom_perm.trimcode} not dominant")
    dpn = ~dom_perm
    dpc = dpn.code
    divv = prod([dpn[i] - dpc[i] if i < len(dpc) else dpn[i] for i in range(len(dpn))])
    return math.factorial(len(dom_perm))//divv

def bruh_size(perm):
    from schubmult.utils.perm_utils import has_bruhat_descent
    stack = [perm]
    seen = set()
    while len(stack) > 0:
        new_perm = stack.pop()
        if new_perm in seen:
            continue
        seen.add(new_perm)
        for i in range(len(new_perm) - 1):
            for j in range(i + 1, len(new_perm)):
                if has_bruhat_descent(new_perm, i, j):
                    stack.append(new_perm.swap(i, j))
    return len(seen)


def potato_sum(n):
    perms = [pp for pp in Permutation.all_permutations(n) if not pp.has_pattern([3,1,2])]
    bs = 0
    for perm in perms:
        md = perm.mul_dominant()
        da_perm = perm * (~md)
        bs = max(bs, bruh_size(da_perm))
    print(f"{bs} pants size vs {len(perms)}")

@cache
def actual_interval_size(perm):
    stack = [perm]
    seen = set()
    while len(stack) > 0:
        new_perm = stack.pop()
        if new_perm in seen:
            continue
        seen.add(new_perm)
        for i in (~new_perm).descents():
            stack.append(~((~new_perm).swap(i, i + 1)))
    return len(seen)

def ass_interval_size(perm):
    #perm = perm
    perms = Permutation.all_permutations(len(perm))
    #return len([perr for perr in perms if len((~perr).code) <= len((~perm).code) and all((~perr).code[i] <= (~perm).code[i] for i in range(len((~perr).code)))])
    return len([perr for perr in perms if perr.inversion_set.issubset(perm.inversion_set)])

def dp_interval_test(n):
    dom_perms = [perm for perm in Permutation.all_permutations(n) if perm.is_dominant]
    for dom_perm in dom_perms:
        assert dominant_interval_size(dom_perm) == actual_interval_size(dom_perm), f"failed for {dom_perm.trimcode},\n{dominant_interval_size(dom_perm)=}\n{actual_interval_size(dom_perm)=}"
        print("Socks")

def dp_interval_test2(n):
    dom_perms = [perm for perm in Permutation.all_permutations(n) if perm.is_dominant]
    for dom_perm in dom_perms:
        assert ass_interval_size(dom_perm) == actual_interval_size(dom_perm), f"failed for {dom_perm.trimcode},\n{ass_interval_size(dom_perm)=}\n{actual_interval_size(dom_perm)=}"
        print("Socks")


def dom_size_diff_test(n):
    perms = Permutation.all_permutations(n)

    for perm in perms:
        cd = [*perm.trimcode]
        dom_perm = perm.minimal_dominant_above()
        diff = 0
        for i in range(len(cd) - 1, 0, -1):
            for j in range(i - 1, -1, -1):
                if cd[i] > cd[j]:
                    cd[i] += 1
                    #diff += 1
        new_perm = uncode(cd)
        diff = new_perm.inv - perm.inv
        assert diff == dom_perm.inv - perm.inv, f"failed {perm.trimcode=} {dom_perm.trimcode=}"
        if new_perm != perm:
            print(f"{new_perm=} {perm=}")

def dom_interval_increasing_test(n):
    perms = Permutation.all_permutations(n)
    for perm in perms:
        sz = dominant_interval_size(perm)
        for d in (~perm).descents():
            perm2 = ~((~perm).swap(d, d + 1))
            sz2 = dominant_interval_size(perm2)
            assert sz2 <= sz, f"failed {perm.trimcode=} {perm2.trimcode=} {sz=} {sz2=}"
        print("Farting stinkbat")

def dom_interval_bound_test(n):
    perms = [popo for popo in Permutation.all_permutations(n) if not popo.has_pattern([3,1,2])]
    pigeon = {}
    w0 = Permutation.w0(n)
    for perm in perms:
        assert (w0 * perm).is_dominant
        size1 = dominant_interval_size(perm)
        size2 = actual_interval_size(perm)
        # size3 = dominant_interval_size(perm * w0)
        # size4 = dominant_interval_size(perm.minimal_dominant_above()*w0)
        # assert size1 <= size2, f"failed {perm.trimcode=} {size1=} > {size2=}"
        # assert size3 >= size4, f"finkas {perm.trimcode=} {size4=} > {size3=}"
        assert size1 == size2
        # dom1 = ~((~perm).minimal_dominant_above())
        # pigeon[dom1] = pigeon.get(dom1, 0) + 1
        # assert all(perm.code[i] <= dom1.code[i] for i in range(min(perm.inv + 1,len(perm)) - 1))
        #diff == dom_perm.inv - perm.inv, f"failed {perm.trimcode=} {dom_perm.trimcode=}"
        # if new_perm != perm:
        #     print(f"{new_perm=} {perm=}")
        print("Stinkbat")
    # for dp, ct in pigeon.items():
    #     assert ct == actual_interval_size(dp)

def bad_ones(v):
    """ <= dom(v) but not <= v """
    dom = v.mul_dominant()
    if v == dom:
        return set()
    stack = [dom]
    seen = set()
    while len(stack) > 0:
        new_perm = stack.pop()
        if new_perm in seen or new_perm.inversion_set.issubset(v.inversion_set):
            continue
        seen.add(new_perm)
        for i in (~new_perm).descents():
            perm2 = ~((~new_perm).swap(i, i + 1))
            stack.append(perm2)
    return seen

def puke(n):
    import itertools
    from schubmult.symbolic import expand
    import math

    perms = Permutation.all_permutations(n)
    mxlgdiff = {}
    for perm1, perm2 in itertools.product(perms, repeat=2):
        if dominant_interval_size(perm2) == 1:
            continue
        md2 = perm2.mul_dominant()
        if perm2 == md2:
            continue
        #spinach
        # spain_pain = uncode([a + b for a,b in itertools.zip_longest(md1.trimcode, md2.trimcode, fillvalue=0)])
        # n1, n2 = dominant_interval_size(spain_pain), 2**(n-1) * math.sqrt(dominant_interval_size(md1) * dominant_interval_size(md2))
        # assert n1 <= n2
        # print(f"{math.log(n2/n1)}")
        purple1 = len({k for k,v in (DSx(perm1)*DSx(perm2, "z")).items() if expand(v) != 0})
        purple2 = len({k for k,v in (DSx(perm1)*DSx(md2, "z")).items() if expand(v) != 0})
        assert purple2 <= purple1 * dominant_interval_size(md2)/dominant_interval_size(perm2), f"failed {perm1.trimcode=} {perm2.trimcode=} {purple1=} {purple2=} {dominant_interval_size(md2)=}"
        logbat = math.log(dominant_interval_size(md2)/dominant_interval_size(perm2))/math.log(n)
        # prt = logbat if logbat > 1 else None
        # # length of code(perm2) is precisely length of code(md2)
        # # code_m(md2) = code_m(perm2)
        # # if code_{m-1}(perm2) >= code_m(perm2), code_{m-1}(md2) = code_{m-1}(perm2)
        # # if code_{i-1}(perm2) <= code_i(perm2), code_{i-1}(md2) = code_i(md2)
        # # if code_{i-1}(perm2) > code_i(perm2), code_{i-1}(md2) = max(code_{i}(md2) + 1, code_{i-1}(perm2))
        # dev_m = 0
        # dev_{i-1} = dev_{i} + c_i - c_{i-1}
        # dev_{i-1} = max(dev_i + c_i - c_{i-1}, 0)
        # max jump at i <= n - i
        # (n - i)^{i - 1}
        # (i - 1)log ( n - i)
        # log(n - i) - (i - 1)/(n - i) = 0
        # # dom is (height)^length
        # # grass >= (height)^(length - 1) 
        # if prt is not None:
        #     print(f"Batbat {prt}")
        assert logbat < 2, f"failed {perm1.trimcode=} {perm2.trimcode=} {purple1=} {purple2=} {dominant_interval_size(md2)=} {dominant_interval_size(perm2)=} {logbat=}"
        print("Pinto")

# def bruhat_interval(w):
#     ret = {w}

#     stack = [w]
#     while len(stack) > 0:
#         w2 = stack.pop()
#         if w2 in ret:
#             continue
#         ret.add(w2)
#         for i in range(len(w2) - 1):
#             for j in range(i + 1,len(w2)):
#                 if has_bruhat_descent(w2, i, j):
#                     stack.append(w2.swap(i,j))
#     return ret

def one32231(n):
    perms = set([(porn,porn*((~porn).minimal_dominant_above())) for porn in Permutation.all_permutations(n) if not porn.has_pattern([3,1,2])])
    print(f"{len(perms)=}")
    print(f"{sorted(set([(bruh_size(pp), actual_interval_size(ppp)) for (ppp,pp) in perms]))}")
    # 312 132 AVOIDANT INTERVAL SIZE
    #for perm in perms:

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    #dom_interval_increasing_test(n)
    puke(n)
    #one32231(n)
    #potato_sum(n)
    #dom_interval_bound_test(n)