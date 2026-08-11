from schubmult import *
from sympy import prod
import math

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


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = [pp for pp in Permutation.all_permutations(n) if not pp.has_pattern([3,1,2]) and pp.inv >= n]
    for perm in perms:
        domsize = dominant_interval_size(perm.mul_dominant())
        actsize = actual_interval_size(perm)
        L = len([c for c in (~perm).trimcode if c != 0])
        diff = n**L#(perm.mul_dominant().inv - perm.inv + 1)
        assert (domsize/actsize) <= diff, f"{perm} {domsize=} {actsize=} {diff=}"