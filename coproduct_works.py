from schubmult import *
from schubmult.rings import *
from schubmult.symbolic import S
from itertools import permutations
import sys
from functools import cache

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

def coprod(val, T):
    if val == val.ring.zero:
        return T.zero
    
    cprd_val = T.zero

    while val != val.ring.zero:
        mx = [k[0].code for k in val.keys() if val[k] != S.Zero ]
        mx.sort(reverse=True)
        cd = mx[0]
        

        mx_key = next(iter([k for k in val.keys() if k[0].code == cd]))
        if len(cd) == 0:
            return cprd_val + T.from_dict({((Permutation([]),mx_key[1]),(Permutation([]),mx_key[1])): val[mx_key] * S.One})
        cd = [*cd]
        fv = cd.pop(0)
        while len(cd) > 1 and cd[-1] == 0:
            cd.pop()
        cf = val[mx_key]
        cprd_val += (T.from_dict({((Permutation([]),0),(Permutation([]),0)): cf*S.One}))* single_coprod(fv, 1, T) * coprod(ASx(uncode(cd), mx_key[1] - 1), T)
        val -= cf * ASx(uncode([fv]),1)*ASx(uncode(cd), mx_key[1] - 1)
        print(val)
        
    return cprd_val


def coprod_back(val, T):
    if val == val.ring.zero:
        return T.zero
    
    cprd_val = T.zero

    while val != val.ring.zero:
        mx = [k[0].code for k in val.keys() if val[k] != S.Zero ]
        mx.sort()
        cd = mx[0]
        

        mx_key = next(iter([k for k in val.keys() if k[0].code == cd]))
        cd = [*cd]
        while len(cd) > 0 and cd[-1] == 0:
            cd.pop()
        
        print(f"{cd=}")
        fv = cd.pop()
        
        cf = val[mx_key]
        if len(cd) == 0:
            cprd_val += cf * single_coprod(fv, mx_key[1], T)
            val -= cf * ASx(uncode([fv]),mx_key[1])
        else:
            cprd_val += cf * coprod_back(ASx(uncode(cd), mx_key[1] - 1), T) * single_coprod(fv, 1, T)
            val -= cf * ASx(uncode(cd), mx_key[1] - 1) * ASx(uncode([fv]),1)
        print(val)
    return cprd_val


if __name__ == "__main__":
    
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
