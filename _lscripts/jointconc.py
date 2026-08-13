"""Is g=log Lambda JOINTLY concave on the partition cone?
If yes (+ g(0)=0) => subadditive => rho<=1, a clean unconditional proof.

Test discrete concavity in all directions:
 (D1) axis e_i:        Lam(c)^2 >= Lam(c-e_i)Lam(c+e_i)                 [redo, all c]
 (D2) e_i+e_j:         Lam(c)^2 >= Lam(c-e_i-e_j)Lam(c+e_i+e_j)
 (D3) e_i-e_j:         Lam(c)^2 >= Lam(c-e_i+e_j)Lam(c+e_i-e_j)   (i!=j)
 (MID) general 2-point midpoint concavity on integer points that stay partitions:
        for partitions p<=q with (p+q) even coordinatewise,
        Lam((p+q)/2)^2 >= Lam(p)Lam(q)
Also directly retest subadditivity rho<=1 for reference.
"""
from math import factorial
from itertools import product as iproduct, combinations
from fractions import Fraction
from schubmult import uncode


def sigma_oneline(code, M):
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    if not code:
        return list(range(1, M+1))
    ol = list(uncode(code))
    return ol + list(range(len(ol)+1, M+1))


def Lam(code, M):
    ol = sigma_oneline(code, M)
    p = 1
    for i in range(M):
        p *= ol[i] - (code[i] if i < len(code) else 0)
    return factorial(M) // p


def is_part(t):
    return all(t[i] >= t[i+1] for i in range(len(t)-1)) and all(x >= 0 for x in t)


def main():
    M = 8
    maxpart, maxlen = 4, 4
    parts = [t for t in iproduct(range(maxpart+1), repeat=maxlen) if is_part(t)]

    def L(t):
        return Lam(tuple(t), M)

    d1f=d2f=d3f=0; d1c=d2c=d3c=0
    d3worst=None; d2worst=None
    for c in parts:
        if max(c) > M-2:
            continue
        for i in range(maxlen):
            cm=list(c); cm[i]-=1; cp=list(c); cp[i]+=1
            if is_part(cm) and is_part(cp) and max(cp)<=M-1:
                d1c+=1
                if L(c)**2 < L(cm)*L(cp): d1f+=1
        for i,j in combinations(range(maxlen),2):
            # D2 e_i+e_j
            cm=list(c); cp=list(c)
            cm[i]-=1; cm[j]-=1; cp[i]+=1; cp[j]+=1
            if is_part(cm) and is_part(cp) and max(cp)<=M-1:
                d2c+=1
                if L(c)**2 < L(cm)*L(cp):
                    d2f+=1
                    r=Fraction(L(cm)*L(cp),L(c)**2)
                    if d2worst is None or r>d2worst[0]: d2worst=(r,c,i,j)
            # D3 e_i-e_j  and  e_j-e_i
            cm=list(c); cp=list(c)
            cm[i]-=1; cm[j]+=1; cp[i]+=1; cp[j]-=1
            if is_part(cm) and is_part(cp) and max(cp)<=M-1:
                d3c+=1
                if L(c)**2 < L(cm)*L(cp):
                    d3f+=1
                    r=Fraction(L(cm)*L(cp),L(c)**2)
                    if d3worst is None or r>d3worst[0]: d3worst=(r,c,i,j)
    print(f"D1 axis e_i:      {d1f} fail / {d1c}")
    print(f"D2 e_i+e_j:       {d2f} fail / {d2c}" + (f"  worst ratio {float(d2worst[0]):.4f} c={d2worst[1]} (i,j)=({d2worst[2]},{d2worst[3]})" if d2worst else ""))
    print(f"D3 e_i-e_j:       {d3f} fail / {d3c}" + (f"  worst ratio {float(d3worst[0]):.4f} c={d3worst[1]} (i,j)=({d3worst[2]},{d3worst[3]})" if d3worst else ""))

    # MID midpoint concavity over partition pairs of equal parity
    midf=midc=0; midworst=None
    small = [t for t in parts if max(t)<=3]
    for p in small:
        for q in small:
            if any((p[k]+q[k])%2 for k in range(maxlen)):
                continue
            mid=tuple((p[k]+q[k])//2 for k in range(maxlen))
            if not is_part(mid) or not is_part(p) or not is_part(q):
                continue
            midc+=1
            if L(mid)**2 < L(p)*L(q):
                midf+=1
                r=Fraction(L(p)*L(q),L(mid)**2)
                if midworst is None or r>midworst[0]: midworst=(r,p,q)
    print(f"MID midpoint:     {midf} fail / {midc}" + (f"  worst {float(midworst[0]):.4f} p={midworst[1]} q={midworst[2]}" if midworst else ""))


if __name__ == "__main__":
    main()
