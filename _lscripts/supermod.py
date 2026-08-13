"""Find the structural fact that ACTUALLY implies Lam(a+b)<=Lam(a)Lam(b).

Componentwise join/meet: (a v b)_i=max, (a ^ b)_i=min.  Note a+b = (a v b)+(a ^ b).

Test:
 (S+) log-supermodular:  Lam(a)Lam(b) <= Lam(a v b) Lam(a ^ b)
 (S-) log-submodular:    Lam(a)Lam(b) >= Lam(a v b) Lam(a ^ b)
 (SUB) subadditive:      Lam(a+b) <= Lam(a) Lam(b)                 [known true; reref]
 (SUPERADD across meet/join sum):  Lam((avb)+(a^b)) <= Lam(a v b) Lam(a ^ b)?  [= SUB on avb,a^b]

KEY candidate chain for a PROOF of SUB:
   If  Lam(x+y) <= Lam(x)Lam(y)  reduces to the case of COMPARABLE x<=y (nested partitions),
   because a v b >= a ^ b are comparable and Lam(a)Lam(b) and Lam(avb)Lam(a^b) compare via (S).
   So test the COMPARABLE (nested) case separately:
 (NEST) for partitions x<=y (coordinatewise), Lam(x+y) <= Lam(x)Lam(y)?  worst ratio?
   and whether SUB for comparable pairs + (S) gives SUB in general.
"""
from math import factorial
from itertools import product as iproduct
from fractions import Fraction
from schubmult import uncode


def sig(code, M):
    code=list(code)
    while code and code[-1]==0: code.pop()
    if not code: return list(range(1,M+1))
    ol=list(uncode(code)); return ol+list(range(len(ol)+1,M+1))

def Lam(code,M):
    ol=sig(code,M); p=1
    for i in range(M): p*= ol[i]-(code[i] if i<len(code) else 0)
    return factorial(M)//p

def is_part(t):
    return all(t[i]>=t[i+1] for i in range(len(t)-1)) and all(x>=0 for x in t)

def main():
    M=8; maxpart,maxlen=4,4
    parts=[t for t in iproduct(range(maxpart+1),repeat=maxlen) if is_part(t)]
    def L(t): return Lam(tuple(t),M)
    sp=sm=0; spc=0; spworst=None; smworst=None
    subf=0; subc=0
    nestf=0; nestc=0; nestworst=None
    for a in parts:
        for b in parts:
            jn=tuple(max(a[i],b[i]) for i in range(maxlen))
            mt=tuple(min(a[i],b[i]) for i in range(maxlen))
            # both a,b partitions => join & meet are partitions
            if is_part(jn) and is_part(mt) and max(jn)<=M-1:
                spc+=1
                lhs=L(a)*L(b); rhs=L(jn)*L(mt)
                if lhs>rhs:
                    sp+=1
                    r=Fraction(lhs,rhs)
                    if spworst is None or r>spworst[0]: spworst=(r,a,b)
                if lhs<rhs:
                    sm+=1
                    r=Fraction(rhs,lhs)
                    if smworst is None or r>smworst[0]: smworst=(r,a,b)
            g=tuple(a[i]+b[i] for i in range(maxlen))
            if max(g)<=M-1:
                subc+=1
                if L(g)>L(a)*L(b): subf+=1
            # nested case
            if all(a[i]<=b[i] for i in range(maxlen)) and max(a[i]+b[i] for i in range(maxlen))<=M-1:
                nestc+=1
                r=Fraction(L(g),L(a)*L(b))
                if L(g)>L(a)*L(b): nestf+=1
                if nestworst is None or r>nestworst[0]: nestworst=(r,a,b)
    print(f"(S+) log-supermodular  Lam(a)Lam(b)<=Lam(avb)Lam(a^b): {sp} fail /{spc}" + (f"  worst {float(spworst[0]):.3f} a={spworst[1]} b={spworst[2]}" if spworst else ""))
    print(f"(S-) log-submodular    Lam(a)Lam(b)>=Lam(avb)Lam(a^b): {sm} fail /{spc}" + (f"  worst {float(smworst[0]):.3f} a={smworst[1]} b={smworst[2]}" if smworst else ""))
    print(f"(SUB) subadditive      Lam(a+b)<=Lam(a)Lam(b):         {subf} fail /{subc}")
    print(f"(NEST) nested a<=b      Lam(a+b)<=Lam(a)Lam(b):         {nestf} fail /{nestc}" + (f"  sup ratio {float(nestworst[0]):.4f} a={nestworst[1]} b={nestworst[2]}" if nestworst else ""))


if __name__=="__main__":
    main()
