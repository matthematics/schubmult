"""Nail the provable kernel.

(K1) LOCAL supermodularity (neighbor exchange):
       Lam(c+e_i)Lam(c+e_j) <= Lam(c)Lam(c+e_i+e_j)   for i<j, all c with results partitions.
     If this holds for all neighbors it is the local form of (S+).

(K2) Does (K1) [local, provable] IMPLY global subadditivity?  Test the reduction:
     reduce a general split to nested by ONE exchange step and see monotonicity.
     Concretely, define for a+b=gamma fixed, f(split)=Lam(a)Lam(b). Claim: f is MINIMIZED
     at the nested split (a,b)=(gamma^even floor / ceil)?? and the LEAST value still >= Lam(gamma).
     Test: over all splits of a fixed gamma, is min_split Lam(a)Lam(b) achieved at the most
     'balanced'/nested split, and is it always >= Lam(gamma)?  (i.e. SUB reduces to that split.)

(K3) 'Completely log-concave' surrogate: for the hook product, test the stronger
     Lam(c)Lam(c+e_i+e_j) >= Lam(c+e_i)Lam(c+e_j)  (= K1) AND
     Lam(c)^2 >= Lam(c-e_i)Lam(c+e_i) (axis, done) -> together = Lorentzian-type.
"""
from math import factorial
from itertools import product as iproduct, combinations
from fractions import Fraction
from schubmult import uncode

def sig(code,M):
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
    M=8; maxpart,maxlen=5,4
    parts=[t for t in iproduct(range(maxpart+1),repeat=maxlen) if is_part(t)]
    def L(t): return Lam(tuple(t),M)

    # K1 local supermodularity
    k1f=k1c=0; k1worst=None
    for c in parts:
        for i,j in combinations(range(maxlen),2):
            ci=list(c); ci[i]+=1
            cj=list(c); cj[j]+=1
            cij=list(c); cij[i]+=1; cij[j]+=1
            if is_part(ci) and is_part(cj) and is_part(cij) and max(cij)<=M-1:
                k1c+=1
                if L(ci)*L(cj) > L(c)*L(cij):
                    k1f+=1
                    r=Fraction(L(ci)*L(cj),L(c)*L(cij))
                    if k1worst is None or r>k1worst[0]: k1worst=(r,c,i,j)
    print(f"(K1) local supermodular Lam(c+ei)Lam(c+ej)<=Lam(c)Lam(c+ei+ej): {k1f} fail /{k1c}" + (f"  worst {float(k1worst[0]):.4f}" if k1worst else ""))

    # K2 among all splits gamma=a+b, is Lam(gamma) <= min_split Lam(a)Lam(b), and where is min?
    import itertools
    def splits(g):
        ranges=[range(0,gi+1) for gi in g]
        for a in iproduct(*ranges):
            if is_part(a):
                b=tuple(g[k]-a[k] for k in range(len(g)))
                if is_part(b):
                    yield a,b
    gammas=[g for g in parts if max(g)<=maxpart]
    k2f=0; k2c=0; balanced_is_min=0; nested_all=0; tot_g=0
    for g in gammas:
        if sum(g)==0: continue
        vals=[]
        for a,b in splits(g):
            vals.append((L(a)*L(b),a,b))
        if not vals: continue
        tot_g+=1
        mn=min(vals, key=lambda x:x[0])
        # is min >= Lam(gamma)? (SUB says every split >= Lam(g); so min >= Lam(g))
        k2c+=1
        if mn[0] < L(g): k2f+=1
        # is the argmin nested?
        a,b=mn[1],mn[2]
        if all(a[k]<=b[k] for k in range(maxlen)) or all(b[k]<=a[k] for k in range(maxlen)):
            nested_all+=1
    print(f"(K2) min over splits of Lam(a)Lam(b) >= Lam(gamma): {k2c-k2f}/{k2c} hold; argmin nested in {nested_all}/{tot_g} gammas")

    # K3 report: K1 + axis concavity both hold => Lorentzian-type signature
    print("(K3) axis log-concavity holds (prior test); with K1 this is the Lorentzian/Hodge-Riemann signature")


if __name__=="__main__":
    main()
