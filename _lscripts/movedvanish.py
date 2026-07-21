"""Classify the moved-vanishing case and find the exact cancellation pairing.

Moved-vanishing: u(k)!=w(k) and u NOT->_k w. Must show COMB-REC RHS = 0.
Split by whether u->_{k-1} w:
  3a: u->_{k-1} w  (DIAG, REC nonzero)
  3b: NOT u->_{k-1} w (DIAG=REC=0; RHS = MONK + QUANT)

For each instance, list every nonzero term with its (label, value, target-perm, q-weight),
and try to find the sign-reversing pairing. Print term-count patterns per subcase.

Also directly TEST two reduction hypotheses:
  (H3b) In 3b: MONK + QUANT == 0   (reflectioncancel without u(k)=w(k)).
  (H3a) In 3a: split reflections into length-INCREASING (raise) vs length-DECREASING (drop).
        Test: DIAG + REC + (all RAISE reflections) == 0     [ordinary/classical part]
          and (all DROP reflections) + QUANT == 0           [quantum part]

Run: conda activate schubmult_312 && python _lscripts/movedvanish.py 5
"""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, expand_func, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym

y = GeneratingSet("y"); z = GeneratingSet("z"); q_gs = GeneratingSet("q")
def q_ab(a,b): return prod([q_gs[s] for s in range(a,b)])
def val(perm,pos): return perm[pos-1]

def enumerate_pieri(u,k,N):
    u=Permutation(u); results={u:S.One}; stack=[(u,frozenset(),N+1,S.One,u.inv)]; seen=set()
    while stack:
        perm,used_a,last_b,qw,clen=stack.pop(); key=(perm,used_a,last_b)
        if key in seen: continue
        seen.add(key)
        for a in range(1,k+1):
            if a in used_a: continue
            for b in range(k+1,N+1):
                if b>last_b: continue
                nperm=perm.swap(a-1,b-1); d=nperm.inv-clen
                if d==1: nqw=qw
                elif d==-2*(b-a)+1: nqw=qw*q_ab(a,b)
                else: continue
                results.setdefault(nperm,nqw); stack.append((nperm,used_a|{a},b,nqw,nperm.inv))
    return results

def Pset(u,w,m): return sorted(val(u,i) for i in range(1,m+1) if val(u,i)==val(w,i))
def c_term(u,w,m,a,Rm):
    if m<0: return S.Zero
    w=Permutation(w)
    if w not in Rm: return S.Zero
    qm=Rm[w]; Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    if deg<0 or deg>nv: return S.Zero
    if nv==0: return qm if deg==0 else S.Zero
    yv=[y[v] for v in Pm]; nc=nv+1-deg; zv=[z[i] for i in range(1,nc+1)]
    return qm*expand_func(FactorialElemSym(deg,nv,yv,zv))

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; fail3b=0; fail3a_ord=0; fail3a_q=0
    sub=Counter()
    for u in perms:
        u=Permutation(u)
        for k in range(2,n):
            Rk=enumerate_pieri(u,k,N); Rkm1=enumerate_pieri(u,k-1,N)
            Rkm2=enumerate_pieri(u,k-2,N) if k-2>=0 else {u:S.One}
            cands=set(Rk)|set(Rkm1)|set(Rkm2)
            for v in list(Rkm1):
                for pp in range(1,N+1):
                    if pp!=k: cands.add(v.swap(k-1,pp-1))
            for w in cands:
                w=Permutation(w)
                if val(u,k)==val(w,k): continue
                if w in Rk: continue  # u NOT ->_k w
                km1 = w in Rkm1
                for p in range(1,k+1):
                    diag=(y[val(w,k)]-z[k-p+1])*c_term(u,w,k-1,p-1,Rkm1)
                    rec=c_term(u,w,k-1,p,Rkm1)
                    quant=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    raise_sum=S.Zero; drop_sum=S.Zero; monk=S.Zero
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        sign=1 if pp>k else -1; lo,hi=min(pp,k),max(pp,k)
                        if dl==1:
                            t=sign*c_term(u,vv,k-1,p-1,Rkm1); raise_sum+=t; monk+=t
                        elif dl==1-2*(hi-lo):
                            t=sign*q_ab(lo,hi)*c_term(u,vv,k-1,p-1,Rkm1); drop_sum+=t; monk+=t
                    rhs=expand(diag+monk+rec+quant)
                    if rhs!=0:
                        continue  # not a valid instance? shouldn't happen
                    tot+=1
                    subcase='3a' if km1 else '3b'
                    nz=(expand(diag)!=0,expand(monk)!=0,expand(rec)!=0,expand(quant)!=0)
                    sub[(subcase,)+nz]+=1
                    if subcase=='3b':
                        if expand(monk+quant)!=0: fail3b+=1
                    else:
                        if expand(diag+rec+raise_sum)!=0: fail3a_ord+=1
                        if expand(drop_sum+quant)!=0: fail3a_q+=1
    print(f"moved-vanishing instances = {tot}")
    print(f"(H3b) MONK+QUANT != 0 in 3b            : {fail3b}")
    print(f"(H3a-ord) DIAG+REC+RAISE != 0 in 3a    : {fail3a_ord}")
    print(f"(H3a-q)   DROP+QUANT != 0 in 3a        : {fail3a_q}")
    print("subcase (DIAG,MONK,REC,QUANT) nonzero patterns:")
    for pt,c in sorted(sub.items(),key=lambda x:-x[1]): print(f"    {pt}: {c}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
