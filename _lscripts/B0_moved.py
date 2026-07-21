"""B0-moved case: u NOT->_k w, u(k)!=w(k). Test whether the vanishing there also
reduces to pieriknotkn1 / leveldichotomy, or genuinely needs computation.

The COMB-REC RHS must be 0. Terms:
  DIAG = (y_{w(k)}-z_{k-p+1}) c_{k-1}(p-1;u,w)
  MONK = sum_{pp!=k} sign(pp-k) kappa c_{k-1}(p-1;u,w t_{k,pp})
  REC  = c_{k-1}(p;u,w)
  QUANT= q_{k-1} c_{k-2}(p-2;u,w)

Sub-questions when u(k)!=w(k) and u NOT->_k w:
  (Q1) Is c_{k-1}(p-1;u,w) [DIAG factor] and c_{k-1}(p;u,w) [REC] ever nonzero?
       i.e. can u ->_{k-1} w when u NOT ->_k w and u(k)!=w(k)?
  (Q2) Which MONK neighbors are nonzero, and do they pair with DIAG/REC/QUANT?

Print the breakdown of which terms are nonzero, to see the cancellation pattern.
Run: conda activate schubmult_312 && python _lscripts/B0_moved.py 5
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
    total=0; failsum=0; patt=Counter()
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
                if val(u,k)==val(w,k): continue      # moved only
                if w in Rk: continue                  # u NOT ->_k w only  (B0 moved)
                for p in range(1,k+1):
                    diag=(y[val(w,k)]-z[k-p+1])*c_term(u,w,k-1,p-1,Rkm1)
                    rec=c_term(u,w,k,p,Rk) if False else S.Zero  # not needed; REC uses level k? No, c_{k-1}(p)
                    rec=c_term(u,w,k-1,p,Rkm1)
                    quant=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    monk=S.Zero
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        sign=1 if pp>k else -1; lo,hi=min(pp,k),max(pp,k)
                        if dl==1: monk+=sign*c_term(u,vv,k-1,p-1,Rkm1)
                        elif dl==1-2*(hi-lo): monk+=sign*q_ab(lo,hi)*c_term(u,vv,k-1,p-1,Rkm1)
                    total_rhs=expand(diag+monk+rec+quant)
                    total+=1
                    if total_rhs!=0: failsum+=1
                    patt[(expand(diag)!=0,expand(monk)!=0,expand(rec)!=0,expand(quant)!=0)]+=1
    print(f"B0-moved instances = {total}")
    print(f"RHS != 0 failures  = {failsum}")
    print("(DIAG,MONK,REC,QUANT) nonzero pattern:")
    for pt,c in sorted(patt.items(),key=lambda x:-x[1]): print(f"    {pt}: {c}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
