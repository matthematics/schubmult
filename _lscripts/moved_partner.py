"""MOVED-B0 QUANT!=0: find the QUANT partner by Q-MONOMIAL matching, and test the
complete decomposition:
   (A) There is a UNIQUE reflection term T* whose contribution, added to QUANT, has
       q-monomial cancellation, such that QUANT + T* == 0.
   (B) The remaining terms (all reflection terms except T*, plus DIAG + REC) sum to 0
       and constitute the CLASSICAL Pieri vanishing identity (provable via Samuel-Molev
       by setting the relevant q's; more precisely they form a q-weighted classical
       vanishing).

Partner identification: T* is the reflection term with v=w t (t=t_{ak} a<k or t_{kb} b>k)
such that v = the unique level-(k-1) neighbour with  u ->_{k-1} v, u NOT ->_{k-2} v, and
v t_{k-1,k} = w' for some relation making q-weight = q_{k-1} q_{k-2}(u,w).

Simplest robust test: among all reflection terms, find those T with q-weight(T) ==
q_{k-1} * q_{k-2}(u,w) as a q-monomial AND c-poly(T) == c_{k-2}(p-2;u,w) poly part.
Check exactly one, and QUANT + T == 0, and remainder == 0.

Run: conda activate schubmult_312 && python _lscripts/moved_partner.py 5
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

def reachable(u,w,k,N): return Permutation(w) in enumerate_pieri(u,k,N)
def Pset(u,w,m): return sorted(val(u,i) for i in range(1,m+1) if val(u,i)==val(w,i))
def qmono(u,w,m,Rm):
    w=Permutation(w)
    return Rm.get(w, None)
def c_term(u,w,m,a,Rm):
    if m<0: return S.Zero
    w=Permutation(w)
    if w not in Rm: return S.Zero
    qm=Rm[w]; Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    if deg<0 or deg>nv: return S.Zero
    if nv==0: return qm if deg==0 else S.Zero
    yv=[y[v] for v in Pm]; nc=nv+1-deg; zv=[z[i] for i in range(1,nc+1)]
    return qm*expand_func(FactorialElemSym(deg,nv,yv,zv))
def c_nz(u,w,m,a,Rm):
    if m<0: return False
    w=Permutation(w)
    if w not in Rm: return False
    Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    return 0<=deg<=nv

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; partner_fail=0; cancel_fail=0; remainder_fail=0
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
                if w in Rk: continue
                for p in range(1,k+1):
                    if not c_nz(u,w,k-2,p-2,Rkm2): continue
                    tot+=1
                    diag=(y[val(w,k)]-z[k-p+1])*c_term(u,w,k-1,p-1,Rkm1)
                    rec=c_term(u,w,k-1,p,Rkm1)
                    quant=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    # QUANT partner: reflection term with v = w t, c_{k-1}(p-1;u,v) nonzero,
                    # u NOT ->_{k-2} v, and v t_{k-1,k} == w (the reflectioncancel signature),
                    # OR more generally q-weight(v)= q_{k-1}(u,v) with the adjacency.
                    partner_terms=[]
                    reflection_total=S.Zero
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: tval=(1 if pp>k else -1)*c_term(u,vv,k-1,p-1,Rkm1)
                        elif dl==1-2*(hi-lo): tval=(1 if pp>k else -1)*q_ab(lo,hi)*c_term(u,vv,k-1,p-1,Rkm1)
                        else: continue
                        reflection_total+=tval
                        # partner signature: v t_{k-1,k} == w  (adjacent recovers w)
                        back=vv.swap(k-2,k-1)
                        if back==w and c_nz(u,vv,k-1,p-1,Rkm1) and not reachable(u,vv,k-2,N):
                            partner_terms.append((pp,tval))
                    if len(partner_terms)!=1: partner_fail+=1; continue
                    pp,tval=partner_terms[0]
                    if expand(quant+tval)!=0: cancel_fail+=1
                    remainder=expand(diag+rec+reflection_total-tval)
                    if remainder!=0: remainder_fail+=1
    print(f"moved-B0 QUANT!=0 instances = {tot}")
    print(f"partner (v t_{{k-1,k}}=w) count != 1 : {partner_fail}")
    print(f"QUANT + partner != 0                : {cancel_fail}")
    print(f"remainder (rest) != 0               : {remainder_fail}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
