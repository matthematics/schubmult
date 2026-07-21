"""Verify the clean dichotomy and the 'exactly two reflections' claim.

Moved case: u NOT->_k w, u(k)!=w(k). pi=u^{-1}w, C=cycle(k), B=max(C).

CLAIMS to check (0 failures expected at n=5):
 (D1) B>k  <=>  u NOT->_{k-1} w   (so DIAG & STRAIGHT vanish)  AND  u NOT->_{k-2} w (QUANT=0).
 (D2) B=k  <=>  the cycle C has a unique index >=k (=k).  (RECURR / QUANT world.)
 (D3) In the B>k world: for every p in 1..k, the contributing position-k reflections
      t_{rk} (wt_{rk} lup w and u->_{k-1} wt_{rk}) with c_{k-1}(p-1;u,wt_{rk})!=0 are
      EXACTLY {B} and {one r*<k, r* in C}, with:
          sign cancel, equal E-signature (deg,nv,Pset), equal q-weight kappa*q_{k-1}.
      i.e. #contributing up-neighbors ==1 (==B) and #contributing down-neighbors ==1.
 (D4) In the B>k world F (whole four-term RHS) == 0 (sanity).

Run: conda activate schubmult_312 && python _lscripts/moved_dichotomy.py 5
"""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, expand_func, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym

y=GeneratingSet("y"); z=GeneratingSet("z"); q_gs=GeneratingSet("q")
def q_ab(a,b): return prod([q_gs[s] for s in range(a,b)])
def val(p,i): return p[i-1]
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
def Pset(u,w,m): return tuple(sorted(val(u,i) for i in range(1,m+1) if val(u,i)==val(w,i)))
def sig(u,w,m,a,Rm):
    if m<0: return None
    w=Permutation(w)
    if w not in Rm: return None
    qm=Rm[w]; Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    if deg<0 or deg>nv: return None
    if nv==0: poly=S.One if deg==0 else S.Zero
    else:
        yv=[y[v] for v in Pm]; nc=nv+1-deg; zv=[z[i] for i in range(1,nc+1)]
        poly=expand_func(FactorialElemSym(deg,nv,yv,zv))
    return (expand(qm), poly, (deg,nv,Pm))
def orbit(perm,start,N):
    c=[start]; cur=perm[start-1]
    while cur!=start: c.append(cur); cur=perm[cur-1]
    return c

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    d1a=0; d1b=0; d3fail=0; d4fail=0
    bgtk_inst=0; beqk_inst=0
    down_count_hist=Counter(); up_count_hist=Counter()
    pair_ok=0
    for u in perms:
        u=Permutation(u)
        for k in range(2,n):
            Rk=enumerate_pieri(u,k,N); Rkm1=enumerate_pieri(u,k-1,N)
            Rkm2=enumerate_pieri(u,k-2,N) if k-2>=0 else {u:S.One}
            cset=set(Rk)|set(Rkm1)|set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1,N+1):
                    if pp!=k: cset.add(v0.swap(k-1,pp-1))
            for w in cset:
                w=Permutation(w)
                if val(u,k)==val(w,k): continue
                if w in Rk: continue
                pi=(~u)*w
                C=orbit(pi,k,N); B=max(C)
                inkm1 = w in Rkm1; inkm2 = w in Rkm2
                if B>k:
                    bgtk_inst+=1
                    # D1: B>k => not ->_{k-1} and not ->_{k-2}
                    if inkm1: d1a+=1
                    if inkm2: d1b+=1
                    # D3: per p exactly one up (=B), one down (r*<k in C)
                    for p in range(1,k+1):
                        ups=[]; downs=[]
                        for pp in range(1,N+1):
                            if pp==k: continue
                            vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                            lo,hi=min(pp,k),max(pp,k)
                            if dl==1: kap=S.One
                            elif dl==1-2*(hi-lo): kap=q_ab(lo,hi)
                            else: continue
                            s=sig(u,vv,k-1,p-1,Rkm1)
                            if not s: continue
                            wt=expand(kap*s[0]); Esig=s[2]; poly=s[1]
                            if pp>k: ups.append((pp,wt,Esig,poly))
                            else: downs.append((pp,wt,Esig,poly))
                        up_count_hist[len(ups)]+=1; down_count_hist[len(downs)]+=1
                        if len(ups)==1 and len(downs)==1:
                            (up_pp,uw,ue,up_poly)=ups[0]; (dn_pp,dw,de,dn_poly)=downs[0]
                            ok = (up_pp==B) and (dn_pp<k) and (dn_pp in C) and (uw==dw) and (ue==de) \
                                 and expand(up_poly-dn_poly)==0
                            if ok: pair_ok+=1
                            else: d3fail+=1
                        elif len(ups)==0 and len(downs)==0:
                            pass
                        else:
                            d3fail+=1
                        # D4: whole F==0
                        F=S.Zero
                        s=sig(u,w,k-1,p-1,Rkm1)
                        if s: F+=(y[val(w,k)]-z[k-p+1])*s[1]
                        s=sig(u,w,k-1,p,Rkm1)
                        if s: F+=s[1]*s[0]/s[0] if False else expand(s[1])  # straight has weight qm too; keep full
                        # rebuild F properly with weights:
                        F=S.Zero
                        s=sig(u,w,k-1,p-1,Rkm1)
                        if s: F+=expand((y[val(w,k)]-z[k-p+1])*s[1]*s[0])
                        s=sig(u,w,k-1,p,Rkm1)
                        if s: F+=expand(s[1]*s[0])
                        s=sig(u,w,k-2,p-2,Rkm2)
                        if s: F+=expand(q_gs[k-1]*s[0]*s[1])
                        for pp in range(1,N+1):
                            if pp==k: continue
                            vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                            lo,hi=min(pp,k),max(pp,k)
                            if dl==1: kap=S.One; sgn=1 if pp>k else -1
                            elif dl==1-2*(hi-lo): kap=q_ab(lo,hi); sgn=1 if pp>k else -1
                            else: continue
                            s=sig(u,vv,k-1,p-1,Rkm1)
                            if s: F+=expand(sgn*kap*s[0]*s[1])
                        if expand(F)!=0: d4fail+=1
                else:
                    beqk_inst+=1
    print(f"B>k instances: {bgtk_inst} | B=k instances: {beqk_inst}")
    print(f"D1 fails (B>k but ->_{{k-1}}): {d1a} | (B>k but ->_{{k-2}}): {d1b}")
    print(f"D3: down-count hist {dict(down_count_hist)} | up-count hist {dict(up_count_hist)}")
    print(f"D3 pair-ok: {pair_ok} | D3 fails: {d3fail}")
    print(f"D4 (F!=0 in B>k world): {d4fail}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
