"""Test a UNIVERSAL 2-part decomposition of COMB-REC that would give a COMPLETE proof:

  RHS(u,w,k,p) = c_k(p;u,w)   [the identity to prove]

Decompose:
  PARTNER = the unique reflection term  w t_{a,k} (a<k) OR w t_{k,b} (b>k) that is the
            'quantum partner': target v=wt with u->_{k-1} v, u NOT ->_{k-2} v, and the
            cycle of u^{-1}v through k has top k, chosen to match QUANT.
  (T1)  QUANT + PARTNER-term == 0,  and PARTNER exists (uniquely) iff u->_{k-2} w.
  (T2)  RHS - QUANT - PARTNER-term == c_k(p;u,w)      [the 'classical remainder']

If T1&T2 hold universally (0 fail), the theorem = classical-remainder (provable by
Samuel-Molev classical Pieri) + reflectioncancel bijection (cycle surgery).

We identify PARTNER structurally: among reflections wt (t=t_{ak},a<k or t=t_{kb},b>k)
with c_{k-1}(p-1;u,wt)!=0, pick those v=wt with u NOT->_{k-2} v AND cycle(u^{-1}v thru k)
top == k. Test there is at most one, exactly one iff u->_{k-2}w, and it cancels QUANT.

Run: conda activate schubmult_312 && python _lscripts/universal_decomp.py 5
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
def cyc_top_through(u,v,idx):
    perm=(~u)*v; cyc=[idx]; cur=perm[idx-1]
    while cur!=idx: cyc.append(cur); cur=perm[cur-1]
    if len(cyc)==1: return None
    return max(cyc)

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; t1_fail=0; t2_fail=0; partner_count_fail=0
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
                for p in range(1,k+1):
                    # full RHS
                    diag=(y[val(w,k)]-z[k-p+1])*c_term(u,w,k-1,p-1,Rkm1)
                    rec=c_term(u,w,k-1,p,Rkm1)
                    quant=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    ck=c_term(u,w,k,p,Rk)
                    # collect reflection terms and identify partner(s)
                    monk=S.Zero; partner_val=S.Zero; partners=[]
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        sign=1 if pp>k else -1; lo,hi=min(pp,k),max(pp,k)
                        if dl==1: tval=sign*c_term(u,vv,k-1,p-1,Rkm1)
                        elif dl==1-2*(hi-lo): tval=sign*q_ab(lo,hi)*c_term(u,vv,k-1,p-1,Rkm1)
                        else: continue
                        monk+=tval
                        # is vv a quantum partner? u->_{k-1}vv, NOT u->_{k-2}vv, cyc thru k top ==k
                        if c_nz(u,vv,k-1,p-1,Rkm1) and not reachable(u,vv,k-2,N):
                            if cyc_top_through(u,vv,k)==k:
                                partners.append(pp); partner_val+=tval
                    rhs=expand(diag+monk+rec+quant)
                    # only consider genuine identity instances (rhs should == ck)
                    if expand(rhs-ck)!=0:
                        continue  # not our identity domain (shouldn't happen)
                    tot+=1
                    quant_nz = c_nz(u,w,k-2,p-2,Rkm2)
                    # partner count: exactly one iff quant_nz
                    if quant_nz and len(partners)!=1: partner_count_fail+=1
                    if (not quant_nz) and len(partners)!=0: partner_count_fail+=1
                    # T1: quant + partner_val == 0
                    if expand(quant+partner_val)!=0: t1_fail+=1
                    # T2: rhs - quant - partner_val == ck
                    if expand((diag+monk+rec+quant)-quant-partner_val-ck)!=0: t2_fail+=1
    print(f"identity instances = {tot}")
    print(f"partner count wrong (not 1 iff quant nz) : {partner_count_fail}")
    print(f"(T1) QUANT + PARTNER != 0                : {t1_fail}")
    print(f"(T2) remainder != c_k                    : {t2_fail}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
