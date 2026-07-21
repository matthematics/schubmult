"""Examine the 20 moved-B0 QUANT!=0 cases where v t_{k-1,k}=w partner signature fails.
Print full detail: u,w,p,k, all reflection terms with their v, dl-type, whether
v t_{k-1,k}==w, u->_{k-1}v, u->_{k-2}v, and which term actually cancels QUANT.
"""
import sys
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

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n)); shown=0
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
                    quant=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    partner_terms=[]
                    all_terms=[]
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: tval=(1 if pp>k else -1)*c_term(u,vv,k-1,p-1,Rkm1); typ="RAISE"
                        elif dl==1-2*(hi-lo): tval=(1 if pp>k else -1)*q_ab(lo,hi)*c_term(u,vv,k-1,p-1,Rkm1); typ="DROP"
                        else: continue
                        if expand(tval)==0: continue
                        back=vv.swap(k-2,k-1)
                        sig = (back==w and c_nz(u,vv,k-1,p-1,Rkm1) and not reachable(u,vv,k-2,N))
                        cancels = (expand(quant+tval)==0)
                        all_terms.append((pp,typ,tuple(vv),sig,cancels,c_nz(u,vv,k-1,p-1,Rkm1),reachable(u,vv,k-2,N),tuple(back)==tuple(w)))
                        if sig: partner_terms.append(pp)
                    if len(partner_terms)==1: continue  # not a failure
                    if shown>=20: return
                    shown+=1
                    print(f"\n#{shown} u={tuple(u)} p={p} k={k} w={tuple(w)}  (partner_count={len(partner_terms)})")
                    print(f"   u(k-1)={val(u,k-1)} u(k)={val(u,k)} w(k-1)={val(w,k-1)} w(k)={val(w,k)}")
                    for (pp,typ,vv,sig,cancels,uk1,uk2,adj) in all_terms:
                        print(f"   pp={pp} {typ} v={vv} sig={sig} cancels_QUANT={cancels} u->k-1={uk1} u->k-2={uk2} vt_(k-1,k)=w?{adj}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
