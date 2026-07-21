"""For moved instances where DIAG+REC != Monk-at-a0, find the ACTUAL reflection w'=w t_{k,pp}
in the weight-M0 group (M0=q_{k-1}(u,w)) with c_{k-1}(p-1;u,w')*(sign,kappa) = -(DIAG+REC),
and report: pp (vs k), its type (raise/drop), P_{k-1}(u,w') vs P_{k-1}(u,w), and whether
P_{k-1}(u,w') == P_{k-1}(u,w) ∪ {w(k)}.
Goal: characterize the DIAG+REC partner reflection correctly for the write-up.
Run: conda activate schubmult_312 && python _lscripts/moved_findpartner.py 5
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
    tot=0; nopartner=0; multi=0
    ppkind=Counter(); psetrel=Counter()
    for u in perms:
        u=Permutation(u)
        for k in range(2,n):
            Rk=enumerate_pieri(u,k,N); Rkm1=enumerate_pieri(u,k-1,N)
            Rkm2=enumerate_pieri(u,k-2,N) if k-2>=0 else {u:S.One}
            cands=set(Rk)|set(Rkm1)|set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1,N+1):
                    if pp!=k: cands.add(v0.swap(k-1,pp-1))
            for w in cands:
                w=Permutation(w)
                if val(u,k)==val(w,k): continue
                if w in Rk: continue
                for p in range(1,k+1):
                    diag=(y[val(w,k)]-z[k-p+1])*c_term(u,w,k-1,p-1,Rkm1)
                    rec=c_term(u,w,k-1,p,Rkm1)
                    dr=expand(diag+rec)
                    if dr==0: continue
                    tot+=1
                    M0=expand(Rkm1[w]) if w in Rkm1 else None
                    found=[]
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: extra=S.One; sgn=(1 if pp>k else -1); ty="RAISE"
                        elif dl==1-2*(hi-lo): extra=q_ab(lo,hi); sgn=(1 if pp>k else -1); ty="DROP"
                        else: continue
                        if vv not in Rkm1: continue
                        # weight of summand
                        wt=expand(extra*Rkm1[vv])
                        if wt!=M0: continue
                        summand=sgn*extra*c_term(u,vv,k-1,p-1,Rkm1)
                        if expand(dr+summand)==0:
                            found.append((pp,ty,vv))
                    if len(found)==0: nopartner+=1; continue
                    if len(found)>1: multi+=1
                    pp,ty,vv=found[0]
                    ppkind[(ty, "pp<k" if pp<k else "pp>k")]+=1
                    Pw=set(Pset(u,w,k-1)); Pv=set(Pset(u,vv,k-1)); wk=val(w,k)
                    if Pv==Pw|{wk}: psetrel["P(v)=P(w)+{w(k)}"]+=1
                    else: psetrel[f"other diff={tuple(sorted(Pv-Pw))}/{tuple(sorted(Pw-Pv))}"]+=1
    print(f"moved instances DIAG+REC!=0 = {tot}")
    print(f"NO weight-matched partner reflection : {nopartner}")
    print(f"MULTIPLE partners : {multi}")
    print("partner (type, pos) histogram:")
    for kk,c in sorted(ppkind.items(),key=lambda x:-x[1]): print(f"    {kk}: {c}")
    print("P-set relation histogram:")
    for kk,c in sorted(psetrel.items(),key=lambda x:-x[1])[:10]: print(f"    {kk}: {c}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
