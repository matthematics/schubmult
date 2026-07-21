"""Inspect the q-monomial group structure of the moved-case RHS F on sample instances.
For each instance print, per q-monomial M, the member terms (type + which perm + (deg,nv,Pset))
so we can SEE that each group is a classical Pieri recursion RHS for a residual config.
Run: conda activate schubmult_312 && python _lscripts/moved_inspect.py 5
"""
import sys
from collections import defaultdict
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
def cinfo(u,w,m,a,Rm):
    """Return (qmon, poly, (deg,nv,Pset)) or None."""
    if m<0: return None
    w=Permutation(w)
    if w not in Rm: return None
    qm=Rm[w]; Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    if deg<0 or deg>nv: return None
    if nv==0:
        return (qm,S.One,(deg,nv,tuple(Pm))) if deg==0 else None
    yv=[y[v] for v in Pm]; nc=nv+1-deg; zv=[z[i] for i in range(1,nc+1)]
    return (qm,expand_func(FactorialElemSym(deg,nv,yv,zv)),(deg,nv,tuple(Pm)))

def run(n, want=6):
    N=2*n; perms=list(Permutation.all_permutations(n)); shown=0
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
                    groups=defaultdict(list); qexpr={}
                    def add(label, qmon, poly, info):
                        kk=str(expand(qmon)); qexpr[kk]=qmon; groups[kk].append((label,poly,info))
                    ci=cinfo(u,w,k-1,p-1,Rkm1)
                    if ci: add("DIAG", ci[0], (y[val(w,k)]-z[k-p+1])*ci[1], ci[2])
                    ci=cinfo(u,w,k-1,p,Rkm1)
                    if ci: add("REC", ci[0], ci[1], ci[2])
                    ci=cinfo(u,w,k-2,p-2,Rkm2)
                    if ci: add("QUANT", q_gs[k-1]*ci[0], ci[1], ci[2])
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: extra=S.One; sgn=(1 if pp>k else -1); ty=f"RAISE(pp={pp})"
                        elif dl==1-2*(hi-lo): extra=q_ab(lo,hi); sgn=(1 if pp>k else -1); ty=f"DROP(pp={pp})"
                        else: continue
                        ci=cinfo(u,vv,k-1,p-1,Rkm1)
                        if ci: add(ty, extra*ci[0], sgn*ci[1], ci[2])
                    # only show instances with >1 q-group and a QUANT present
                    if len(groups)<2: continue
                    if not any("QUANT" in [t[0] for t in g] for g in groups.values()): continue
                    shown+=1
                    print(f"\n=== u={tuple(u)} w={tuple(w)} k={k} p={p} ===")
                    for kk,mem in groups.items():
                        tot=expand(sum(m[1] for m in mem))
                        print(f"  q-monomial {kk}:  sum={'0' if tot==0 else tot}")
                        for lab,poly,info in mem:
                            print(f"      {lab:14s} (deg,nv,P)={info}")
                    if shown>=want: return
    print("done")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
