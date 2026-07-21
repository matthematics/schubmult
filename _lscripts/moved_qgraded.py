"""MOVED case: q-graded vanishing. Group every atomic term by its q-monomial factor.
Claim: within each fixed q-monomial, the factorial-elem polynomials sum to 0.
This would reduce the quantum vanishing to (shifted) CLASSICAL Pieri identities.
Report: instances where SOME q-group is nonzero (would refute the graded claim).
Also report distribution of #distinct q-monomials per instance.
Run: conda activate schubmult_312 && python _lscripts/moved_qgraded.py 5
"""
import sys
from collections import Counter, defaultdict
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, expand_func, prod, sympify
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
def c_parts(u,w,m,a,Rm):
    """Return (qmonomial_expr, poly_expr) so that c_term = qmon*poly, or None if zero."""
    if m<0: return None
    w=Permutation(w)
    if w not in Rm: return None
    qm=Rm[w]; Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    if deg<0 or deg>nv: return None
    if nv==0:
        return (qm, S.One) if deg==0 else None
    yv=[y[v] for v in Pm]; nc=nv+1-deg; zv=[z[i] for i in range(1,nc+1)]
    poly=expand_func(FactorialElemSym(deg,nv,yv,zv))
    return (qm, poly)

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; badgroup=0; ngroups=Counter()
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
                    groups=defaultdict(lambda: S.Zero)
                    any_term=False
                    # DIAG
                    cp=c_parts(u,w,k-1,p-1,Rkm1)
                    if cp: groups[str(expand(cp[0]))]+= (y[val(w,k)]-z[k-p+1])*cp[1]*cp[0]/cp[0]*0 + (y[val(w,k)]-z[k-p+1])*cp[1]; any_term=True
                    # careful: store poly only, keyed by qmon
                    # redo cleanly below
                    groups=defaultdict(lambda: S.Zero); qmon_expr={}
                    def add(qmon, poly):
                        key=str(expand(qmon)); qmon_expr[key]=qmon; groups[key]+=poly
                    cp=c_parts(u,w,k-1,p-1,Rkm1)
                    if cp: add(cp[0], (y[val(w,k)]-z[k-p+1])*cp[1]); any_term=True
                    cp=c_parts(u,w,k-1,p,Rkm1)
                    if cp: add(cp[0], cp[1]); any_term=True
                    cp=c_parts(u,w,k-2,p-2,Rkm2)
                    if cp: add(q_gs[k-1]*cp[0], cp[1]); any_term=True
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: extra=S.One; sgn=(1 if pp>k else -1)
                        elif dl==1-2*(hi-lo): extra=q_ab(lo,hi); sgn=(1 if pp>k else -1)
                        else: continue
                        cp=c_parts(u,vv,k-1,p-1,Rkm1)
                        if cp: add(extra*cp[0], sgn*cp[1]); any_term=True
                    if not any_term: continue
                    tot+=1
                    nz=0
                    for key,poly in groups.items():
                        if expand(poly)!=0: nz+=1; badgroup_local=True
                    ngroups[len([1 for key,poly in groups.items() if expand(poly)!=0 or True])]+=1
                    if nz>0: badgroup+=1
    print(f"moved instances = {tot}")
    print(f"instances with a NONZERO q-group (refutes graded vanishing) = {badgroup}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
