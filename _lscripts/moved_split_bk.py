"""MOVED case (u NOT->_k w, u(k)!=w(k)): test the b>k vs a<k partition.
Group UP  = DIAG + REC + sum_{b>k} reflection(t_{kb})          -> ?= 0
Group DN  = sum_{a<k} reflection(t_{ak}) + QUANT               -> ?= 0
Also test full RHS == 0 (sanity). Enumerate ALL moved non-reachable w (not just QUANT!=0).
Run: conda activate schubmult_312 && python _lscripts/moved_split_bk.py 5
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

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; upf=0; dnf=0; fullf=0
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
                    quant=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    up=S.Zero; dn=S.Zero
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: tval=(1 if pp>k else -1)*c_term(u,vv,k-1,p-1,Rkm1)
                        elif dl==1-2*(hi-lo): tval=(1 if pp>k else -1)*q_ab(lo,hi)*c_term(u,vv,k-1,p-1,Rkm1)
                        else: continue
                        if pp>k: up+=tval
                        else: dn+=tval
                    if diag==0 and rec==0 and quant==0 and up==0 and dn==0: continue
                    tot+=1
                    if expand(diag+rec+up)!=0: upf+=1
                    if expand(dn+quant)!=0: dnf+=1
                    if expand(diag+rec+up+dn+quant)!=0: fullf+=1
    print(f"moved non-reachable instances (some term nonzero) = {tot}")
    print(f"UP group (DIAG+REC+sum_b>k) != 0 : {upf}")
    print(f"DN group (sum_a<k + QUANT)   != 0 : {dnf}")
    print(f"FULL RHS != 0                     : {fullf}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
