"""MOVED-B0 with QUANT!=0: find the exact partner mechanism.

Case 3 moved: u NOT->_k w, u(k)!=w(k), and QUANT = q_{k-1} c_{k-2}(p-2;u,w) != 0
(so u ->_{k-2} w). We must show DIAG+UP+DOWN+REC+QUANT = 0.

Split the reflection sum into raise (q-free) and drop (q-carrying) parts.
Test the hypothesis:
  (H1) DIAG + UP_raise + DOWN_raise + REC = 0     [classical remainder vanishes]
  (H2) UP_drop + DOWN_drop + QUANT = 0            [quantum part vanishes]
If both hold with 0 failures, moved-B0 = classical-vanishing (H1) + quantum-cancel (H2).

BUT c_{k-1} terms carry their own q, so H1/H2 may not split cleanly. So ALSO test a
finer split by exact q-monomial: group every term by its q-monomial and check each
q-graded piece sums to 0.

Also, for the QUANT partner: identify v=wt (t=t_{ak} a<k, or t=t_{kb} b>k) such that
  c_{k-1}(p-1;u,v) has q-weight q_{k-1}(u,v) with q_{k-1}*q_{k-2}(u,w) == weight,
i.e. the term that carries the SAME q-monomial as QUANT. Report how many such partners.

Run: conda activate schubmult_312 && python _lscripts/moved_quant.py 5
"""
import sys
from collections import Counter, defaultdict
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
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; h1_fail=0; h2_fail=0; graded_fail=0
    partner_hist=Counter()
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
                if val(u,k)==val(w,k): continue     # moved only
                if w in Rk: continue                 # u NOT ->_k w
                for p in range(1,k+1):
                    if not c_nz(u,w,k-2,p-2,Rkm2): continue  # QUANT != 0 only
                    tot+=1
                    diag=(y[val(w,k)]-z[k-p+1])*c_term(u,w,k-1,p-1,Rkm1)
                    rec=c_term(u,w,k-1,p,Rkm1)
                    quant=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    up_raise=S.Zero; up_drop=S.Zero; dn_raise=S.Zero; dn_drop=S.Zero
                    partners=0
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        cval=c_term(u,vv,k-1,p-1,Rkm1)
                        if dl==1:
                            if pp>k: up_raise+=cval
                            else: dn_raise-=cval
                        elif dl==1-2*(hi-lo):
                            if pp>k: up_drop+=q_ab(lo,hi)*cval
                            else: dn_drop-=q_ab(lo,hi)*cval
                    # H1: classical remainder
                    h1=expand(diag+up_raise+dn_raise+rec)
                    # H2: quantum part
                    h2=expand(up_drop+dn_drop+quant)
                    if h1!=0: h1_fail+=1
                    if h2!=0: h2_fail+=1
                    # graded by q-monomial
                    full=expand(diag+up_raise+up_drop+dn_raise+dn_drop+rec+quant)
                    if full!=0: graded_fail+=1
    print(f"moved-B0 QUANT!=0 instances = {tot}")
    print(f"(H1) DIAG+UP_raise+DOWN_raise+REC != 0 : {h1_fail}")
    print(f"(H2) UP_drop+DOWN_drop+QUANT != 0      : {h2_fail}")
    print(f"full RHS != 0                          : {graded_fail}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
