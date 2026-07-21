"""Verify the load-bearing claim for the moved-case write-up:
In the moved case (u(k)!=w(k), u NOT->_k w), the combination DIAG+REC equals a SINGLE
length-increasing Monk coefficient c_{k-1}(p-1; u, w') where w' = w t_{k,a0},
a0 = u^{-1}(w(k)) (the unique position where value w(k) agrees with u), via the
single-variable recurrence (elemrecursion) adding y_{w(k)}; and:
  - the step w->w' is length-INCREASING (raise),
  - q_{k-1}(u,w') == q_{k-1}(u,w)  (same q-weight M0),
  - so DIAG+REC + [that Monk summand with its sign(a0-k)] == 0 inside F_{M0}.
Report failures of each sub-claim. Also confirm: after removing DIAG,REC and this w',
the LEFTOVER of the whole RHS F still == 0 (so the rest is a pure reflection/quant identity).
Run: conda activate schubmult_312 && python _lscripts/moved_diagrec.py 5
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
    tot=0; fexist=0; fincr=0; fweight=0; feq=0; fleft=0; a0_lt=0; a0_gt=0
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
                    if expand(diag+rec)==0:
                        continue  # nothing to pair (both zero)
                    tot+=1
                    # a0 = position where value w(k) agrees with u: u(a0)=w(k)
                    wk=val(w,k)
                    a0=None
                    for i in range(1,N+2):
                        if val(u,i)==wk: a0=i; break
                    if a0 is None or a0==k:
                        fexist+=1; continue
                    wp=w.swap(k-1,a0-1)  # w t_{k,a0}
                    lo,hi=min(a0,k),max(a0,k)
                    incr = (wp.inv==w.inv+1)
                    isdrop = (wp.inv==w.inv-2*(hi-lo)+1)
                    if a0<k: a0_lt+=1
                    else: a0_gt+=1
                    if not (incr or isdrop):
                        fincr+=1; continue  # not a valid Monk step at all
                    kappa = S.One if incr else q_ab(lo,hi)
                    # q-weight of the Monk SUMMAND = kappa * q_{k-1}(u,wp); should equal q_{k-1}(u,w)
                    if wp not in Rkm1:
                        fexist+=1; continue
                    if expand(kappa*Rkm1[wp]-Rkm1[w])!=0: fweight+=1
                    # DIAG+REC == (signed, weighted) Monk SUMMAND at a0?
                    monk=kappa*c_term(u,wp,k-1,p-1,Rkm1)   # bare summand magnitude
                    if expand(diag+rec-monk)!=0: feq+=1
                    # leftover: full RHS minus (diag+rec) minus [sign(a0-k)*monk] == 0?
                    # build full reflection sum + quant
                    refl=S.Zero
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); ddl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if ddl==1: tv=(1 if pp>k else -1)*c_term(u,vv,k-1,p-1,Rkm1)
                        elif ddl==1-2*(hi-lo): tv=(1 if pp>k else -1)*q_ab(lo,hi)*c_term(u,vv,k-1,p-1,Rkm1)
                        else: continue
                        refl+=tv
                    quant=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    fullF=diag+rec+refl+quant
                    # the a0 monk summand carries sign(a0-k) and weight kappa:
                    signed_monk=(1 if a0>k else -1)*monk
                    leftover=fullF-(diag+rec)-signed_monk
                    if expand(leftover)!=0: fleft+=1
    print(f"moved instances with DIAG+REC!=0 = {tot}")
    print(f"a0 missing/==k                 : {fexist}")
    print(f"step w->w t_(k,a0) neither raise nor drop : {fincr}")
    print(f"weight(summand) != q_(k-1)(u,w)  : {fweight}")
    print(f"DIAG+REC != signed-weighted Monk summand : {feq}")
    print(f"leftover (F - DIAG-REC - signed monk) != 0 : {fleft}")
    print(f"a0<k count={a0_lt}, a0>k count={a0_gt}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
