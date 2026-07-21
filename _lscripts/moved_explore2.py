"""Detailed orbit-order structure of the cycle C containing k, and ALL contributing
position-k reflections, classified. Goal: pin the exact involution r<->B.

For each moved instance (u NOT->_k w, u(k)!=w(k)), pi=u^{-1}w, C=cycle(k).
Write C in pi-orbit order starting at B=max(C): [B=c0, c1, ..., c_{L-1}], pi(c_i)=c_{i+1}.
Let j be s.t. c_j = k.  a0 = pi(k) = c_{j+1 mod L} = u^{-1}(w(k)).
For each p and each contributing Monk neighbor r (wt_{rk} lup w, u->_{k-1} wt_{rk}):
  record (orbit-index i with c_i=r [or 'notC'], sign(r-k), kappa is drop?, weight-class label,
          RECURR or CANCEL, E-sig).
Aggregate: for CANCEL down-neighbors, what is orbit-index i of r, and relationship to j?

Run: conda activate schubmult_312 && python _lscripts/moved_explore2.py 5
"""
import sys
from collections import Counter, defaultdict
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
def orbit_from(perm,start,N):
    c=[start]; cur=perm[start-1]
    while cur!=start: c.append(cur); cur=perm[cur-1]
    return c

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    cancel_down_i_minus_j=Counter()   # orbit index i of down-r minus j (k's index)
    cancel_down_isucc=Counter()       # is down-r == c_{j+1}? c_{j-1}? etc
    recurr_a0_check=Counter()
    twoidx_ge_k=Counter()
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
                orb=orbit_from(pi,max(orbit_from(pi,k,N)),N)  # start at B=max of C
                L=len(orb); B=orb[0]; j=orb.index(k)
                geqk=[x for x in orb if x>=k]
                twoidx_ge_k[len(geqk)]+=1
                for p in range(1,k+1):
                    terms=[]
                    s=sig(u,w,k-1,p-1,Rkm1)
                    if s: terms.append((s[0], expand((y[val(w,k)]-z[k-p+1])*s[1]), None, "DIAG", None))
                    s=sig(u,w,k-1,p,Rkm1)
                    if s: terms.append((s[0], s[1], s[2], "REC", None))
                    s=sig(u,w,k-2,p-2,Rkm2)
                    if s: terms.append((expand(q_gs[k-1]*s[0]), s[1], s[2], "QUANT", None))
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: kap=S.One; sgn=1 if pp>k else -1
                        elif dl==1-2*(hi-lo): kap=q_ab(lo,hi); sgn=1 if pp>k else -1
                        else: continue
                        s=sig(u,vv,k-1,p-1,Rkm1)
                        if s: terms.append((expand(kap*s[0]), expand(sgn*s[1]), s[2], "MONK", pp))
                    groups=defaultdict(list)
                    for M,poly,Esig,kind,rpos in terms: groups[M].append((poly,Esig,kind,rpos))
                    for M,g in groups.items():
                        kinds=[x[2] for x in g]
                        if "DIAG" in kinds or "REC" in kinds:
                            # RECURR: single monk r=a0=c_{j+1}
                            for poly,Esig,kind,rpos in g:
                                if kind=="MONK":
                                    a0=pi[k-1]
                                    recurr_a0_check[(rpos==a0, rpos==orb[(j+1)%L])]+=1
                            continue
                        bysig=defaultdict(list)
                        for poly,Esig,kind,rpos in g: bysig[Esig].append((poly,kind,rpos))
                        for Esig,lst in bysig.items():
                            if len(lst)!=2: continue
                            (p1,k1,r1),(p2,k2,r2)=lst
                            if {k1,k2}=={"MONK"}:
                                rr=r1 if r1<k else r2
                                i_r=orb.index(rr) if rr in orb else -99
                                cancel_down_i_minus_j[i_r-j]+=1
                                # relationship guesses
                                cancel_down_isucc[( rr==orb[(j-1)%L], rr==orb[(j+1)%L], i_r<j, i_r>j )]+=1
    print("size of {indices>=k} in cycle C (should be 2 whenever a CANCEL pair exists):", dict(twoidx_ge_k))
    print("CANCEL down-r orbit-index minus j (k's index):", dict(cancel_down_i_minus_j))
    print("CANCEL down-r relationship (is c_{j-1}, is c_{j+1}, i<j, i>j):", dict(cancel_down_isucc))
    print("RECURR monk r (==a0=pi(k), ==c_{j+1}):", dict(recurr_a0_check))

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
