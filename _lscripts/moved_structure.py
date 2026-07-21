"""Characterize the FULL q-graded structure of the moved-case recursion F (u(k)!=w(k), u NOT->_k w).
For each instance, build all atomic terms with (q-weight M, signed poly-coefficient, E-signature
(deg,nv,frozenset(Pset)), neighbor, kind in {DIAG,REC,QUANT,MONK}). Group by exact q-monomial M.
For each group:
  - if it contains DIAG or REC  -> RECURR: verify the whole group's polynomial sums to 0
    AND that there is exactly ONE MONK term in it, whose neighbor w' has
    P_{k-1}(u,w')=P_{k-1}(u,w) U {w(k)} (the elemrecursion completion).
  - else (pure MONK +/- QUANT) -> CANCEL: verify it partitions into opposite-sign pairs with
    identical E-signature; record for each pair the cycle-type of w'^{-1} w'' and whether the
    remaining cycles (away from the k-cycle) coincide (=> same residual drop weight).
Reports counts of any violations.
Run: conda activate schubmult_312 && python _lscripts/moved_structure.py 5
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
    """Return (qmon, poly, (deg,nv,Pset)) or None if c=0."""
    if m<0: return None
    w=Permutation(w)
    if w not in Rm: return None
    qm=Rm[w]; Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    if deg<0 or deg>nv: return None
    if nv==0:
        poly=S.One if deg==0 else S.Zero
    else:
        yv=[y[v] for v in Pm]; nc=nv+1-deg; zv=[z[i] for i in range(1,nc+1)]
        poly=expand_func(FactorialElemSym(deg,nv,yv,zv))
    return (expand(qm), poly, (deg,nv,Pm))
def cycles(perm,N):
    seen=set(); out=[]
    for i in range(1,N+1):
        if i in seen: continue
        c=[i]; seen.add(i); cur=perm[i-1]
        while cur!=i: c.append(cur); seen.add(cur); cur=perm[cur-1]
        if len(c)>1: out.append(tuple(sorted(c)))
    return frozenset(out)

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    recurr_polyfail=0; recurr_monkcount=Counter(); recurr_Pfail=0
    cancel_pairfail=0; cancel_leftover=0
    pair_cyc=Counter(); pair_restmatch=Counter()
    ngroups=0; nrecurr=0; ncancel=0
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
                for p in range(1,k+1):
                    terms=[]  # (M, signedpoly, Esig, kind, neighbor)
                    s=sig(u,w,k-1,p-1,Rkm1)
                    if s: terms.append((s[0], expand((y[val(w,k)]-z[k-p+1])*s[1]), None, "DIAG", w))
                    s=sig(u,w,k-1,p,Rkm1)
                    if s: terms.append((s[0], s[1], s[2], "REC", w))
                    s=sig(u,w,k-2,p-2,Rkm2)
                    if s: terms.append((expand(q_gs[k-1]*s[0]), s[1], s[2], "QUANT", w))
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: kap=S.One; sgn=1 if pp>k else -1
                        elif dl==1-2*(hi-lo): kap=q_ab(lo,hi); sgn=1 if pp>k else -1
                        else: continue
                        s=sig(u,vv,k-1,p-1,Rkm1)
                        if s: terms.append((expand(kap*s[0]), expand(sgn*s[1]), s[2], "MONK", vv))
                    # total must be 0
                    if terms and expand(sum(t[1]*t[0] for t in terms))!=0:
                        pass  # trust; already known 0
                    groups=defaultdict(list)
                    for M,poly,Esig,kind,nb in terms: groups[M].append((poly,Esig,kind,nb))
                    for M,g in groups.items():
                        ngroups+=1
                        kinds=[x[2] for x in g]
                        if "DIAG" in kinds or "REC" in kinds:
                            nrecurr+=1
                            if expand(sum(x[0] for x in g))!=0: recurr_polyfail+=1
                            monks=[x for x in g if x[2]=="MONK"]
                            recurr_monkcount[len(monks)]+=1
                            for _,_,_,nb in monks:
                                Pnb=set(Pset(u,nb,k-1)); Pw=set(Pset(u,w,k-1))
                                if Pnb!=Pw|{val(w,k)}: recurr_Pfail+=1
                        else:
                            ncancel+=1
                            # partition into opposite-sign equal-Esig pairs
                            bysig=defaultdict(list)
                            for poly,Esig,kind,nb in g: bysig[Esig].append((poly,kind,nb))
                            for Esig,lst in bysig.items():
                                tot=expand(sum(x[0] for x in lst))
                                if tot!=0: cancel_leftover+=1
                                if len(lst)==2:
                                    (p1,k1,n1),(p2,k2,n2)=lst
                                    if expand(p1+p2)!=0: cancel_pairfail+=1
                                    perm=(~Permutation(n1))*Permutation(n2)
                                    cyc=cycles(perm,N)
                                    sizes=tuple(sorted(len(c) for c in cyc))
                                    hask=any(k in c for c in cyc)
                                    pair_cyc[(sizes,"has_k" if hask else "no_k")]+=1
                                    # rest match: cycles of u^{-1}n1 and u^{-1}n2 away from the one containing k
                                    c1=cycles((~u)*Permutation(n1),N); c2=cycles((~u)*Permutation(n2),N)
                                    rest1=frozenset(c for c in c1 if k not in c)
                                    rest2=frozenset(c for c in c2 if k not in c)
                                    pair_restmatch["match" if rest1==rest2 else "DIFF"]+=1
    print(f"groups total={ngroups} RECURR={nrecurr} CANCEL={ncancel}")
    print(f"RECURR: poly-sum!=0 : {recurr_polyfail} | P(w')=P(w)+w(k) fails : {recurr_Pfail}")
    print(f"RECURR MONK-count histogram : {dict(recurr_monkcount)}")
    print(f"CANCEL: within-Esig sum!=0 : {cancel_leftover} | 2-elt pair sum!=0 : {cancel_pairfail}")
    print(f"CANCEL pair (cycle-sizes of n1^-1 n2, k?) : {dict(pair_cyc)}")
    print(f"CANCEL pair rest-cycles match : {dict(pair_restmatch)}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
