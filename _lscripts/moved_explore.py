"""Discover the explicit involution rule for Case (ii) of lemma:vanishing / prop:qreflection.

Setup: u NOT->_k w, u(k)!=w(k). Build the four-term recursion F's atomic terms, group by
exact q-monomial M. For Case (ii) groups (no DIAG/REC/straight), we must pair the Monk
neighbors (+ QUANT) into opposite-sign equal-E-signature pairs.

For each pure-Monk pair (down r<k, up s>k) we record the RELATIONSHIP to pi=u^{-1}w:
  - C = cycle of pi containing k; its sorted elements; top B=max(C)
  - the >=k indices of C (call them k=i0<i1<...<it=B)
  - where r, s sit relative to C and these indices
  - candidate rule guesses (s==B? r,s adjacency in C? etc.)

For QUANT pairs we record a* (the down partner) and its cycle role.

Run: conda activate schubmult_312 && python _lscripts/moved_explore.py 5
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
    if nv==0:
        poly=S.One if deg==0 else S.Zero
    else:
        yv=[y[v] for v in Pm]; nc=nv+1-deg; zv=[z[i] for i in range(1,nc+1)]
        poly=expand_func(FactorialElemSym(deg,nv,yv,zv))
    return (expand(qm), poly, (deg,nv,Pm))

def cycle_of(perm,elt,N):
    c=[elt]; cur=perm[elt-1]
    while cur!=elt: c.append(cur); cur=perm[cur-1]
    return c  # in orbit order starting at elt

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    rule_s_is_B=Counter(); rule_r_in_C=Counter(); rule_s_in_C=Counter()
    downs_ups_hist=Counter()
    s_position=Counter()  # where does s sit among >=k indices of C
    r_position=Counter()
    examples=[]
    quant_astar=Counter()
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
                C=sorted(cycle_of(pi,k,N))
                B=max(C); geqk=[x for x in C if x>=k]
                for p in range(1,k+1):
                    terms=[]  # (M, signedpoly, Esig, kind, rpos)
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
                        if "DIAG" in kinds or "REC" in kinds: continue
                        # Case (ii). partition by Esig
                        bysig=defaultdict(list)
                        for poly,Esig,kind,rpos in g: bysig[Esig].append((poly,kind,rpos))
                        for Esig,lst in bysig.items():
                            if len(lst)!=2: continue
                            (p1,k1,r1),(p2,k2,r2)=lst
                            kinds2={k1,k2}
                            if kinds2=={"MONK"}:
                                # down/up
                                rr = r1 if r1<k else r2
                                ss = r1 if r1>k else r2
                                downs_ups_hist[(rr<k, ss>k)]+=1
                                rule_s_is_B[ss==B]+=1
                                rule_r_in_C[rr in C]+=1
                                rule_s_in_C[ss in C]+=1
                                if geqk:
                                    s_position[ (geqk.index(ss) if ss in geqk else -1, len(geqk)) ]+=1
                                    r_position[ (C.index(rr) if rr in C else -1) ]+=1
                                if len(examples)<25:
                                    examples.append((tuple(u),tuple(w),k,p,rr,ss,C,B,geqk))
                            elif "QUANT" in kinds2:
                                astar = r1 if k1=="MONK" else r2
                                quant_astar[(astar<k, astar in C, astar==B)]+=1
    print("down/up hist (r<k, s>k):", dict(downs_ups_hist))
    print("rule s==B :", dict(rule_s_is_B))
    print("rule r in C :", dict(rule_r_in_C))
    print("rule s in C :", dict(rule_s_in_C))
    print("s position among >=k indices of C (idx,len):", dict(s_position))
    print("r position in C (index in sorted C, -1 if not):", dict(r_position))
    print("QUANT partner a* (a*<k, a* in C, a*==B):", dict(quant_astar))
    print("--- examples (u,w,k,p,r,s,C,B,geqk) ---")
    for e in examples: print(e)

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
