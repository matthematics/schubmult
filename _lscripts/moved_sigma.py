"""Reduce Case (ii) B>k to a pure q-weight identity Sigma=0 and probe its structure.

B>k world: pi=u^{-1}w, C=cycle(k) written in pi-orbit order from B=max(C): [c0=B,...,c_{L-1}],
c_j=k. Contributing reflections t_{rk} with u->_{k-1} wt_{rk} are among r in {B=c0, c1,...,c_{j-1}}.
We PROVED P_{k-1} and n_{k-1} are identical for ALL these neighbors, so a common factorial-elem
polynomial E factors out and F = E * Sigma with
    Sigma = sum over valid r  sign(r-k) * kappa_{rk} * q_{k-1}(u, wt_{rk}).
This script (length-valid reflections only; p-independent):
 (S1) valid r are all in {c0=B, c1,...,c_{j-1}} (i.e. i in 0..j-1);
 (S2) P_{k-1}(u,wt_{rk}) == P_{k-1}(u,w) and n_{k-1} match, for every valid r;
 (S3) exactly ONE valid up (r=B) and exactly ONE valid down (r=c_{i*}, 1<=i*<=j-1) when nonempty;
 (S4) their weights are equal:  kappa_{Bk} q_{k-1}(u,wt_{Bk}) == kappa_{c_i*,k} q_{k-1}(u,wt_{c_i*,k});
 (S5) Sigma == 0 as a q-polynomial.
Also records orbit-index i* of the down partner and the length-type (drop?) of each.

Run: conda activate schubmult_312 && python _lscripts/moved_sigma.py 5
"""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet

q_gs=GeneratingSet("q")
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
def Pset(u,w,m): return frozenset(val(u,i) for i in range(1,m+1) if val(u,i)==val(w,i))
def orbit(perm,start,N):
    c=[start]; cur=perm[start-1]
    while cur!=start: c.append(cur); cur=perm[cur-1]
    return c

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    s1fail=0; s2fail=0; s3fail=0; s4fail=0; s5fail=0
    up_cnt=Counter(); dn_cnt=Counter(); istar=Counter(); nonempty=0
    for u in perms:
        u=Permutation(u)
        for k in range(2,n):
            Rk=enumerate_pieri(u,k,N); Rkm1=enumerate_pieri(u,k-1,N)
            cset=set(Rkm1)
            for v0 in list(Rkm1):
                for pp in range(1,N+1):
                    if pp!=k: cset.add(v0.swap(k-1,pp-1))
            cset|=set(Rk)
            for w in cset:
                w=Permutation(w)
                if val(u,k)==val(w,k): continue
                if w in Rk: continue
                pi=(~u)*w
                C=orbit(pi,k,N); B=max(C)
                if B<=k: continue   # only B>k world
                orb=orbit(pi,B,N); L=len(orb); j=orb.index(k)
                Pw=Pset(u,w,k-1); nw=(k-1)-len(Pw)
                allowed_r=set(orb[0:j])  # c0..c_{j-1}
                ups=[]; downs=[]
                Sigma=S.Zero
                for pp in range(1,N+1):
                    if pp==k: continue
                    vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                    lo,hi=min(pp,k),max(pp,k)
                    if dl==1: kap=S.One
                    elif dl==1-2*(hi-lo): kap=q_ab(lo,hi)
                    else: continue
                    if vv not in Rkm1: continue     # u->_{k-1} vv
                    # S1
                    if pp not in allowed_r: s1fail+=1
                    # S2
                    if Pset(u,vv,k-1)!=Pw or (k-1)-len(Pset(u,vv,k-1))!=nw: s2fail+=1
                    wt=expand(kap*Rkm1[vv]); sgn=1 if pp>k else -1
                    Sigma+=sgn*wt
                    if pp>k: ups.append((pp,wt))
                    else: downs.append((pp,wt,orb.index(pp) if pp in orb else -1))
                up_cnt[len(ups)]+=1; dn_cnt[len(downs)]+=1
                if ups or downs:
                    nonempty+=1
                    if not(len(ups)==1 and len(downs)==1): s3fail+=1
                    else:
                        if expand(ups[0][1]-downs[0][1])!=0: s4fail+=1
                        istar[downs[0][2]-j]+=1   # i* - j (negative)
                if expand(Sigma)!=0: s5fail+=1
    print(f"B>k nonempty instances: {nonempty}")
    print(f"S1 (valid r not in c0..c_{{j-1}}): {s1fail}")
    print(f"S2 (P/n mismatch): {s2fail}")
    print(f"S3 (not exactly 1 up +1 down): {s3fail} | up-count {dict(up_cnt)} down-count {dict(dn_cnt)}")
    print(f"S4 (up weight != down weight): {s4fail}")
    print(f"S5 (Sigma != 0): {s5fail}")
    print(f"down partner orbit-index i*-j : {dict(istar)}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
