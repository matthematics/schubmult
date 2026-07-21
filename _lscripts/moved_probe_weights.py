"""Probe the up (r=B) vs down (r=c_i*) partner in the B>k world:
compare the two-piece cycle structure, descents, and q-weights, to find WHY they match
and what pins i*. Print a handful of fully-worked instances.

pi=u^{-1}w, C=cycle(k), orbit from B: [c0=B,...,c_{L-1}], c_j=k.
up partner v_B: pieces K0=(k;c1..c_{j-1}), A0=(B;c_{j+1}..c_{L-1}).
down partner v_{c_i*}: pieces K=(k;c_{i*+1}..c_{j-1}), A=(B;c1..c_i*,c_{j+1}..c_{L-1}).
For each we show: which c_i (1<=i<=j-1) are length-valid (u->_{k-1}), their kappa and q-weight.
Descents: desc_u(cycle) = positions a in cycle (a<=k-1 side) that are descents (u(a)>u(prev)).
"""
import sys
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
def orbit(perm,start,N):
    c=[start]; cur=perm[start-1]
    while cur!=start: c.append(cur); cur=perm[cur-1]
    return c

def run(n, want=12):
    N=2*n; perms=list(Permutation.all_permutations(n)); shown=0
    for u in perms:
        u=Permutation(u)
        for k in range(2,n):
            Rkm1=enumerate_pieri(u,k,N) and enumerate_pieri(u,k-1,N)
            Rkm1=enumerate_pieri(u,k-1,N); Rk=enumerate_pieri(u,k,N)
            cset=set(Rkm1)
            for v0 in list(Rkm1):
                for pp in range(1,N+1):
                    if pp!=k: cset.add(v0.swap(k-1,pp-1))
            for w in cset:
                w=Permutation(w)
                if val(u,k)==val(w,k): continue
                if w in Rk: continue
                pi=(~u)*w; C=orbit(pi,k,N); B=max(C)
                if B<=k: continue
                orb=orbit(pi,B,N); L=len(orb); j=orb.index(k)
                valids=[]
                for i in range(0,j):  # i=0 up (r=B), i=1..j-1 down
                    r=orb[i]; vv=w.swap(k-1,r-1); dl=w.inv-vv.inv
                    lo,hi=min(r,k),max(r,k)
                    if dl==1: kap=S.One
                    elif dl==1-2*(hi-lo): kap=q_ab(lo,hi)
                    else: continue
                    if vv not in Rkm1: continue
                    valids.append((i,r,expand(kap),expand(Rkm1[vv]),expand(kap*Rkm1[vv])))
                if len(valids)<2: continue
                # show
                print(f"u={tuple(u)} w={tuple(w)} k={k} B={B} j={j}")
                print(f"  C orbit from B: {orb}")
                for (i,r,kap,qq,wt) in valids:
                    tag="UP(B)" if i==0 else f"DOWN c_{i}"
                    print(f"    i={i} r={r} [{tag}] kappa={kap} q_{{k-1}}={qq} weight={wt}")
                shown+=1
                if shown>=want: return

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
