"""MOVED-B0 QUANT!=0: verify P-set/n/q-weight matching between partner v and w,
and identify the exact structural relation between v and w (since v t_{k-1,k} != w here).

partner v = w t_{pp,k}, pp<k, u->_{k-1}v, NOT u->_{k-2}v.
Check:  P_{k-1}(u,v) == P_{k-2}(u,w)  ;  n_{k-1}(u,v) == n_{k-2}(u,w)+1.
Also: apply pieriknotkn1 to v at level k-1: unique b>k-1 with v t_{k-1,b} <| v and
u->_{k-2}(v t_{k-1,b}). Record that b and w' = v t_{k-1,b}; is w' == w? if not, what is it?
Record the cycle-top-through-k of u^{-1}v.
Run: conda activate schubmult_312 && python _lscripts/moved_match.py 5
"""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet

q_gs = GeneratingSet("q")
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
def nk(u,w,m): return sum(1 for i in range(1,m+1) if val(u,i)!=val(w,i))
def c_nz(u,w,m,a,Rm):
    if m<0: return False
    w=Permutation(w)
    if w not in Rm: return False
    Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    return 0<=deg<=nv
def cyc_top(u,v,idx):
    perm=(~u)*v; cyc=[idx]; cur=perm[idx-1]
    while cur!=idx: cyc.append(cur); cur=perm[cur-1]
    return max(cyc) if len(cyc)>1 else None

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; pfail=0; nfail=0; btop=Counter(); wprime_is_w=0; wprime_not=0
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
                    if not c_nz(u,w,k-2,p-2,Rkm2): continue
                    # find intrinsic partner
                    vv=None; pp_star=None
                    for pp in range(1,k):
                        cand=w.swap(k-1,pp-1); dl=w.inv-cand.inv
                        lo,hi=pp,k
                        if dl==1 or dl==1-2*(hi-lo):
                            if c_nz(u,cand,k-1,p-1,Rkm1) and not reachable(u,cand,k-2,N):
                                vv=cand; pp_star=pp; break
                    if vv is None: continue
                    tot+=1
                    P1=Pset(u,vv,k-1); P2=Pset(u,w,k-2)
                    if P1!=P2: pfail+=1
                    if nk(u,vv,k-1)!=nk(u,w,k-2)+1: nfail+=1
                    # pieriknotkn1 on v at level k-1
                    top=cyc_top(u,vv,k-1)
                    btop[top]+=1
                    wprime=vv.swap(k-2, (top-1) if top else k-1)  # v t_{k-1, top}
                    if top is not None and Permutation(wprime)==w: wprime_is_w+=1
                    else: wprime_not+=1
    print(f"moved-B0 QUANT!=0 (partner found) = {tot}")
    print(f"P-set mismatch P_{{k-1}}(u,v)!=P_{{k-2}}(u,w) : {pfail}")
    print(f"n mismatch                                  : {nfail}")
    print(f"cycle-top-through-(k-1) of u^-1 v histogram (relative to k): ")
    krel=Counter()
    # recompute relative not easy here; just print raw
    for t,c in sorted(btop.items(), key=lambda x:-x[1])[:10]: print(f"    top={t}: {c}")
    print(f"v t_{{k-1,top}} == w : {wprime_is_w} ; != w : {wprime_not}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
