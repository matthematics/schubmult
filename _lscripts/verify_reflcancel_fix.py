"""Verify the facts needed for the CORRECTED reflectioncancel Step 2-3 (no v* t_{k-1,k}=w).
For all reachable (u,w,m): n_m == m - |P_m|.
For reflectioncancel contributing instances (u(k)=w(k), u->_{k-2}w, contributing a*<k, v*=w t_{a*,k}):
  (A) cycle of pi=u^{-1}w containing k-1 has top k-1 (when w moves k-1); a* in that cycle
  (B) P_{k-1}(u,v*) == P_{k-2}(u,w)   and   n_{k-1}(u,v*) == n_{k-2}(u,w)+1
  (C) weight identity:  q_{a*,k} * q_{k-1}(u,v*)  ==  q_{k-1} * q_{k-2}(u,w)
  (D) t_{a*,k} is length-DECREASING (a drop)
Run: conda activate schubmult_312 && python _lscripts/verify_reflcancel_fix.py 5
"""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet

q_gs = GeneratingSet("q")
def q_ab(a,b): return prod([q_gs[s] for s in range(a,b)])   # a<=..<b
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
def Ppos(u,w,m): return [i for i in range(1,m+1) if val(u,i)==val(w,i)]
def nval(u,w,m,Rm):
    # n_m via the chain: length drop count = number of drops = derived from q-weight monomial degree/2? 
    # Instead compute from the coefficient convention m - |P_m|:
    return m-len(Pset(u,w,m))
def c_nonzero(u,w,m,a,Rm):
    if m<0: return False
    w=Permutation(w)
    if w not in Rm: return False
    Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    return 0<=deg<=nv
def cycle_containing(u,w,idx):
    perm=(~u)*w; cyc=[idx]; cur=perm[idx-1]
    while cur!=idx: cyc.append(cur); cur=perm[cur-1]
    return tuple(cyc) if len(cyc)>1 else None

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; A_fail=0; B_fail=0; C_fail=0; D_fail=0
    for u in perms:
        u=Permutation(u)
        for k in range(2,n):
            Rkm1=enumerate_pieri(u,k-1,N)
            Rkm2=enumerate_pieri(u,k-2,N) if k-2>=0 else {u:S.One}
            for w in list(Rkm2):
                w=Permutation(w)
                if val(u,k)!=val(w,k): continue
                for p in range(1,k+1):
                    if not c_nonzero(u,w,k-2,p-2,Rkm2): continue
                    astar=None; vv=None
                    for a in range(1,k):
                        vt=w.swap(k-1,a-1); dl=w.inv-vt.inv; lo,hi=a,k
                        if dl==1 or dl==1-2*(hi-lo):
                            if c_nonzero(u,vt,k-1,p-1,Rkm1):
                                astar=a; vv=vt; break
                    if astar is None: continue
                    tot+=1
                    # (D) t_{a*,k} drop or raise; kappa accordingly
                    dl=w.inv-vv.inv
                    is_drop = (dl==1-2*(k-astar))
                    is_raise = (dl==1)
                    kappa = q_ab(astar,k) if is_drop else S.One
                    if not (is_drop or is_raise): D_fail+=1
                    # (A)
                    if val(w,k-1)!=val(u,k-1):
                        cyc=cycle_containing(u,w,k-1)
                        if cyc is None or max(cyc)!=k-1 or astar not in cyc: A_fail+=1
                    else:
                        if astar!=k-1: A_fail+=1
                    # (B)
                    if Pset(u,vv,k-1)!=Pset(u,w,k-2): B_fail+=1
                    if (k-1-len(Pset(u,vv,k-1)))!=(k-2-len(Pset(u,w,k-2)))+1: B_fail+=1
                    # (C)
                    lhs=expand(kappa*Rkm1[vv])
                    rhs=expand(q_gs[k-1]*Rkm2[w])
                    if lhs!=rhs: C_fail+=1
    print(f"contributing instances = {tot}")
    print(f"(A) cycle top k-1 & a* in cycle : {A_fail} fail")
    print(f"(B) P & n matching              : {B_fail} fail")
    print(f"(C) weight q_{{a*,k}}*q_{{k-1}}(u,v*)==q_{{k-1}}*q_{{k-2}}(u,w): {C_fail} fail")
    print(f"(D) t_{{a*,k}} is a drop         : {D_fail} fail")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
