"""AUDIT the current reflectioncancel .tex proof claims. Critical: does the proof's
assertion  v* t_{k-1,k} = w  actually hold?  (Earlier RA_via_pieriknotkn1 part (c)
suggested it FAILS in 336 cases.)  If false, the written proof is flawed and must
be redone or reverted to computational verification.

For every QUANT-nonzero Case-A instance (u(k)=w(k), u->_{k-2} w), find the unique
contributing a*<k with v*=w t_{a*,k}, u->_{k-1} v*.  Then check:
   (C1) v* t_{k-1,k} == w                          [the .tex proof's core claim]
   (C2) a* == k-1  <=>  w fixes k-1                 [when adjacency holds]
   (C3) if w moves k-1: cycle of u^{-1}w containing k-1 has TOP == k-1
   (C4) a* is a position IN that cycle (w moves) / a*==k-1 (w fixes)
   (C5) v* has: u->_{k-1} v*, NOT u->_{k-2} v*, and cycle of u^{-1}v* thru k has top k

Run: conda activate schubmult_312 && python _lscripts/audit_reflectioncancel.py 5
"""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, prod
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
    tot=0
    c1_fail=0  # v* t_{k-1,k}==w
    c3_fail=0  # cycle top (w moves)
    c4_fail=0  # a* in cycle / a*=k-1
    c5_fail=0
    wfix_count=0; wmove_count=0; astar_eq_km1_when_move=0
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
                    wfix = (val(w,k-1)==val(u,k-1))
                    # C1
                    back=vv.swap(k-2,k-1)  # v* t_{k-1,k}
                    if back!=w: c1_fail+=1
                    if wfix:
                        wfix_count+=1
                        if astar!=k-1: c4_fail+=1
                    else:
                        wmove_count+=1
                        if astar==k-1: astar_eq_km1_when_move+=1
                        cyc=cycle_containing(u,w,k-1)
                        if cyc is None or max(cyc)!=k-1: c3_fail+=1
                        if cyc is not None and astar not in cyc and astar!=k-1: c4_fail+=1
                    # C5
                    r1=(vv in Rkm1); r2=reachable(u,vv,k-2,N)
                    cyc2=cycle_containing(u,vv,k)
                    if not (r1 and not r2 and cyc2 is not None and max(cyc2)==k): c5_fail+=1
    print(f"QUANT-nonzero contributing instances = {tot}")
    print(f"(C1) v* t_{{k-1,k}} == w              : {c1_fail} FAILURES  <-- if >0, current .tex proof is WRONG")
    print(f"(C3) w-moves: cycle(k-1) top == k-1   : {c3_fail} failures")
    print(f"(C4) a* in cycle(k-1) / a*=k-1        : {c4_fail} failures")
    print(f"(C5) v* props (u->_{{k-1}},not k-2,top k): {c5_fail} failures")
    print(f"    w fixes k-1: {wfix_count} | w moves k-1: {wmove_count} (of which a*=k-1: {astar_eq_km1_when_move})")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
