"""Subcase 3a: u(k)!=w(k), u->_{k-1} w, u NOT->_k w. Find the sign-reversing
pairing among ATOMIC terms {DIAG, REC, each reflection term, QUANT} that proves =0.

For each 3a instance, collect atomic terms as (label, symbolic-value). Then greedily
find pairs summing to 0, and report the residual. Also tabulate the pairing 'type'
(which label pairs with which) to reveal the rule.

Run: conda activate schubmult_312 && python _lscripts/movedvanish3a.py 5
"""
import sys
from collections import Counter
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

def run(n, show=14):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; residual_fail=0; shown=0
    pair_types=Counter()
    for u in perms:
        u=Permutation(u)
        for k in range(2,n):
            Rk=enumerate_pieri(u,k,N); Rkm1=enumerate_pieri(u,k-1,N)
            Rkm2=enumerate_pieri(u,k-2,N) if k-2>=0 else {u:S.One}
            for w in list(Rkm1):   # u->_{k-1} w
                w=Permutation(w)
                if val(u,k)==val(w,k): continue
                if w in Rk: continue  # u NOT ->_k w
                for p in range(1,k+1):
                    atoms=[]  # (label, value)
                    d=(y[val(w,k)]-z[k-p+1])*c_term(u,w,k-1,p-1,Rkm1)
                    if expand(d)!=0: atoms.append(("DIAG",expand(d)))
                    r=c_term(u,w,k-1,p,Rkm1)
                    if expand(r)!=0: atoms.append(("REC",expand(r)))
                    qy=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    if expand(qy)!=0: atoms.append(("QUANT",expand(qy)))
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        sign=1 if pp>k else -1; lo,hi=min(pp,k),max(pp,k)
                        if dl==1:
                            t=sign*c_term(u,vv,k-1,p-1,Rkm1)
                            if expand(t)!=0: atoms.append((f"R{pp}{'up' if pp>k else 'dn'}raise",expand(t)))
                        elif dl==1-2*(hi-lo):
                            t=sign*q_ab(lo,hi)*c_term(u,vv,k-1,p-1,Rkm1)
                            if expand(t)!=0: atoms.append((f"R{pp}{'up' if pp>k else 'dn'}drop",expand(t)))
                    if not atoms: continue
                    tot+=1
                    # greedy sign-reversing pairing
                    used=[False]*len(atoms)
                    pairs=[]
                    for i in range(len(atoms)):
                        if used[i]: continue
                        for j in range(i+1,len(atoms)):
                            if used[j]: continue
                            if expand(atoms[i][1]+atoms[j][1])==0:
                                used[i]=used[j]=True
                                pairs.append((atoms[i][0],atoms[j][0]))
                                break
                    residual=sum((atoms[i][1] for i in range(len(atoms)) if not used[i]), S.Zero)
                    if expand(residual)!=0:
                        residual_fail+=1
                        if shown<show:
                            print(f"RESIDUAL u={tuple(u)} p={p} k={k} w={tuple(w)}")
                            for i in range(len(atoms)):
                                mark="" if not used[i] else " [paired]"
                                print(f"    {atoms[i][0]}: {atoms[i][1]}{mark}")
                            shown+=1
                    for a,b in pairs:
                        # normalize label roles
                        def role(lbl):
                            if lbl in ("DIAG","REC","QUANT"): return lbl
                            return "Rraise" if "raise" in lbl else "Rdrop"
                        pair_types[tuple(sorted((role(a),role(b))))]+=1
    print(f"3a instances (nonzero atoms) = {tot}")
    print(f"residual (unpaired) != 0     = {residual_fail}")
    print("pairing type histogram (role--role):")
    for pt,c in sorted(pair_types.items(),key=lambda x:-x[1]): print(f"    {pt}: {c}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
