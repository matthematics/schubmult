"""MOVED case: does the full multiset of atomic terms admit a sign-reversing involution?
Atomic terms: DIAG, REC, QUANT, and each reflection term (indexed by pp).
Check: can we pair them so each pair sums to 0 (equal magnitude opposite sign)?
Necessary+sufficient for a perfect matching by value: multiset of expanded values,
grouped, each distinct nonzero value V must appear as often as -V.
Report: instances where NO such pairing exists (odd one out / value with no negative).
Run: conda activate schubmult_312 && python _lscripts/moved_involution.py 5
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

def can_pair(vals):
    # vals: list of sympy exprs (already expanded, nonzero). Return True if multiset
    # partitions into pairs summing to 0.
    key=Counter(str(expand(v)) for v in vals)
    neg=Counter(str(expand(-v)) for v in vals)
    # for each value, count(V) must equal count(-V)
    for s in list(key):
        vexpr=None
    # rebuild: map str->expr
    exprs={}
    for v in vals:
        exprs[str(expand(v))]=expand(v)
    for s,cnt in key.items():
        ns=str(expand(-exprs[s]))
        if key.get(ns,0)!=cnt:
            return False
    return True

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; nopair=0; examples=[]
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
                    terms=[]
                    d=(y[val(w,k)]-z[k-p+1])*c_term(u,w,k-1,p-1,Rkm1)
                    if expand(d)!=0: terms.append(d)
                    r=c_term(u,w,k-1,p,Rkm1)
                    if expand(r)!=0: terms.append(r)
                    qt=q_gs[k-1]*c_term(u,w,k-2,p-2,Rkm2)
                    if expand(qt)!=0: terms.append(qt)
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: tval=(1 if pp>k else -1)*c_term(u,vv,k-1,p-1,Rkm1)
                        elif dl==1-2*(hi-lo): tval=(1 if pp>k else -1)*q_ab(lo,hi)*c_term(u,vv,k-1,p-1,Rkm1)
                        else: continue
                        if expand(tval)!=0: terms.append(tval)
                    if not terms: continue
                    tot+=1
                    if not can_pair(terms):
                        nopair+=1
                        if len(examples)<5: examples.append((u,w,k,p,len(terms)))
    print(f"moved instances (nonzero terms) = {tot}")
    print(f"instances with NO exact sign-reversing pairing = {nopair}")
    for e in examples: print("   no-pair example:",e)

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
