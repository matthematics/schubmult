"""Classify every q-monomial group of the moved-case RHS F.
Hypothesis (clean dichotomy):
  (A) the q^0 group (M=1) is a CLASSICAL Pieri recursion RHS: it contains DIAG+REC
      and only RAISE reflections (no DROP, no QUANT); it vanishes by elemrecursion / classical Pieri.
  (B) every q^{!=1} group vanishes and contains NO DIAG/REC (only QUANT and/or reflections);
      moreover within it the multiset of (deg,nv,Pset) is balanced by sign (each poly cancels).
Report violations of:
  - M=1 group contains a DROP or QUANT  (viol_A1)
  - M!=1 group contains DIAG or REC     (viol_B1)
  - any group whose signed (deg,nv,P)-multiset does NOT cancel exactly (viol_bal)
  - any group with nonzero polynomial sum (viol_sum)  [should be 0, already known]
Also histogram group sizes for M=1 vs M!=1.
Run: conda activate schubmult_312 && python _lscripts/moved_classify.py 5
"""
import sys
from collections import defaultdict, Counter
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
def cinfo(u,w,m,a,Rm):
    if m<0: return None
    w=Permutation(w)
    if w not in Rm: return None
    Pm=Pset(u,w,m); nm=m-len(Pm); deg=a-nm; nv=m-nm
    if deg<0 or deg>nv: return None
    if nv==0:
        return (S.One,(deg,nv,tuple(Pm))) if deg==0 else None
    yv=[y[v] for v in Pm]; nc=nv+1-deg; zv=[z[i] for i in range(1,nc+1)]
    return (expand_func(FactorialElemSym(deg,nv,yv,zv)),(deg,nv,tuple(Pm)))

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    tot=0; vA1=0; vB1=0; vbal=0; vsum=0
    sizes_M1=Counter(); sizes_Mo=Counter()
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
                    groups=defaultdict(list)  # key-> list of (label, sign, poly, info)
                    def add(label,qmon,sign,poly,info):
                        groups[str(expand(qmon))].append((label,sign,poly,info))
                    ci=cinfo(u,w,k-1,p-1,Rkm1)
                    if ci: add("DIAG", q_gs[0]**0*Rkm1[w] if False else Rkm1[w], 1, (y[val(w,k)]-z[k-p+1])*ci[0], ci[1])
                    # need q-monomial of DIAG/REC = q_{k-1}(u,w) = Rkm1[w]
                    ci=cinfo(u,w,k-1,p,Rkm1)
                    if ci: add("REC", Rkm1[w], 1, ci[0], ci[1])
                    ci=cinfo(u,w,k-2,p-2,Rkm2)
                    if ci: add("QUANT", q_gs[k-1]*Rkm2[w], 1, ci[0], ci[1])
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: extra=S.One; sgn=(1 if pp>k else -1); ty="RAISE"
                        elif dl==1-2*(hi-lo): extra=q_ab(lo,hi); sgn=(1 if pp>k else -1); ty="DROP"
                        else: continue
                        ci=cinfo(u,vv,k-1,p-1,Rkm1)
                        if ci: add(ty, extra*Rkm1[vv], sgn, ci[0], ci[1])
                    # need DIAG q-monomial fix: DIAG uses Rkm1[w]
                    # (re-add DIAG properly)
                    # -- remove the bogus DIAG add above and redo
                    # rebuild groups cleanly:
                    groups=defaultdict(list)
                    ci=cinfo(u,w,k-1,p-1,Rkm1)
                    if ci: add("DIAG", Rkm1[w], 1, (y[val(w,k)]-z[k-p+1])*ci[0], ci[1])
                    ci=cinfo(u,w,k-1,p,Rkm1)
                    if ci: add("REC", Rkm1[w], 1, ci[0], ci[1])
                    ci=cinfo(u,w,k-2,p-2,Rkm2)
                    if ci: add("QUANT", q_gs[k-1]*Rkm2[w], 1, ci[0], ci[1])
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: extra=S.One; sgn=(1 if pp>k else -1); ty="RAISE"
                        elif dl==1-2*(hi-lo): extra=q_ab(lo,hi); sgn=(1 if pp>k else -1); ty="DROP"
                        else: continue
                        ci=cinfo(u,vv,k-1,p-1,Rkm1)
                        if ci: add(ty, extra*Rkm1[vv], sgn, ci[0], ci[1])
                    if not groups: continue
                    tot+=1
                    for kk,mem in groups.items():
                        labels=[m[0] for m in mem]
                        isM1 = (kk=="1")
                        if isM1: sizes_M1[len(mem)]+=1
                        else: sizes_Mo[len(mem)]+=1
                        if isM1 and any(l in ("DROP","QUANT") for l in labels): vA1+=1
                        if (not isM1) and any(l in ("DIAG","REC") for l in labels): vB1+=1
                        # signed (deg,nv,P) multiset balance
                        bal=Counter()
                        for lab,sgn,poly,info in mem: bal[info]+=sgn
                        if any(v!=0 for v in bal.values()): vbal+=1
                        if expand(sum(sgn*poly for lab,sgn,poly,info in mem))!=0: vsum+=1
    print(f"moved instances = {tot}")
    print(f"viol_A1 (M=1 group has DROP/QUANT)          : {vA1}")
    print(f"viol_B1 (M!=1 group has DIAG/REC)           : {vB1}")
    print(f"viol_bal (group signed (deg,nv,P) unbalanced): {vbal}")
    print(f"viol_sum (group poly sum != 0)              : {vsum}")
    print(f"M=1 group sizes  : {dict(sizes_M1)}")
    print(f"M!=1 group sizes : {dict(sizes_Mo)}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
