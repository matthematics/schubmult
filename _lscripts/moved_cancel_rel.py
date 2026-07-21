"""For each CANCEL-type q-group in the moved case, identify the structural relation between
the two contributing endpoint permutations (w1, w2) whose Pieri coefficients are the SAME
factorial-elem polynomial with opposite sign. Endpoints: for a reflection term the endpoint is
w t_{k,pp}; for QUANT the "endpoint" is w at level k-2.
Record: w1^{-1} w2 as a permutation (should be a single transposition?), and the (positions).
Goal: show CANCEL pairs are related by a single transposition t_{a,b} with matched P/n, i.e.
the reflectioncancel mechanism (adjacent-in-value swap preserving P-set & multiplicity).
Run: conda activate schubmult_312 && python _lscripts/moved_cancel_rel.py 5
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

def transp_positions(w1,w2):
    """If w1^{-1} w2 is a transposition, return the two positions swapped, else None."""
    diff=[i for i in range(1,max(len(w1),len(w2))+3) if val(w1,i)!=val(w2,i)]
    if len(diff)!=2: return None
    a,b=diff
    if val(w1,a)==val(w2,b) and val(w1,b)==val(w2,a): return (a,b)
    return None

def run(n):
    N=2*n; perms=list(Permutation.all_permutations(n))
    rel=Counter()
    for u in perms:
        u=Permutation(u)
        for k in range(2,n):
            Rkm1=enumerate_pieri(u,k-1,N)
            Rkm2=enumerate_pieri(u,k-2,N) if k-2>=0 else {u:S.One}
            Rk=enumerate_pieri(u,k,N)
            cands=set(Rk)|set(Rkm1)|set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1,N+1):
                    if pp!=k: cands.add(v0.swap(k-1,pp-1))
            for w in cands:
                w=Permutation(w)
                if val(u,k)==val(w,k): continue
                if w in Rk: continue
                for p in range(1,k+1):
                    groups=defaultdict(list)  # qmon -> (label, endpoint_perm, endpoint_level, info)
                    def add(label,qmon,ep,lev,info):
                        groups[str(expand(qmon))].append((label,ep,lev,info))
                    ci=cinfo(u,w,k-1,p-1,Rkm1)
                    if ci: add("DIAG", Rkm1[w], w, k-1, ci[1])
                    ci=cinfo(u,w,k-1,p,Rkm1)
                    if ci: add("REC", Rkm1[w], w, k-1, ci[1])
                    ci=cinfo(u,w,k-2,p-2,Rkm2)
                    if ci: add("QUANT", q_gs[k-1]*Rkm2[w], w, k-2, ci[1])
                    for pp in range(1,N+1):
                        if pp==k: continue
                        vv=w.swap(k-1,pp-1); dl=w.inv-vv.inv
                        lo,hi=min(pp,k),max(pp,k)
                        if dl==1: extra=S.One; ty="RAISE"
                        elif dl==1-2*(hi-lo): extra=q_ab(lo,hi); ty="DROP"
                        else: continue
                        ci=cinfo(u,vv,k-1,p-1,Rkm1)
                        if ci: add(ty, extra*Rkm1[vv], vv, k-1, ci[1])
                    for kk,mem in groups.items():
                        if len(mem)!=2: continue
                        # is it a CANCEL (same info)?
                        if mem[0][3]!=mem[1][3]: continue
                        (l1,e1,lv1,_),(l2,e2,lv2,_)=mem
                        tp=transp_positions(e1,e2)
                        # describe the transposition relative to k
                        if tp is None:
                            rel[("NOT-TRANSP", frozenset((l1,l2)))]+=1
                        else:
                            a,b=tp
                            # positions relative to k
                            def relpos(x): return "k" if x==k else ("<k" if x<k else ">k")
                            rel[(frozenset((l1,l2)), f"pos({relpos(a)},{relpos(b)})", f"|b-a|={b-a}")]+=1
    print("CANCEL-pair endpoint relations (w1^-1 w2 as transposition, positions rel to k):")
    for key,c in sorted(rel.items(), key=lambda x:-x[1]):
        print(f"  {key}: {c}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
