"""Pin the down-partner r* rule and the weight identity via descents.
B>k world. up v_B pieces K0=(k; c1..c_{j-1}) [top k], A0=(B; c_{j+1}..c_{L-1}) [top B].
down v_{r*} pieces K=(k; c_{i*+1}..c_{j-1}), A=(B; c1..c_{i*}, c_{j+1}..c_{L-1}).

desc_u(cycle written (a_p,...,a_1,b), b=top): positions a_l with u(a_l)>u(a_{l-1}) (a_0=b).
q-weight of a piece with top b: prod over descents a of q_{a,b}=prod_{s=a}^{b-1} q_s.

Print, for worked instances: desc & q of K0,A0 (up) and K,A (down); kappa of the down reflection;
and verify kappa_down * q_{k-1}(down) == q_{k-1}(up). Also report where r*=c_{i*} sits vs desc(K0).
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
def piece_desc_q(u,elems_top):
    # elems_top: list [top=b, then a_1, a_2, ... a_p] in the (b,a_1,...,a_p) = orbit order pi
    # written as cycle (a_p,...,a_1,b): a_l descent if u(a_l)>u(a_{l-1}), a_0=b
    b=elems_top[0]; seq=elems_top[1:]  # a_1..a_p in order
    prev=b; descs=[]
    for a in seq:
        if val(u,a)>val(u,prev): descs.append(a)
        prev=a
    qw=prod([q_ab(a,b) for a in descs]) if descs else S.One
    return descs, expand(qw)

def run(n, want=14):
    N=2*n; perms=list(Permutation.all_permutations(n)); shown=0
    allok=True; checked=0; fails=0
    for u in perms:
        u=Permutation(u)
        for k in range(2,n):
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
                # valid reflections
                ups=[]; downs=[]
                for i in range(0,j):
                    r=orb[i]; vv=w.swap(k-1,r-1); dl=w.inv-vv.inv
                    lo,hi=min(r,k),max(r,k)
                    if dl==1: kap=S.One
                    elif dl==1-2*(hi-lo): kap=q_ab(lo,hi)
                    else: continue
                    if vv not in Rkm1: continue
                    (ups if i==0 else downs).append((i,r,expand(kap),expand(Rkm1[vv])))
                if not(len(ups)==1 and len(downs)==1): continue
                iU,rU,kapU,qU = ups[0]; iD,rD,kapD,qD = downs[0]
                # up pieces: K0 = [k, c1..c_{j-1}], A0=[B, c_{j+1}..c_{L-1}]
                K0=[k]+orb[1:j]; A0=[B]+orb[j+1:L]
                # down pieces: K=[k, c_{iD+1}..c_{j-1}], A=[B, c1..c_iD, c_{j+1}..c_{L-1}]
                Kd=[k]+orb[iD+1:j]; Ad=[B]+orb[1:iD+1]+orb[j+1:L]
                dK0,qK0=piece_desc_q(u,K0); dA0,qA0=piece_desc_q(u,A0)
                dKd,qKd=piece_desc_q(u,Kd); dAd,qAd=piece_desc_q(u,Ad)
                # identity check: kappaU=1 expected (up length-incr), and kapD*qD==qU
                checked+=1
                lhs=expand(kapU*qU); rhs=expand(kapD*qD)
                if lhs!=rhs: fails+=1
                if shown<want:
                    print(f"u={tuple(u)} w={tuple(w)} k={k} B={B} j={j} orb={orb}")
                    print(f"  UP r={rU} kap={kapU} q_up={qU} | K0={K0} descK0={dK0} qK0={qK0} ; A0={A0} descA0={dA0} qA0={qA0}")
                    print(f"  DN r={rD}(c_{iD}) kap={kapD} q_dn={qD} | Kd={Kd} descKd={dKd} qKd={qKd} ; Ad={Ad} descAd={dAd} qAd={qAd}")
                    print(f"     check kapU*q_up={lhs} == kapD*q_dn={rhs} : {lhs==rhs}")
                    shown+=1
    print(f"checked pairs={checked} weight-identity fails={fails}")

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 5)
