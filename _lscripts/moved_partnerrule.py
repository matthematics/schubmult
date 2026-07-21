"""Identify the DOWN partner a of the unique UP-neighbour t_{kB} by cycle surgery.

Moved case u NOT->_k w, u(k)!=w(k). pi=u^{-1}w. C=cycle thru k in NORMAL FORM (a_p,...,a_1,B),
top B=max(C), and k=a_j for some j; all other a_i<k (C has exactly two indices >=k: k,B).

up neighbour v_up = w t_{kB}; pi t_{kB} splits C into
   C'  = (a_{j-1},...,a_1,B)        [top B, the 'low' side that keeps B]
   C'' = (a_p,...,a_{j+1},k)        [top k, the 'high-index' side that keeps k]
down neighbour v_dn = w t_{a,k}, a=a_i (i!=j). Print for each pure Case(ii) bucket:
  the normal-form C, j, the up split (C',C''), the partner a=a_i and its orbit index i,
  and the two split cycles of pi t_{a_i k}. Look for the rule linking a to (C',C'').
Also print P_{k-1},n_{k-1} for up and down to confirm equality, and q-weights.
"""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet

q_gs = GeneratingSet("q")
def q_ab(a, b): return prod([q_gs[s] for s in range(a, b)])
def val(perm, pos): return perm[pos-1]

def enumerate_pieri(u, k, N):
    u = Permutation(u); results = {u: S.One}
    stack = [(u, frozenset(), N+1, S.One, u.inv)]; seen = set()
    while stack:
        perm, used_a, last_b, qw, clen = stack.pop(); key = (perm, used_a, last_b)
        if key in seen: continue
        seen.add(key)
        for a in range(1, k+1):
            if a in used_a: continue
            for b in range(k+1, N+1):
                if b > last_b: continue
                nperm = perm.swap(a-1, b-1); d = nperm.inv - clen
                if d == 1: nqw = qw
                elif d == -2*(b-a)+1: nqw = qw*q_ab(a, b)
                else: continue
                results.setdefault(nperm, nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results

def Pset(u, w, m): return frozenset(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def npart(u, w, m): return m - len(Pset(u, w, m))

def cycles_of(perm, N):
    seen = set(); cyc = []
    for x in range(1, N+1):
        if x in seen or perm[x-1] == x: continue
        c = [x]; seen.add(x); cur = perm[x-1]
        while cur != x: c.append(cur); seen.add(cur); cur = perm[cur-1]
        cyc.append(c)
    return cyc

def normal_form(cyc):
    """write cycle list as (a_p,...,a_1,b) with b=max at the END (top last)."""
    b = max(cyc); i = cyc.index(b)
    rot = cyc[i:] + cyc[:i]   # starts with b: [b, x1, x2, ...] meaning b->x1->x2..->b
    # as (a_p,...,a_1,b): the predecessors. In our convention cycle (a_p,...,a_1,b) means
    # a_p->a_{p-1}->...->a_1->b->a_p. So starting from b: b->a_p->a_{p-1}->...->a_1->b.
    # rot = [b, a_p, a_{p-1},...,a_1]; so a-sequence (a_1..a_p) = reversed(rot[1:]).
    a_seq = list(reversed(rot[1:]))  # a_1, a_2, ..., a_p
    return b, a_seq

def run(n, want=25):
    N = 2*n; perms = list(Permutation.all_permutations(n)); shown = 0
    rule_hits = Counter()
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N); Rkm1 = enumerate_pieri(u, k-1, N); Rkm2 = enumerate_pieri(u, k-2, N)
            cands = set(Rk) | set(Rkm1) | set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue
                pi = (~u)*w
                # find up neighbour B
                Cs = [c for c in cycles_of(pi, N) if k in c]
                if not Cs: continue
                C = Cs[0]; B = max(C)
                if B <= k: continue
                b_, a_seq = normal_form(C)  # b_=B, a_seq=[a_1,...,a_p]
                if k not in a_seq: continue
                j = a_seq.index(k) + 1  # k=a_j (1-based)
                for p in range(1, k+1):
                    from collections import defaultdict
                    buckets = defaultdict(list); kind = defaultdict(set)
                    if w in Rkm1:
                        qk = str(expand(Rkm1[w])); kind[qk].add('ds')
                    if w in Rkm2:
                        qk = str(expand(q_gs[k-1]*Rkm2[w])); kind[qk].add('quant')
                    for pp in range(1, N+1):
                        if pp == k: continue
                        vv = w.swap(k-1, pp-1); dl = w.inv - vv.inv
                        lo, hi = min(pp, k), max(pp, k)
                        if dl == 1: extra = S.One
                        elif dl == 1-2*(hi-lo): extra = q_ab(lo, hi)
                        else: continue
                        if vv not in Rkm1: continue
                        nm = npart(u, vv, k-1); deg = (p-1)-nm; nv = (k-1)-nm
                        if not (0 <= deg <= nv): continue
                        buckets[str(expand(extra*Rkm1[vv]))].append((pp, vv))
                    for qk, mem in buckets.items():
                        if 'ds' in kind[qk] or 'quant' in kind[qk]: continue
                        ups = [m for m in mem if m[0] > k]; downs = [m for m in mem if m[0] < k]
                        if not (len(ups) == 1 and len(downs) == 1): continue
                        (sB, v_up) = ups[0]; (a, v_dn) = downs[0]
                        i = a_seq.index(a)+1 if a in a_seq else -1
                        # hypothesized rule candidates
                        rule_hits[('i==j-1', i == j-1)] += 1
                        if shown < want:
                            shown += 1
                            print(f"u={tuple(u)} w={tuple(w)} k={k} p={p} | C=(a_p..a_1,B): "
                                  f"a_seq={a_seq} B={B} j={j} k=a_{j} | up=t_(k,{sB}) "
                                  f"down=t_({a},k) a=a_{i}")
    print(f"\nrule i==j-1 : {dict(rule_hits)}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
