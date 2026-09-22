import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
exec(open("/home/matthematics/schubmult/_lscripts/_diamond_pipes.py").read().split("stats = Counter()")[0])
count = 0
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 6 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(1, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                LT = L(word, p)
                if LT is None:
                    continue
                m = md; wa = arr(w)
                Lpipes = {wa[m - 1], wa[m]}
                B = bumps(word, seq, p); BL = bumps(LT, seq, p)
                for (a, b), bT in B.items():
                    if bT == LT:
                        continue
                    share = len({wa[a - 1], wa[b - 1]} & Lpipes)
                    tr = transport(word, LT, a, b)
                    LbT = L(bT, p)
                    if LbT is not None and not (tr in BL and BL[tr] == LbT) and LbT != LT and LbT in BL.values() and share == 0:
                        r = [r for r, x in BL.items() if x == LbT][0]
                        wl = arr(perm_of(LT)); wb = arr(perm_of(bT)); wlb = arr(perm_of(LbT))
                        print(f"T={rows[:p]} p={p} m={m} w={''.join(map(str,wa[:8]))} beta=({a},{b}) pipes=({wa[a-1]},{wa[b-1]})  LT w={''.join(map(str,wl[:8]))}  transport={tr}  beta'={r} beta' pipes=({wl[r[0]-1]},{wl[r[1]-1]})  bT w={''.join(map(str,wb[:8]))} L(bT) w={''.join(map(str,wlb[:8]))}")
                        count += 1
                        if count >= 12:
                            raise SystemExit
