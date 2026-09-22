"""Print failures where bumping isomorphic hw pairs at the same position gives non-isomorphic results."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
exec(open("/home/matthematics/schubmult/_lscripts/_cell_iso.py").read().split("stats = Counter()")[0])
shown = 0
for key, lst in groups.items():
    p = key[0]
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            (wa, sa, u), (wb, sb, v) = lst[i], lst[j]
            RA, RB = root_tab(wa), root_tab(wb)
            if set(RA) != set(RB): continue
            cA, cB = cells(wa), cells(wb)
            for c in RA:
                ra, rb = RA[c], RB[c]
                if not (ra[0] > p and is_cover(u, *ra) and rb[0] > p and is_cover(v, *rb)): continue
                na, SA = bump_cells(wa, cA[c]); nb, SB = bump_cells(wb, cB[c])
                if na is None or nb is None or not valid_rc(na, sa) or not valid_rc(nb, sb): continue
                TA = RCGraph.from_reduced_compatible(list(na), list(sa)); TB = RCGraph.from_reduced_compatible(list(nb), list(sb))
                if not (is_hw(TA, p) and is_hw(TB, p)): continue
                if extwt(TA, p) != extwt(TB, p):
                    print(f"p={p} key={key}  T={wa} rows={sa} u={arr(u,7)}  T'={wb} rows={sb} v={arr(v,7)}  cell={c} roots={ra},{rb}  -> {na} extwt={extwt(TA,p)} ; {nb} extwt={extwt(TB,p)}")
                    shown += 1
                    if shown >= 6: raise SystemExit
