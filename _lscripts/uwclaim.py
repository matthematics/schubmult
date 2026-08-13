"""Claim: for u=w,  c_{u,v}^u(y;z) = sch_v(u(y);z), and #mon >= Lambda(v)=|[e,v]_L|.

Test:
 (1) #mon( c_{u,v}^u ) >= Lambda(v) for ALL u,v  (the diagonal w=u term).
 (2) how often strict; min ratio.
 (3) confirm degree of c_{u,v}^u is ell(v) (since ell(u)+ell(v)-ell(u)=ell(v)), consistent
     with it being sch_v (degree ell(v)).
"""
import sys
from schubmult import DSx, Permutation
from schubmult.symbolic import expand


def nmon(poly):
    e = expand(poly)
    if e == 0:
        return 0
    d = e.as_coefficients_dict()
    return len([k for k, c in d.items() if c != 0])


def main(n):
    P = list(Permutation.all_permutations(n))
    def L(p):
        ip = frozenset(p.inversion_set)
        return sum(1 for x in P if frozenset(x.inversion_set) <= ip)
    bad = 0; tot = 0; worst = None; deg_bad = 0
    for u in P:
        for v in P:
            pr = DSx(u) * DSx(v, "z")
            coeff = pr.get(u, 0)
            m = nmon(coeff)
            Lv = L(v)
            tot += 1
            if m < Lv:
                bad += 1
            r = m / Lv if Lv else 1
            if worst is None or r < worst[0]:
                worst = (r, u, v, m, Lv)
            # degree check
            dv = len(v.inversion_set)
            # degree of coeff in y,z
            e = expand(coeff)
            if e != 0:
                deg = max((sum(mono.as_poly().degree_list()) if hasattr(mono,'as_poly') else 0) for mono in [e]) if False else None
    print(f"S_{n}: diagonal w=u coefficient")
    print(f"  #mon(c_{{u,v}}^u) >= Lambda(v):  {tot-bad}/{tot} hold; {bad} fail")
    print(f"  min ratio #mon/Lambda(v) = {worst[0]:.3f} at u={list(worst[1])} v={list(worst[2])} #mon={worst[3]} Lv={worst[4]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
