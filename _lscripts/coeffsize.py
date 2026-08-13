"""NEW TOPIC: representation size of the double coefficient c_{u,v}^w(y;z).
Count monomials of each coefficient and compare to candidate weak-order interval sizes.

DSx(u)*DSx(v,'z') -> dict {w: c_{u,v}^w(y;z)}.  #monomials = terms of expand(coeff).
degree d(w) = ell(u)+ell(v)-ell(w).

Candidate interval sizes to compare against #mon(c_{u,v}^w):
  Lambda(v), Lambda(u), Lambda(w),
  |[e,w]_R|-ish via Lambda,
  binom(d + (#vars-1), ...) trivial,
  |[v,w]_Bruhat|, |[u,w]_Bruhat|.
Just print a table first to see the shape of the data.
"""
import sys
from itertools import permutations
from schubmult import DSx, Permutation, uncode
from schubmult.symbolic import expand


def nmon(poly):
    e = expand(poly)
    if e == 0:
        return 0
    d = e.as_coefficients_dict()
    return len([k for k,c in d.items() if c != 0])


def ell(p):
    return p.number_of_inversions() if hasattr(p,'number_of_inversions') else len(p.inversion_set)


def main(n):
    P = list(Permutation.all_permutations(n))
    def L(p):
        ip = frozenset(p.inversion_set)
        return sum(1 for x in P if frozenset(x.inversion_set) <= ip)
    rows = []
    for u in P:
        for v in P:
            pr = DSx(u)*DSx(v,"z")
            for w,coeff in pr.items():
                m = nmon(coeff)
                if m == 0:
                    continue
                du = len(u.inversion_set); dv=len(v.inversion_set); dw=len(w.inversion_set)
                deg = du+dv-dw
                rows.append((u,v,w,m,deg,L(u),L(v),L(w)))
    # print a sample where deg>0 and m>1
    rows.sort(key=lambda r:-r[3])
    print(f"S_{n}: {len(rows)} nonzero coeffs; showing largest #mon")
    print(f"{'u':<12}{'v':<12}{'w':<12}{'#mon':>5}{'deg':>4}{'L(u)':>5}{'L(v)':>5}{'L(w)':>5}")
    for u,v,w,m,deg,Lu,Lv,Lw in rows[:25]:
        print(f"{str(list(u)):<12}{str(list(v)):<12}{str(list(w)):<12}{m:>5}{deg:>4}{Lu:>5}{Lv:>5}{Lw:>5}")
    # check simple relations
    print("\nchecks over all coeffs:")
    ge_Lv = sum(1 for r in rows if r[3] >= r[6])
    print(f"  #mon >= L(v):  {ge_Lv}/{len(rows)}")
    ge_Lu = sum(1 for r in rows if r[3] >= r[5])
    print(f"  #mon >= L(u):  {ge_Lu}/{len(rows)}")


if __name__=="__main__":
    main(int(sys.argv[1]) if len(sys.argv)>1 else 3)
