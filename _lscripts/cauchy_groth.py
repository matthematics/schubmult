from schubmult import *
from schubmult.rings.schubert.grothendieck_ring import GrothendieckRing
from schubmult.abc import x, y
from schubmult.symbolic.common_polys import grothendieck_poly
from functools import cache

#$$F_k(x;y) = \sum_{i=0}^{k} f_{k-i}(x)(y_1^i + y_1^{i+1})$$

def _elem_perm(p, k):
    return uncode([0] * (k - p) + [1] * p)

rx = GrothendieckRing(x, beta=1)
ry = GrothendieckRing(y, beta=1)

@cache
def _elem_sym(p, k):
    return rx(_elem_perm(p, k)).expand()

@cache
def _elem_sym_s(p, k):
    return Sx(_elem_perm(p, k)).expand()


def elem_sym_by_cauchy(k):
    sm = 0
    # sm = rx(_elem_perm(k, k)).expand()
    for i in range(k + 1):
        sm += ((1 + y[1])**k) * _elem_sym_s(i, k) * ((y[1]/(1 + y[1]))**(k - i))

    return sm.expand().simplify()

def elem_sym_by_cauchy2(k):
    sm = 0
    # sm = rx(_elem_perm(k, k)).expand()
    #z = 1 / y[i]
    for i in range(k + 1):
        sm += (1 + y[1]) * _elem_sym(k - i, k) * (y[1]/(1 + y[1]))

    return sm


if __name__ == "__main__":
    from symengine import expand
    import sys

    n = int(sys.argv[1])
    for k in range(1, n):
        poly1 = elem_sym_by_cauchy(k)
        poly2 = grothendieck_poly(_elem_perm(k, k), x, y, beta=1)
        assert (poly1 - poly2).expand() == 0, f"Failed for k={k}: {poly1=}, {poly2=}"
        print(f"Success for k={k}")