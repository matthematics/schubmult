import math
from functools import cache
from schubmult import *
from schubmult.symbolic import *
from schubmult.symbolic.common_polys import grothendieck_poly

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    # for k in range(1, n + 1):
    #     schub_elem = S.Zero
    #     for i in range(k, n + 1):
    #         schub_elem += sum([Gx._beta**(i - k) * math.comb(i - 1, k - 1) * Sx(uncode([0] * (n - i) + [1]* i ))])
    #     assert expand(schub_elem.as_polynomial() - Gx(uncode([0] * (n - k) + [1] * k)).as_polynomial()) == S.Zero, f"Failed for n={n}, k={k}: {expand(schub_elem.as_polynomial() - Gx(uncode([0] * (n - k) + [1] * k)).as_polynomial())=}"
    def _groth_elem_sym(p, k, *args):
        beta = Gx._beta
        schub_elem = Sx.zero
        if p > k or k < 0:
            return schub_elem
        if p == 0:
            return Sx([])
        for i in range(p, k + 1):
            schub_elem += sum([beta**(i - p) * math.comb(i - 1, p - 1) * Sx(uncode([0] * (k - i) + [1]* i ))])
        return schub_elem

    r = NilHeckeRing(Sx.genset)

    @cache
    def strip_isobaric(index, length, genset=Sx.genset, beta=Gx._beta):
        from schubmult.abc import E
        from schubmult.rings.schubert.double_schubert_ring import DoubleSchubertRing
        # from schubmult.symoblic.poly.variables import CustomGeneratingSet
        operator = r(uncode([0] * (index - 1) + [length]))
        #ring = DoubleSchubertRing(genset, CustomGeneratingSet([0, ]))
        poly = (beta**length)*E(length, length, genset[index + 1:], [-beta**(-1)])
        schub = Sx([]).ring.from_expr(poly)#.subs({coeff_genset[i]: -coeff_genset[i] for i in range(index, index + length + 2)})
        return operator*schub
    
    def n_isobaric(perm):
        if perm.inv == 0:
            return r.one
        bacon = perm.trimcode
        result = prod([strip_isobaric(i + 1, bacon[i]) for i in range(len(bacon)) if bacon[i] != 0])
        return result
        
    coeff_genset = tuple([0 for _ in range(25)])
    w0 = Permutation.w0(n)
    for perm in perms:
        if perm.inv == 0:
            continue
        w0 = Permutation.w0(len(perm))
        schub_elem = Sx(w0)
        bacon = ((~perm)*w0).trimcode
        for i in range(len(bacon) - 1, -1, -1):
            if bacon[i] != 0:
                schub_elem = strip_isobaric(i + 1, bacon[i]).apply(schub_elem)
        assert Gx.from_expr(schub_elem.as_polynomial())==Gx(perm), f"Failed for {perm}: {Gx.from_expr(schub_elem.as_polynomial())=}, {Gx(perm)=}"
        print("Success ", perm)

## LEFT ISOBARIC DIV DIFF 1 - beta + beta y_i partial_i