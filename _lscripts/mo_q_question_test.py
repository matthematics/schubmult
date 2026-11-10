from schubmult import *
from schubmult.rings.schubert import *
from schubmult.abc import *
from schubmult.mult.quantum_double import *
from schubmult.symbolic import *
import numpy as np

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    w0 = Permutation.w0(n)
    dring = DSx([]).ring
    qn1 = ChangeOfVariables(z, q, lambda index, gset: gset[index]**(-1))
    QSy = QuantumSingleSchubertRing(y, q_var = qn1)
    sm = sum([QSx(perm).as_polynomial() * QSy(perm*w0).as_polynomial() for perm in perms])
    Sy = SingleSchubertRing(y)
    #da_full_elem = dring.from_expr(sm)
    da_elem = Sx.from_expr(sm).simplify()
    #(Sx([])*da_full_elem).simplify()
    #print(da_elem.simplify())
    w00 = Permutation.w0(n - 1)
    # #3w00 = w0
    da_elem = Sx.from_dict({k * w00: v for k, v in da_elem.items() if (k * w00).inv == k.inv - w00.inv})
    #print(da_elem.simplify())
    as_tensor = (Sx @ Sx @ Sy).zero
    for perm, coeff in da_elem.items():
        q_dict = factor_out_q(coeff)
        for pow_of_q, val in q_dict.items():
            if expand(pow_of_q) == 0:
                continue
            yploy = Sy.from_expr(val)
            # for perm2, coeff2 in yploy.items():
            #     perm22 = w00 * perm2
            #     if perm22.inv != perm2.inv - w00.inv:
            #         continue
            as_tensor += ((pow_of_q) * Sx([]))@ Sx(perm) @ yploy

    # #print(da_elem.expand(deep=False).simplify())
    print(as_tensor)