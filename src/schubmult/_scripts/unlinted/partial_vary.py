from schubmult import *
from schubmult.abc import x
from schubmult.utils.perm_utils import mu_A

if __name__ == "__main__":
    n = 6
    A = [1, 3, 5]
    perms = Permutation.all_permutations(n)
    # $$\sch_w(x;y)^A = \sum_{u\leq_L w} \sch_u(x)c_{v,w_0(n)_B}^{uw_0(n)w}\sch_v(y)$$
    top_elem = 1
    for i in range(1, n):
        if i in A:
            for j in range(n - i):
                top_elem *= (Sx.from_expr(x[j + 1]) @ Sx([]) + Sx([]) @ Sx.from_expr(x[A.index(i) + 1]))
        else:
            top_elem *= Sx(uncode([1] * (n - i))) @ Sx([])
    w0 = Permutation.w0(n)
    cd = w0.trimcode
    assert len(cd) == n - 1
    w0B = mu_A(cd, [j for j in range(n - 1) if j + 1 not in A])
    w0A = mu_A(cd, [j for j in range(n - 1) if j + 1 in A])
    w0_B = uncode(w0B)
    w0_A = uncode(w0A)
    print(w0B)
    print(w0A)
    assert w0_A.inv + w0_B.inv == w0.inv
    assert Sx(w0_A) * Sx(w0_B) == Sx(w0)
    print(top_elem)
    for (u, v), coeff in top_elem.items():
        vv = v * w0_A
        prd = Sx(vv) * Sx(w0_B)
        assert (prd).get(u, 0) == coeff, f"Failed for u={u}, v={v}, vv={vv}, coeff={coeff}, prd={prd}"
    for perm in perms:
        divdiff_elem = (~perm) * w0
        perm_schub = (Sx@Sx).from_dict({(u * (~divdiff_elem), v): coeff for (u, v), coeff in top_elem.items() if (u*(~divdiff_elem)).inv == u.inv - divdiff_elem.inv})
        print(f"Perm = {perm}")
        print(perm_schub)
        makeitup_schub = 0
        
        for u in perms:
            if (u * (~divdiff_elem)).inv != u.inv - divdiff_elem.inv:
                continue
            prd = Sx(w0*u)*Sx(w0_B)
            for v, coeff in prd.items():
                if len(v) > n:
                    continue
                if ((w0 * v) * (~w0_A)).inv != w0_A.inv - (w0 * v).inv:
                    continue
                makeitup_schub += coeff * Sx(u * (~divdiff_elem)) @ Sx((w0 * v) * (~w0_A))
        assert perm_schub.almosteq(makeitup_schub), f"Failed for {perm}:\n{perm_schub}\n{makeitup_schub}"
        print("Hooray")