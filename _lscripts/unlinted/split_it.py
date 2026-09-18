from schubmult import *
from schubmult.utils.perm_utils import mu_A

def mu_A_perm(mu, A):
    cd = mu.trimcode
    return uncode(mu_A(cd, A))

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    w0 = Permutation.w0(n)
    grass_perms = [perm for perm in Permutation.all_permutations(n) if len(perm.descents()) == 1]
    for perm in perms:
        w = perm
        for r in range(1, len(perm.trimcode)):
            for A in itertools.combinations(range(len(perm.trimcode)), r):
                B = [i for i in range(n) if i not in A]
                mu = Permutation.w0(n)
                muA = mu_A_perm(mu, A)
                muB = mu_A_perm(mu, B)
                cprd = Sx(perm).coproduct(*[i + 1 for i in A])
                for (u, v), coeff in cprd.items():
                    prd = Sx(u*muA) * Sx(v*muB)
                    assert prd.get(w * mu, 0) == coeff, f"Failed for {perm}, {A}, {B}, {u}, {v}, {muA}, {muB}, got prd={prd}"
                    for i in range(len(A) - 1):
                        if A[i + 1] == A[i] + 1 and i in u.descents():
                            assert A[i] in w.descents()
                            u2 = u.swap(i, i + 1)
                            the_descent1 = (~(w * mu)) * (w.swap(A[i], A[i+1]) * mu) 
                            the_descent2 = (~(u * muA)) * (u2 * muA)
                            assert the_descent1 == the_descent2, f"What a pack of potatos, failed for {perm}, {A}, {B}, {u}, {v}, {muA}, {muB}, got prd={prd} and descent1={the_descent1} and descent2={the_descent2}"
                            prd2 = Sx(u2*muA) * Sx(v*muB)
                            #assert (w.swap(A[i], A[i+1]) * mu).inv == (w*mu).inv + 1
                            assert prd2.get(w.swap(A[i], A[i+1]) * mu, 0) == coeff, f"Failed for {w=}, {u=} {v=} {i=} {mu=} {muA=} {A=}, {~mu=} {v}, {muA}, {muB}, got prd={prd2}"
                            # == prd, f"Failed for {perm}, {A}, {B}, {u}, {v}, {muA}, {muB}, got prd={prd}"
        print("Love cars")


        # if perm1.is_dominant or perm2.is_dominant:
        #     continue
        # k = next(iter(perm2.descents())) + 1
        # # need k grass of n - 2
        # the_mun1 = [n - 2] * k

        # cd = [*w0.trimcode]
        # cd.insert(1, n - 2)

        # w00 = ~uncode(cd)
        # w00_cd = [*cd]
        # print(f"{cd=}")

        # cd = list(cd[:2]) + list(the_mun1) + list(cd[2:])
        # w0_2 = ~uncode(cd)
        # w0_2_cd = list(cd)
        # uw0 = perm1 * w00
        
        
        # mu = ~uncode(the_mun1)
        # vmu = perm2 * mu
        # assert vmu.inv == mu.inv - perm2.inv, f"Failed for {perm1}, {perm2}, {mu}, {vmu}"
        # v = perm2
        # prd = Sx(uw0) * Sx(vmu)
        # assert Sx(w00) * Sx(mu) == Sx(w0_2), f"Failed for {w00.trimcode}, {mu.trimcode}, {w0_2.trimcode}"
        # for ww0_2, coeff in prd.items():
        #     w = ww0_2 * (~w0_2)
        #     u = perm1
        #     if w.inv != u.inv + v.inv:
        #         continue
        #     w00 = ~uncode(w00_cd)
        #     w0_2 = ~uncode(w0_2_cd)
        #     print(f"Doing {u=} {v=} {w=} for {perm1=} {perm2=}")
        #     #new_vmu = vmu
        #     to_left = True
        #     assert u.inv + v.inv == w.inv, f"Failed for {perm1}, {perm2}, {w}, {u}, {vmu}"
        #     while True:
        #         i = max(u.descents(), default=-1)
        #         print(f"got {i=} for {u=}")
        #         if i == -1:
        #             break
        #         if i == 1:
        #             if len(u.descents()) == 1:
        #                 break
        #             else:
        #                 print("CHANGE w00 and w0_2")
        #                 w00_cd[1] = n - 1
        #                 w00 = ~uncode(w00_cd)
        #                 w0_2_cd[1] = n - 1
        #                 w0_2 = ~uncode(w0_2_cd)
        #         if i == 0:
        #             assert w[0] > w[1], f"Failed for {perm1}, {perm2}, {w}"
        #             u = u.swap(0,1)
        #             w = w.swap(0,1)
                    
        #         else:
        #             assert w[i + k + 1] > w[i + k + 2], f"Failed for {u=}, {v=}, {w=}, {i=}, {k=}"
        #             u = u.swap(i, i + 1)
        #             w = w.swap(i + k + 1, i + k + 2)
        #         assert u.inv + v.inv == w.inv, f"Failed for {perm1}, {perm2}, {w}, {u}, {vmu}"
        #     tryprod = (Sx(u * w00) * Sx(vmu))
        #     coeff2 = tryprod.get(w * w0_2, 0)
        #     assert coeff2 == coeff, f"Failed for {perm1}, {perm2}, {u*w00=}, {w* w0_2=} {vmu=} {ww0_2=} coeff={coeff}, {coeff2=} prd={prd} tryprod={tryprod}"
        #     print("Happy sptaula")
