from schubmult import Permutation, elem_sym_perms, elem_sym_perms_q, q_vector

def h_vector(q_vector):
    h = []
    val = 0
    for i in range(len(q_vector)):
        val2 = q_vector[i]
        if val2 < val:
            break
        if val2 > val:
            h += [i + 1]
        val = val2
    return tuple(h)

def l_vector(q_vector):
    """Find l_j = last position where d equals j (where d decreases from j to j-1)."""
    l = []
    val = 0
    for i in range(len(q_vector)):
        val2 = q_vector[i]
        if val2 < val:
            # Record the PREVIOUS position (i) as the last position with value val
            l += [i]
        val = val2
    return tuple(reversed(l))

def tau_d(d):
    lv = l_vector(d)
    hv = h_vector(d)

    tau = [None] * len(d)
    for i in range(len(lv)):
        if lv[i] - i >= len(d):
            tau += [None] * (lv[i] - i - len(d) + 1)
        tau[lv[i] - i] = hv[i]
    return Permutation.from_partial(tau)

def phi_d(d):
    hv = h_vector(d)
    lv = l_vector(d)

    phi = [None] * len(d)
    for i in range(len(hv)):
        if lv[i] - 1 - i >= len(d):
            phi += [None] * (lv[i] - 1 - i - len(d) + 1)
        phi[lv[i] - 1 - i] = hv[i]
    return Permutation.from_partial(phi)

# Test
if __name__ == "__main__":
    from schubmult import *
    #from schubmult.utils.schub_lib import elem_sym_perms_q
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    for k in range(1, n):
        for p in range(1, k + 1):
            for u in perms:
                print(f"{p=} {k=} {u=}")
                qlist = elem_sym_perms_q(u, p, k)
                for w, wdiff, qpow in qlist:
                    d = list(q_vector(qpow)) + [0] * (n - len(q_vector(qpow)))
                    print(f"{phi_d(d)=} {tau_d(d)=}")
                    n_i = k
                    pieri_spot = n_i - d[k - 1] if k - 1 < len(d) else n_i
                    eperms = elem_sym_perms(u * tau_d(d), pieri_spot, pieri_spot)
                    good = False
                    for up, udiff in eperms:
                        if up == w * phi_d(d):
                            good = True
                            break
                    if not good:
                        print("Test failed for n =", n, "k =", k, "p =", p)
                        print("u =", u, "w =", w, "d =", d)
                        print("w * phi_d(d) =", w * phi_d(d))
                        print("u * tau_d(d) =", u * tau_d(d))
                        print("pieri_spot =", pieri_spot)
                        print("elem_sym_perms(u * tau_d(d), pieri_spot, pieri_spot) =")
                        for up, udiff in eperms:
                            print(" ", up, "diff =", udiff)
                    else:
                        print("Test passed for n =", n, "k =", k, "p =", p, "u =", u, "w =", w, "d =", d)