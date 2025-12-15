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
    l = []
    val = 0
    for i in range(len(q_vector)):
        val2 = q_vector[i]
        if val2 < val:
            l += [i + 1]
        val = val2
    l += [len(q_vector) + 1]
    return tuple(l)


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
                for u, udiff, qpow in qlist:
                    d = list(q_vector(qpow))
                    print(f"  {u=} {udiff=} {d=}")
                    print(f"    h = {h_vector(d)}")
                    print(f"    l = {l_vector(d)}")