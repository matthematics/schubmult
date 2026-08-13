"""What is actually true for N^up(u,v) vs Lambda of various 'hull' constructions?
Compare: code-sum dominant, componentwise-max dominant, and the true right-weak-order join.
Also confirm code(join) != code(mdom u)+code(mdom v)."""
import sys
from schubmult import DSx, Permutation, uncode
from schubmult.symbolic import expand


def main(n):
    P = list(Permutation.all_permutations(n))
    inv = {w: frozenset(w.inversion_set) for w in P}
    allinv = {w: inv[w] for w in P}

    def mdom(v):
        return Permutation(list(~uncode((~v).theta())))

    # Lambda(w) = size of left-weak-order ideal [e,w]_L = # x with Inv(x) subset Inv(w)
    def Lam(w):
        iw = inv[w]
        return sum(1 for x in P if inv[x] <= iw)

    def N(u, v):
        pr = DSx(u) * DSx(v, "z")
        return len({k for k, c in pr.items() if expand(c) != 0})

    def dom_from_code(code):
        # dominant perm with given partition code, in large enough S_m
        return uncode(list(code))

    bad_sum = bad_max = bad_join = 0
    tight_sum = tight_max = tight_join = None
    code_neq = 0
    checked = 0
    for u in P:
        a = list(mdom(u).code)
        mu = mdom(u)
        for v in P:
            b = list(mdom(v).code)
            mv = mdom(v)
            Nup = N(u, mdom(v))
            # code-sum dominant
            L = max(len(a), len(b))
            aa = a + [0]*(L-len(a)); bb = b + [0]*(L-len(b))
            csum = [aa[i]+bb[i] for i in range(L)]
            cmax = [max(aa[i], bb[i]) for i in range(L)]
            Lsum = Lam(dom_from_code(csum))
            Lmax = Lam(dom_from_code(cmax))
            if Nup > Lsum:
                bad_sum += 1
            if Nup > Lmax:
                bad_max += 1
            r = Nup/Lsum
            if tight_sum is None or r > tight_sum[0]:
                tight_sum = (r, u, v, Nup, Lsum)
            r = Nup/Lmax
            if tight_max is None or r > tight_max[0]:
                tight_max = (r, u, v, Nup, Lmax)
            checked += 1
    print(f"S_{n}: checked {checked} pairs")
    print(f"  N^up <= Lambda(code-SUM dom):  {bad_sum} failures; tightest {tight_sum[0]:.3f} u={tight_sum[1]} v={tight_sum[2]} Nup={tight_sum[3]} L={tight_sum[4]}")
    print(f"  N^up <= Lambda(code-MAX dom):  {bad_max} failures; tightest {tight_max[0]:.3f} u={tight_max[1]} v={tight_max[2]} Nup={tight_max[3]} L={tight_max[4]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
