"""Test the inversion-set / weak-order-ideal formulation of  Lambda(a+b) <= Lambda(a)Lambda(b).

Everything for dominant sigma_lambda lives in a big S_M. Work with codes.
Key facts to test:

(1) INVERSION SPLIT: Inv(sigma_gamma) 'splits' as an alpha-block and a beta-block.
    For dominant sigma with code c, Inv has exactly c_i entries with first coord i.
    With c=gamma=alpha+beta, does Inv(sigma_gamma) decompose so that restricting a
    weak-order ideal element to the two blocks lands in [e,sigma_alpha]_L x [e,sigma_beta]_L?

(2) PAIRING MAP INJECTIVITY: define
       Phi: [e,sigma_gamma]_L -> [e,sigma_alpha]_L x [e,sigma_beta]_L
    by w |-> (w_a, w_b) where code_i(w_a)=min(code_i(w), alpha_i),
                              code_i(w_b)=code_i(w)-code_i(w_a)  (the overflow into beta).
    Is Phi well-defined (image really lands in the two ideals) and INJECTIVE?
    If injective => Lambda(gamma) <= Lambda(alpha)Lambda(beta).  Test both.

We enumerate [e,sigma_c]_L = { permutations w in S_M : Inv(w) subset Inv(sigma_c) }
via their codes: w is in the ideal iff code_i(w) <= c_i for all i AND w is a valid perm.
BUT not every code dominated by c gives Inv(w) subset Inv(sigma_c); the ideal in LEFT
weak order is {w : w <=_L sigma_c}. We enumerate honestly by inversion sets.
"""
from itertools import permutations, product as iproduct
from schubmult import uncode


def dom_oneline(code, M):
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    if not code:
        base = list(range(1, M+1))
        return base
    ol = list(uncode(code))
    ol = ol + list(range(len(ol)+1, M+1))
    return ol


def inv_set(ol):
    n = len(ol)
    return frozenset((i, j) for i in range(n) for j in range(i+1, n) if ol[i] > ol[j])


def left_ideal(code, M):
    """All w in S_M with w <=_L sigma_code, i.e. Inv(w) subset Inv(sigma_code)."""
    target = inv_set(dom_oneline(code, M))
    out = []
    for p in permutations(range(1, M+1)):
        if inv_set(p) <= target:
            out.append(p)
    return out, target


def is_part(t):
    return all(t[i] >= t[i+1] for i in range(len(t)-1)) and all(x >= 0 for x in t)


def code_of(ol):
    n = len(ol)
    return tuple(sum(1 for j in range(i+1, n) if ol[j] < ol[i]) for i in range(n))


def main():
    M = 6
    maxpart, maxlen = 3, 3
    parts = [t for t in iproduct(range(maxpart+1), repeat=maxlen) if is_part(t)]

    bad_wd = 0     # well-defined failures
    bad_inj = 0    # injectivity failures
    tested = 0
    worst = []
    for a in parts:
        for b in parts:
            L = maxlen
            g = tuple(a[i]+b[i] for i in range(L))
            if max(g) > M-1:
                continue
            Ig = inv_set(dom_oneline(g, M))
            Ia = inv_set(dom_oneline(a, M))
            Ib = inv_set(dom_oneline(b, M))
            ideal_g, _ = left_ideal(g, M)
            ideal_a_set = set(left_ideal(a, M)[0])
            ideal_b_set = set(left_ideal(b, M)[0])
            tested += 1
            images = {}
            wd_fail = inj_fail = 0
            for w in ideal_g:
                cw = code_of(w)
                cw = tuple(list(cw) + [0]*(L - len(cw)))[:L] if len(cw) < L else tuple(cw[:L])
                # split code
                ca = tuple(min(cw[i], a[i]) for i in range(L))
                cb = tuple(cw[i] - ca[i] for i in range(L))
                # reconstruct perms from codes
                wa = tuple(dom_oneline(ca, M)) if False else tuple(list(uncode(list(ca))) + list(range(len(list(uncode(list(ca))))+1, M+1))) if any(ca) else tuple(range(1, M+1))
                wb = tuple(list(uncode(list(cb))) + list(range(len(list(uncode(list(cb))))+1, M+1))) if any(cb) else tuple(range(1, M+1))
                if wa not in ideal_a_set or wb not in ideal_b_set:
                    wd_fail += 1
                key = (wa, wb)
                if key in images:
                    inj_fail += 1
                images[key] = w
            if wd_fail:
                bad_wd += 1
            if inj_fail:
                bad_inj += 1
                worst.append((a, b, g, len(ideal_g), inj_fail, wd_fail))
    print(f"tested {tested} (alpha,beta) pairs in S_{M}")
    print(f"  well-defined (image in both ideals) FAILURES: {bad_wd}")
    print(f"  injective FAILURES:                           {bad_inj}")
    for a, b, g, ng, inj, wd in worst[:12]:
        print(f"    a={a} b={b} g={g} |ideal_g|={ng} inj_collisions={inj} wd_fail={wd}")


if __name__ == "__main__":
    main()
