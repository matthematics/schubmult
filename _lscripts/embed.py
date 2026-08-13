"""Two honest structural questions for  Lambda(a+b) <= Lambda(a)Lambda(b):

(A) FOREST ORDER-EMBEDDING: is T_alpha order-embedded in T_{alpha+beta}?
    T_lambda: nodes 1..N (N=len sigma_lambda), poset x<=y iff x<=y and sigma(x)<=sigma(y).

(B) RESTRICTION INJECTION: The clean route. If T_alpha and T_beta are BOTH order-embedded
    in T_gamma on COMPLEMENTARY node sets covering all of T_gamma, then
      LinExt(T_gamma) -> LinExt(T_alpha) x LinExt(T_beta),  w |-> (w|_A, w|_B)
    is injective (a linear order on gamma-nodes restricts to linear orders on each part,
    and the pair determines the merge up to the cross-comparabilities, which are FORCED
    by T_gamma) -- giving the bound. Test whether such a complementary embedding exists.
"""
from itertools import permutations, product as iproduct
from schubmult import uncode


def dom_oneline(code):
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    if not code:
        return [1]
    return list(uncode(code))


def leq(code):
    ol = dom_oneline(code)
    N = len(ol)
    M = [[ (x <= y and ol[x] <= ol[y]) for y in range(N)] for x in range(N)]
    return N, M


def embeds(codeA, codeG):
    Na, MA = leq(codeA); Ng, MG = leq(codeG)
    if Na > Ng:
        return None
    for f in permutations(range(Ng), Na):
        good = all(MA[x][y] == MG[f[x]][f[y]] for x in range(Na) for y in range(Na))
        if good:
            return f
    return None


def is_part(t):
    return all(t[i] >= t[i+1] for i in range(len(t)-1)) and all(x >= 0 for x in t)


def main():
    parts = [t for t in iproduct(range(4), repeat=3) if is_part(t) and any(t)]
    tot = embA = embBoth = 0
    failsA = []
    for a in parts:
        for b in parts:
            L = 3
            g = tuple(a[i]+b[i] for i in range(L))
            if sum(1 for _ in dom_oneline(g)) > 7:
                continue
            tot += 1
            fa = embeds(a, g)
            fb = embeds(b, g)
            if fa is not None:
                embA += 1
            if fa is not None and fb is not None:
                embBoth += 1
            else:
                failsA.append((a, b, g, fa is not None, fb is not None))
    print(f"pairs: {tot}")
    print(f"  T_alpha order-embeds in T_gamma:   {embA}/{tot}")
    print(f"  BOTH embed:                        {embBoth}/{tot}")
    print(f"  failures: {len(failsA)}")
    for a, b, g, ea, eb in failsA[:15]:
        print(f"    a={a} b={b} g={g}  embA={ea} embB={eb}")


if __name__ == "__main__":
    main()
