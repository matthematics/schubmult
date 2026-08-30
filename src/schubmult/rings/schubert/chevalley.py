r"""Lenart--Postnikov :math:`K_T`-Chevalley formula in type :math:`A`.

Implements Theorem 6.1 / Proposition 14.5 of Lenart--Postnikov, *Affine Weyl
groups in K-theory and representation theory* (arXiv:math/0309207)::

    e^lambda . [O_u] = sum_{w, mu} c^{lambda, mu}_{u, w} x^mu [O_w]

    c^{lambda, mu}_{u, w} = sum_J (-1)^{n(J)}

summed over subsets ``J = {j_1 < ... < j_s}`` of a fixed lambda-chain such that

(a) ``u > u r_{j_1} > ... > u r_{j_1} ... r_{j_s} = w`` is a saturated
    decreasing chain in Bruhat order, and
(b) ``-mu = u r_{j_1} ... r_{j_s} (-lambda)``,

with ``n(J)`` the number of negative roots among ``beta_{j_1}, ..., beta_{j_s}``.

Weights are integer vectors in the ``epsilon`` basis; a root ``(a, b)`` denotes
``eps_a - eps_b`` and is positive iff ``a < b``.
"""

from collections import defaultdict
from fractions import Fraction
from functools import cache

from schubmult.combinatorics.permutation import Permutation

__all__ = [
    "fundamental_weight",
    "kt_chevalley_coefficients",
    "lambda_chain",
]


def fundamental_weight(k, n):
    """``omega_k = eps_1 + ... + eps_k`` as a length-``n`` vector."""
    return tuple(1 if i < k else 0 for i in range(n))


def _act(perm, vec):
    """Apply ``perm`` to a weight vector: ``perm(eps_i) = eps_{perm(i)}``."""
    out = [0] * len(vec)
    for i in range(1, len(vec) + 1):
        out[perm[i - 1] - 1] = vec[i - 1]
    return out


@cache
def lambda_chain(weight, n):
    """Reduced ``lambda``-chain for ``weight`` in ``A_{n-1}``, via Prop. 6.7.

    Returns a tuple of ``(beta, root, k)``: ``root`` is the positive root
    ``alpha`` of the affine reflection ``r_j = s_{alpha, k}``, and ``beta`` is
    the signed root ``b(r_j)`` whose sign determines ``(-1)^{n(J)}``.
    """
    entries = []
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            pairing = weight[a - 1] - weight[b - 1]
            if pairing > 0:
                ks = range(0, -pairing, -1)
            elif pairing < 0:
                ks = range(1, -pairing + 1)
            else:
                continue
            # (omega_m, alpha_ab^vee) = [a <= m] - [b <= m]
            omegas = tuple((1 if a <= m else 0) - (1 if b <= m else 0) for m in range(1, n))
            for k in ks:
                sort_key = tuple(Fraction(x, pairing) for x in (-k, *omegas))
                entries.append((sort_key, (a, b), k))
    entries.sort(key=lambda entry: entry[0])
    return tuple(((a, b) if k <= 0 else (b, a), (a, b), k) for _, (a, b), k in entries)


def kt_chevalley_coefficients(u, weight, n=None):
    """Chevalley coefficients for ``e^weight . [O_u]``.

    Returns ``{(w, mu): coefficient}`` with ``mu`` an integer weight vector.
    """
    u = Permutation(u)
    weight = tuple(weight)
    n = max(n or 0, len(weight), len(u))
    weight += (0,) * (n - len(weight))

    chain = lambda_chain(weight, n)
    length = len(chain)
    results = defaultdict(int)

    def recurse(j, w, v, translation, sign):
        if j == length:
            # mu = -u r_{j_1}...r_{j_s}(-lambda) = w(lambda) - u(translation)
            shift = _act(u, translation)
            mu = _act(w, weight)
            results[(w, tuple(x - y for x, y in zip(mu, shift)))] += sign
            return
        recurse(j + 1, w, v, translation, sign)
        beta, (a, b), k = chain[j]
        stepped = w.swap(a - 1, b - 1)
        if stepped.inv != w.inv - 1:
            return
        # appending r_j = s_{alpha, k} on the right sends (v, tau) to (v s_alpha, v(k alpha) + tau)
        if k:
            translation = list(translation)
            translation[v[a - 1] - 1] += k
            translation[v[b - 1] - 1] -= k
            translation = tuple(translation)
        recurse(j + 1, stepped, v.swap(a - 1, b - 1), translation, -sign if beta[0] > beta[1] else sign)

    recurse(0, u, Permutation([]), (0,) * n, 1)
    return {key: value for key, value in results.items() if value}
