"""Probabilistic zero testing for the multiplication kernels (``probabilistic=True``).

The double kernels are layered dynamic programs whose states carry unexpanded partial sums. The
signed terms cancel, and a state whose sum is zero as a polynomial -- but not structurally -- keeps
fanning out: every descendant is dead work and its leaves are output terms with coefficient zero
(6496 of them in a product of two `S_8` classes with 246 nonzero coefficients). Nothing short of
expanding can see this exactly; polynomial identity testing is only known in randomized form.

The kernels therefore carry a *shadow* of every partial sum: its value at ``points`` random points of
`\\mathbb{F}_p`, `p = 2^{61} - 1`, computed from the children's shadows as each node is built (two
modular multiplications per term, no tree walks). A state whose shadow vanishes at every point is
dropped at the end of its level. By Schwartz-Zippel a nonzero polynomial of total degree `d` vanishes
at a random point with probability at most `d / p`, so with the default two points a product of
`N` states is wrong with probability at most `N (d / p)^2` -- around `10^{-26}` for anything the
kernels can finish. The sample points are drawn afresh for every call.

:class:`ShadowEvaluator` evaluates the *inputs* of a kernel (the coefficients of ``perm_dict`` and
the factorial elementary symmetric polynomials, both small) and is what the C++ kernel calls back.
"""

import random

PRIME = (1 << 61) - 1


class ShadowEvaluator:
    """Evaluate SymEngine expressions at ``points`` random points of `\\mathbb{F}_p`, memoized by node.

    Calling it returns a tuple of ``points`` integers. Symbols get their random values on first
    sight and keep them for the lifetime of the evaluator; expressions must be built from integers,
    rationals, symbols, ``+``, ``*`` and integer powers (a negative power whose base vanishes at a
    sample point raises ``ZeroDivisionError``).
    """

    def __init__(self, points=2, seed=None):
        self.points = points
        self._rng = random.Random(seed)
        self._memo = {}
        self._symbols = {}

    def __call__(self, expr):
        return tuple(self._eval(expr, i) for i in range(self.points))

    def _eval(self, e, i):
        key = (e, i)
        v = self._memo.get(key)
        if v is not None:
            return v
        if e.is_Integer:
            v = int(e) % PRIME
        elif e.is_Symbol:
            vals = self._symbols.get(e)
            if vals is None:
                vals = self._symbols[e] = [self._rng.randrange(1, PRIME) for _ in range(self.points)]
            v = vals[i]
        elif e.is_Add:
            v = 0
            for a in e.args:
                v += self._eval(a, i)
            v %= PRIME
        elif e.is_Mul:
            v = 1
            for a in e.args:
                v = v * self._eval(a, i) % PRIME
        elif e.is_Pow:
            base, exp = e.args
            if not exp.is_Integer:
                raise TypeError(f"cannot evaluate {e}: non-integer exponent")
            b, k = self._eval(base, i), int(exp)
            if k < 0 and b == 0:
                raise ZeroDivisionError(f"{base} vanishes at the sample point")
            v = pow(b, k, PRIME)
        elif e.is_Rational:
            v = int(e.p) * pow(int(e.q), -1, PRIME) % PRIME
        else:
            raise TypeError(f"cannot evaluate {e} ({type(e).__name__}) for probabilistic zero testing")
        self._memo[key] = v
        return v
