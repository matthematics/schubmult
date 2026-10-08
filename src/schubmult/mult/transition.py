"""Ordinary Schubert products by transition + Monk, and the cost-model hybrid with `schubmult_py`.

This is the algorithm of Buch's ``lrcalc``: expand one factor into monomials with the
Lascoux--Schützenberger transition recursion ``S_w = x_r S_v + sum_i S_{v t_{ir}}`` (``r`` the last
descent of ``w``), then multiply the monomials into the other factor one variable at a time by Monk's
rule, as a Horner scheme over the Schubert basis.  Its cost is governed by the number of pipe dreams
of the expanded factor, where the v-path kernel's is governed by the size of its v-path structure;
the two are complementary, and `schubmult_py_hybrid` picks between them with a cost model fitted on
exhaustive timings over ``S_7`` (see ``cpp/kernel_transition.h``).

The compiled kernels in ``schubmult_cpp`` do the work; the Python implementations here are the
reference and the fallback for permutations beyond the extension's ``MAXN``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.mult import _accel
from schubmult.mult.single import schubmult_py, single_variable

if TYPE_CHECKING:
    from schubmult._typing import PermCoeffDict, PermLike

__all__ = ["mult_monomials_py", "pipe_dream_count", "schubmult_py_hybrid", "schubmult_py_transition", "transition_monomials"]

Monomial = tuple[int, ...]
"""Exponent vector ``(a_1, a_2, ...)`` of ``x_1^{a_1} x_2^{a_2} ...``, trailing zeros trimmed."""


def _trim(m: list[int]) -> Monomial:
    while m and m[-1] == 0:
        m.pop()
    return tuple(m)


def transition_monomials(w: PermLike) -> dict[Monomial, int]:
    """The monomial expansion of ``S_w`` as ``{exponent vector: coefficient}`` by the transition recursion."""
    arr = list(Permutation(w))
    out: dict[Monomial, int] = {}
    _trans_rec(arr, out)
    return out


def _trans_rec(w: list[int], out: dict[Monomial, int]) -> None:
    n = len(w)
    r = next((i for i in range(n - 1, 0, -1) if w[i - 1] > w[i]), 0)  # last descent, 1-indexed
    if r == 0:
        out[()] = out.get((), 0) + 1
        return
    s = r + 1
    while s < n and w[r - 1] > w[s]:
        s += 1
    v = list(w)
    v[r - 1], v[s - 1] = v[s - 1], v[r - 1]
    sub: dict[Monomial, int] = {}
    _trans_rec(v, sub)
    for m, c in sub.items():
        mm = list(m) + [0] * max(0, r - len(m))
        mm[r - 1] += 1
        key = tuple(mm)
        out[key] = out.get(key, 0) + c
    vr, last = v[r - 1], 0
    for i in range(r - 1, 0, -1):
        vi = v[i - 1]
        if last < vi < vr:
            last = vi
            nxt = list(v)
            nxt[i - 1], nxt[r - 1] = nxt[r - 1], nxt[i - 1]
            _trans_rec(nxt, out)


def mult_monomials_py(perm_dict: PermCoeffDict, monomials: dict[Monomial, int]) -> PermCoeffDict:
    """``(sum_u coeff_u S_u) * (sum of monomials)`` in the Schubert basis, Horner over the variables with Monk's rule."""
    ret: PermCoeffDict = {}
    _horner([(m, c) for m, c in monomials.items() if c], max((len(m) for m in monomials), default=0), perm_dict, ret)
    return {k: c for k, c in ret.items() if c != 0}


def _horner(terms: list[tuple[Monomial, int]], maxvar: int, perm_dict: PermCoeffDict, out: PermCoeffDict) -> None:
    if not terms:
        return
    if maxvar == 0:
        c = sum(t[1] for t in terms)
        for u, cu in perm_dict.items():
            out[u] = out.get(u, 0) + c * cu
        return
    lower = [t for t in terms if len(t[0]) < maxvar]
    upper = [(_trim([*m[: maxvar - 1], m[maxvar - 1] - 1]), c) for m, c in terms if len(m) == maxvar]
    res: PermCoeffDict = {}
    _horner(upper, max((len(m) for m, _ in upper), default=0), perm_dict, res)
    for u, c in single_variable(res, maxvar).items():
        out[u] = out.get(u, 0) + c
    _horner(lower, max((len(m) for m, _ in lower), default=0), perm_dict, out)


def pipe_dream_count(w: PermLike) -> int:
    """``S_w(1, ..., 1)``, the number of pipe dreams (RC graphs) of ``w``."""
    if _accel.available:
        ret = _accel.pipe_dream_count(w)
        if ret is not None:
            return ret
    from schubmult.combinatorics.rc_graph import RCGraph

    return RCGraph.count_rc_graphs(Permutation(w))


def schubmult_py_transition(perm_dict: PermCoeffDict, v: PermLike) -> PermCoeffDict:
    """``(sum_u coeff_u S_u) * S_v`` by expanding ``S_v`` into monomials and multiplying them in by Monk's rule.

    Same contract as `schubmult_py`; faster when ``v`` has few pipe dreams relative to the size of
    its v-path structure (many descents, large theta entries), slower otherwise.
    """
    if _accel.available:
        ret = _accel.schubmult_py_transition(perm_dict, v)
        if ret is not None:
            return ret
    return mult_monomials_py(perm_dict, transition_monomials(v))


def _vpath_count(v: Permutation) -> int:
    from schubmult.utils.schub_lib import compute_vpathdicts

    th = list((~v).theta())
    while th and th[-1] == 0:
        th.pop()
    if not th:
        return 0
    vpd = compute_vpathdicts(th, v * uncode(th))
    return sum(len(steps) for layer in vpd for steps in layer.values())


def _prefers_transition(pd: int, nvp: int) -> bool:
    # fitted on exhaustive S_7 and S_8..S_12 family timings; see hybrid_choose in cpp/kernel_transition.h
    return (1 + pd) ** 0.6 < 0.247 * (1 + nvp) ** 1.1


def schubmult_py_hybrid(perm_dict: PermCoeffDict, v: PermLike) -> PermCoeffDict:
    """``(sum_u coeff_u S_u) * S_v`` by whichever of `schubmult_py` and `schubmult_py_transition` a cost model predicts to be faster.

    For a single ``u`` either factor may be expanded or recursed on; for several only ``v`` is.
    """
    if _accel.available:
        ret = _accel.schubmult_py_hybrid(perm_dict, v)
        if ret is not None:
            return ret
    v = Permutation(v)
    terms = [(u, c) for u, c in perm_dict.items() if c != 0]
    if not terms:
        return {}
    if len(terms) == 1:
        u, c = terms[0]
        candidates = [(Permutation(u), v), (v, Permutation(u))]  # (expanded or recursed on, other)
        a, b = min(candidates, key=lambda p: _vpath_count(p[0]))
        e = min(candidates, key=lambda p: pipe_dream_count(p[0]))
        if _prefers_transition(pipe_dream_count(e[0]), _vpath_count(a)):
            ret = mult_monomials_py({e[1]: 1}, transition_monomials(e[0]))
        else:
            ret = schubmult_py({b: 1}, a)
        return {w: c * val for w, val in ret.items()} if c != 1 else ret
    if _prefers_transition(pipe_dream_count(v), _vpath_count(v)):
        return schubmult_py_transition(perm_dict, v)
    return schubmult_py(perm_dict, v)
