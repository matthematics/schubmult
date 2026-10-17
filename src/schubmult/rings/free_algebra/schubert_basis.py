"""`SchubertBasis`: the free-algebra basis dual to Schubert polynomials.

A key is ``(perm, numvars)``: the element dual to ``S_perm`` viewed as a polynomial in
exactly ``numvars`` variables (so ``numvars >= max_descent(perm)``). Under the
word/monomial pairing this is the ``SchubertPolyBasis`` of the polynomial algebra.

The product is the separated-descents product (`SeparatedDescentsRing`): ``(u, p) * (v, q)``
places ``u`` in the first ``p`` variables and ``v`` in the next ``q``, giving a ``(w, p + q)``
expansion. ``transition_word`` expands a key into words via the SEM (elementary symmetric)
factorization of ``S_perm``, and the other ``transition_*`` methods reach the remaining
bases either directly or by way of the word basis. ``ASx`` is the standard instance.
"""

from functools import cache

from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.symbolic import Add, Integer, Mul, S, is_of_func_type, sympify, sympify_sympy
from schubmult.utils._lazy import LazyAttr
from schubmult.utils.perm_utils import add_perm_dict, mu_A

from ..printing import SepDescSchubPoly
from ..schubert.schubert_ring import DSx, Sx
from ..schubert.separated_descents import SeparatedDescentsRing
from .free_algebra_basis import FreeAlgebraBasis

FactorialElemSym = LazyAttr("schubmult.symbolic.symmetric_polynomials", "FactorialElemSym")

splugSx = SeparatedDescentsRing(Sx([]).ring)
ADSx = SeparatedDescentsRing(DSx([]).ring)


class SchubertBasis(FreeAlgebraBasis):
    """Free-algebra basis dual to Schubert polynomials; keys are ``(Permutation, numvars)``.

    See the module docstring. Products go through the separated-descents Schubert ring,
    and transitions to the word basis go through elementary symmetric function
    decompositions.
    """

    @classmethod
    def is_key(cls, x):
        """Whether ``x`` is ``(perm,)`` or ``(perm, numvars)`` with ``perm`` a permutation/list/tuple."""
        return (len(x) == 1 and isinstance(x[0], Permutation | list | tuple)) or (len(x) == 2 and isinstance(x[0], Permutation | list | tuple) and isinstance(x[1], int))

    @classmethod
    def from_rc_graph(cls, rc_graph):
        """The key ``(rc_graph.perm, len(rc_graph))``: an RC graph's permutation with its row count as ``numvars``."""
        return {(rc_graph.perm, len(rc_graph)): 1}

    @classmethod
    def as_key(cls, x):
        """Normalize to ``(Permutation, numvars)``; if ``numvars`` is omitted it defaults to the last descent."""
        if len(x) == 1:
            perm = Permutation(x[0])
            return (perm, 0) if len(perm.descents()) == 0 else (perm, max(perm.descents()) + 1)
        return (Permutation(x[0]), x[1])

    @classmethod
    @cache
    def product(cls, key1, key2, coeff=S.One):
        """Separated-descents product: ``(u, p) * (v, q)`` with ``u`` in the first ``p`` variables and
        ``v`` in the next ``q``, computed in `SeparatedDescentsRing`.
        """
        return dict(coeff * splugSx(*cls.as_key(key1)) * splugSx(*cls.as_key(key2)))

    zero_monom = (Permutation([]), 0)

    @classmethod
    def skew_element(cls, w, u, n):
        """The skew element ``S_w / S_u`` in ``n`` variables: the dual of multiplying by ``S_u``,
        computed with the descent-side kernel ``schubmult_py_down`` and truncated to permutations
        fitting in ``n`` variables.
        """
        from schubmult.mult.single import schubmult_py_down

        if u.inv > 0 and max(u.descents()) >= n:
            return {}
        dct = schubmult_py_down({w: S.One}, u)
        ret = {}
        for perm, v in dct.items():
            if v != S.Zero and (perm.inv == 0 or max(perm.descents()) < n):
                ret[(perm, n)] = v
        return ret

    @classmethod
    @cache
    def coproduct(cls, key):
        """Coproduct of ``(perm, numvars)`` (dual to polynomial multiplication): expand to words,
        apply the word coproduct, and convert each tensor factor back to Schubert keys.
        """
        from ...utils._mul_utils import _tensor_product_of_dicts_first
        from .word_basis import WordBasis

        dct = cls.transition_word(*key)
        res = {}
        wbasis = WordBasis
        for key_word, v in dct.items():
            dct2 = wbasis.coproduct(key_word, v)
            for (k1, k2), v2 in dct2.items():
                dct0 = wbasis.transition_schubert(k1)
                dct1 = wbasis.transition_schubert(k2)
                res = add_perm_dict(res, {k: v0 * v2 for k, v0 in _tensor_product_of_dicts_first(dct0, dct1).items()})
        return res

    @classmethod
    @cache
    def transition_grothendieck(cls, perm, numvars):
        """Expand ``(perm, numvars)`` in the `GrothendieckBasis`, by taking the co-BPD of every RC graph
        of ``perm * w0`` and collecting the resulting permutations.
        """
        from schubmult import BPD, RCGraph
        if perm.inv == 0:
            return {(Permutation([]), numvars): S.One}
        n = len(perm)
        pw0 = perm * Permutation.w0(n)
        if perm.inv == 0:
            return {(Permutation([]), numvars): S.One}
        dct = {}
        #for rc in RCGraph.all_rc_graphs(pw0, n):
        for rc in RCGraph.all_rc_graphs(pw0, n):
            bpd = BPD.from_rc_graph(rc)
            cobpd = bpd.co_bpd()
            the_perm = cobpd.perm
            if the_perm.max_descent <= numvars:
                dct[(the_perm, numvars)] = dct.get((the_perm, numvars), 0) + 1
        return dct

    @classmethod
    def transition_schubert_schur(cls, *x):
        """Expand ``(perm, numvars)`` in the `SchubertSchurBasis`: split off the variables beyond
        ``numvars`` via a Schubert coproduct against a dominant permutation, yielding
        ``(partition, perm', numvars)`` keys.
        """
        perm, numvars = x
        extra = len(perm) - numvars

        if extra <= 0:
            return {(tuple([0] * numvars), perm, numvars): 1}
        dom = uncode([numvars] * extra + list(range(numvars - 1, 0, -1)))
        tosplit = perm * dom
        dct = Sx(tosplit).coproduct(*list(range(1, extra + 1)))
        w0 = uncode(list(range(numvars - 1, 0, -1)))
        w0s = uncode([numvars] * extra)
        dct2 = {}
        for (lambd, perm0), v in dct.items():
            perm1 = perm0 * w0
            lambd2 = tuple((lambd * (~w0s)).trimcode)
            dct2[(lambd2, perm1)] = v
        return {(k[0], *([0] * (numvars - len(k[0]))), k[1]): v for k, v in dct2.items()}

    @classmethod
    def transition_schur_elementary(cls, *x):
        """Expand ``(perm, numvars)`` in the `SchurElementaryBasis` (a word-like tuple paired with a partition)."""
        from schubmult.utils.perm_utils import p_trans

        from ..polynomial_algebra import MonomialBasis, Schub
        perm, numvars = x
        # extra = len(perm) - numvars

        # if extra <= 0:
        #     return cls.transition_elementary(perm, numvars)
        # dom = uncode([numvars] * extra + list(range(numvars - 1, 0, -1)))
        # tosplit = perm * dom
        # dct = Sx(tosplit).coproduct(*list(range(1, extra + 1)))
        # w0 = uncode(list(range(numvars - 1, 0, -1)))
        # w0s = uncode([numvars] * extra)
        # dct2 = {}
        # for (lambd, perm0), v in dct.items():
        #     the_words = Schub(perm0 * w0, numvars - 1).change_basis(MonomialBasis)
        #     #elem = cls.transition_elementary(perm0 * w0, numvars - 1)
        #     lambd2 = tuple((lambd * (~w0s)).trimcode)
        #     if len(lambd2) < numvars:
        #         raise ValueError(f"Unexpected lambd2 {lambd2} from lambd {lambd} and numvars {numvars}")
        #     #     if len(lambd2)
        #     #     lambd2 = (0,) * (numvars - len(lambd2)) + lambd2
        #     if len(lambd2) > numvars:
        #         raise ValueError(f"Unexpected lambd2 {lambd2} from lambd {lambd} and numvars {numvars}")
        #     for comp, coeff in the_words.items():
        #         print(comp)
        #         key = (tuple(reversed([a for i, a in enumerate(comp)])), lambd2)
        #         dct2[key] = dct2.get(key, 0) + coeff * v
        # return dct2
        mu = p_trans(list(range(numvars - 1, 0, -1)))
        extra = len(perm) - 1 - len(mu)
        if len(mu) < len(perm) - 1:
            mu = ([numvars] * (extra)) + mu
        muw0 = uncode(mu)
        dct = Sx(perm * muw0).coproduct(*list(range(1, extra + 1)))#.change_basis(MonomialBasis)
        w0s = uncode([numvars] * extra)
        dct2 = {}
        #w0 = uncode(list(range(numvars - 1, 0, -1)))
        for (lambd, perm0), v in dct.items():
            the_words = Schub(perm0, numvars - 1).change_basis(MonomialBasis)
            lambd2 = tuple((lambd * (~w0s)).trimcode)
            if len(lambd2) == 0:
                lambd2 = (0,) * numvars
            for tup, v2 in the_words.items():
                new_tup = tuple(reversed([numvars - 1 - i - tup[i] for i in range(len(tup))]))
                dct2[(new_tup, lambd2)] = dct2.get((new_tup, lambd2), 0) + v * v2
                #ret[((tuple(reversed(new_tup[-numvars + 1 :])), *sorted(new_tup[: -numvars + 1])), numvars)] = v
        return dct2

    @classmethod
    def transition_elementary(cls, perm, numvars):
        """Expand ``(perm, numvars)`` in the `ElementaryBasis`: the coefficient of ``Elem(key)`` is the
        coefficient of ``S_perm`` in the elementary product ``E_key`` (see `ElementaryBasis.schubert_block`).
        """
        from .elementary_basis import ElementaryBasis

        return ElementaryBasis.transition_from_schubert(perm, numvars)

    @classmethod
    def transition_separated_descents(cls, k, *x):
        """Expand ``(perm, numvars)`` in the level-``k`` `SeparatedDescentsBasis` via a Schubert coproduct
        splitting the last ``k - 1`` variables, yielding ``(perm_left, perm_right, numvars)`` keys.
        """
        perm, numvars = x
        dom0_code = list(range(len(perm) - 1, 0, -1))
        dom = uncode(dom0_code)
        spot = len(dom0_code) - k + 2
        tosplit = perm * (~dom)
        dct = Sx(tosplit).coproduct(*list(range(spot, len(perm))))
        w0 = uncode(mu_A(dom.code, list(range(spot - 1, len(perm)))))
        w0s = uncode(mu_A(dom.code, list(range(spot - 1))))
        dct2 = {}
        for (perm0, perm1), v in dct.items():
            perm0_out = perm0 * (w0)
            perm1_out = perm1 * (w0s)
            dct2[(perm1_out, perm0_out, numvars)] = v
        return dct2

    @classmethod
    def transition_jbasis(cls, perm, n):
        """Expand ``(perm, n)`` in the `JBasis`: a code with no zeros is already a J key; leading zeros are
        peeled off (each contributing a factor ``t``, currently ``1``), and anything else goes via words.
        """
        from .j_basis import JBasis
        from .word_basis import WordBasis

        t = S.One
        if len(perm.trimcode) < n:
            return FreeAlgebraBasis.compose_transition(WordBasis.transition(JBasis), cls.transition_word(perm, n))
        if 0 not in perm.trimcode:
            return {tuple(perm.trimcode): S.One}
        leading_zeros = 0
        codecode = [*perm.trimcode]
        for a in perm.trimcode:
            if a == 0:
                leading_zeros += 1
                codecode = codecode[1:]
            else:
                break
        if 0 not in codecode:
            return {tuple(codecode): t ** leading_zeros}
        if leading_zeros > 0:
            return {k: v * t**leading_zeros for k, v in FreeAlgebraBasis.compose_transition(WordBasis.transition_jbasis, cls.transition_word(uncode(codecode), n - leading_zeros)).items()}
        return FreeAlgebraBasis.compose_transition(WordBasis.transition_jbasis, cls.transition_word(perm, n))

    @classmethod
    def dual_basis(cls):
        """``SchubertPolyBasis``: Schubert polynomials are the dual basis under the word/monomial pairing."""
        from ..polynomial_algebra.schubert_poly_basis import SchubertPolyBasis
        return SchubertPolyBasis

    @classmethod
    @cache
    def transition(cls, other_basis):
        """Return the key -> ``{key: coeff}`` function into ``other_basis``.

        Direct routes exist for the word, elementary, Schubert-Schur, Schur-elementary,
        composition-Schubert, separated-descents, and Grothendieck bases; everything else
        is reached by going through the word basis first.
        """
        from .composition_schubert_basis import CompositionSchubertBasis
        from .elementary_basis import ElementaryBasis
        from .forest_basis import ForestBasis
        from .fundamental_slide_basis import FundamentalSlideBasis
        from .glide_basis import GlideBasis
        from .grothendieck_basis import GrothendieckBasis
        from .grove_basis import GroveBasis
        from .j_basis import JBasis
        from .jt_basis import JTBasis
        from .key_basis import KeyBasis
        from .lascoux_basis import LascouxBasis
        from .monomial_slide_basis import MonomialSlideBasis
        from .schubert_schur_basis import SchubertSchurBasis
        from .schur_elementary_basis import SchurElementaryBasis
        from .word_basis import WordBasis
        from .z_basis import ZBasis

        if other_basis == SchubertBasis:
            return lambda x: {x: S.One}
        if other_basis == SchurElementaryBasis:
            return lambda x: cls.transition_schur_elementary(*x)
        if other_basis == CompositionSchubertBasis:
            return lambda x: {CompositionSchubertBasis.as_key(x): S.One}
        if other_basis == ElementaryBasis:
            return lambda x: cls.transition_elementary(*x)
        if other_basis == SchubertSchurBasis:
            return lambda x: cls.transition_schubert_schur(*x)
        if other_basis == WordBasis:
            return lambda x: cls.transition_word(*x)
        if other_basis.__name__ == "_SeparatedDescentsBasis":
            return lambda x: cls.transition_separated_descents(other_basis.k, *x)
        if other_basis == ZBasis or other_basis == JTBasis or other_basis == JBasis or other_basis == MonomialSlideBasis or other_basis == ForestBasis or other_basis == KeyBasis or other_basis == FundamentalSlideBasis:
            return lambda x: FreeAlgebraBasis.compose_transition(WordBasis.transition(other_basis), cls.transition_word(*x))
        if isinstance(other_basis, type) and issubclass(other_basis, GrothendieckBasis):
            return lambda x: cls.transition_grothendieck(*x)
        if other_basis == GroveBasis or other_basis == GlideBasis or other_basis == LascouxBasis:
            return lambda x: FreeAlgebraBasis.compose_transition(WordBasis.transition(other_basis), cls.transition_word(*x))
        raise NotImplementedError(f"Transition from SchubertBasis to {other_basis} is not implemented.")

    @classmethod
    @cache
    def old_transition_word(cls, perm, numvars):
        """Transition to WordBasis via SEM basis (legacy implementation)."""
        res = {}
        expr = Sx(perm * ~uncode(list(range(perm.inv + numvars, perm.inv, -1)))).in_SEM_basis().expand()
        args = expr.args
        if not isinstance(expr, Add):
            args = [expr]
        for arg in args:
            tup = list(range(perm.inv + numvars, perm.inv, -1))
            coeff = S.One
            if is_of_func_type(sympify(arg), FactorialElemSym):
                arg = sympify_sympy(arg)
                tup[perm.inv + numvars - arg.numvars] = arg.numvars - arg.degree
            elif isinstance(arg, Mul):
                for arg0 in arg.args:
                    if is_of_func_type(sympify(arg0), FactorialElemSym):
                        arg0 = sympify_sympy(arg0)
                        tup[perm.inv + numvars - arg0.numvars] = arg0.numvars - arg0.degree
                    else:
                        coeff = Integer(arg0)
            else:
                coeff = Integer(arg)
            tup = tuple(tup)
            res[tup] = res.get(tup, S.Zero) + coeff
        return res

    @classmethod
    @cache
    def transition_word(cls, perm, numvars):
        """Expand ``(perm, numvars)`` in the word basis.

        Multiplies ``perm`` by the inverse of the dominant permutation with code
        ``(inv + numvars, ..., inv + 1)`` (``inv = perm.inv``) and writes the result in the SEM
        basis with a custom ``elem_func`` that records each factor ``e_p(x_1..x_k)`` as the word
        with ``k - p`` in position ``numvars - k + inv``. Products of such words in the monomial
        polynomial algebra add letterwise, so the resulting polynomial *is* the word expansion.
        """
        from ..polynomial_algebra._core import PA

        def word_elem(p, k, *args):  # noqa: ARG001
            """Encode ``e_p(x_1..x_k)`` as a single word with entry ``k - p``; 0/1 for the degenerate cases."""
            import numpy as np
            if p > k or p < 0:
                return 0
            if p == 0:
                return 1
            vec = np.zeros(numvars, dtype=int)
            vec[numvars - k + perm.inv] = k - p
            return PA.from_dict({tuple(vec.tolist()): 1})
        #ret = {k: v for k, v in Sx(perm * ~uncode(list(range(perm.inv + numvars, perm.inv, -1)))).in_SEM_basis(elem_func=word_elem)#.items()}
        result = Sx(perm * ~uncode(list(range(perm.inv + numvars, perm.inv, -1)))).in_SEM_basis(elem_func=word_elem)
        if numvars == 0:
            if int(result) == 0:
                return {}
            return {(): result}
        return dict(result)

    @classmethod
    def printing_term(cls, k):
        """Display symbol for ``(perm, numvars)`` (the separated-descents ``SepDescSchubPoly`` form)."""
        return SepDescSchubPoly(cls.as_key(k), None, None)
