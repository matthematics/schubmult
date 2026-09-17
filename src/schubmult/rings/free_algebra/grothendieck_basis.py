"""`GrothendieckBasis`: the free-algebra basis dual to Grothendieck polynomials.

Keys are ``(perm, numvars)`` as in `SchubertBasis`; the key is dual to ``G_perm`` in ``numvars``
variables (``GrothendieckPolyBasis`` on the polynomial side). The change of
basis to `SchubertBasis` is the transpose of the Grothendieck-to-Schubert expansion and is
computed combinatorially: enumerate unreduced BPDs of ``perm * w0``, take the co-BPD, and keep
the reduced ones, with sign ``(-1)^(inv(perm) - inv(result))``. Products and all other
transitions route through `SchubertBasis`. ``AGx`` is the standard instance.
"""

from functools import cache

from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class GrothendieckBasis(FreeAlgebraBasis):
    """Free-algebra basis dual to Grothendieck polynomials; keys are ``(Permutation, numvars)``.

    See the module docstring. Products and transitions go through `SchubertBasis`.
    """

    zero_monom = (Permutation([]), 0)

    # @classmethod
    # def with_beta(cls, beta):
    #     """Return a configured subclass with a fixed class-level beta."""
    #     return type("GrothendieckBasis", (cls,), {"beta": beta, "zero_monom": cls.zero_monom})

    @classmethod
    def is_key(cls, x):
        """Whether ``x`` is ``(perm,)`` or ``(perm, numvars)``."""
        return (len(x) == 1 and isinstance(x[0], Permutation | list | tuple)) or (len(x) == 2 and isinstance(x[0], Permutation | list | tuple) and isinstance(x[1], int))

    @classmethod
    def as_key(cls, x):
        """Normalize to ``(Permutation, numvars)``; ``numvars`` defaults to the last descent."""
        if len(x) == 1:
            perm = Permutation(x[0])
            return (perm, 0) if len(perm.descents()) == 0 else (perm, max(perm.descents()) + 1)
        return (Permutation(x[0]), x[1])

    @classmethod
    @cache
    def transition_schubert(cls, perm, numvars):
        """Expand ``(perm, numvars)`` in `SchubertBasis`.

        For each unreduced BPD of ``perm * w0`` whose co-BPD is reduced with permutation ``u``
        fitting in ``numvars`` variables, contributes ``(-1)^(inv(perm) - inv(u))`` to ``(u, numvars)``.
        """
        from schubmult.combinatorics.bpd import BPD

        n = len(perm)
        pw0 = perm * Permutation.w0(n)

        dct = {}
        if perm.max_descent > numvars:
            return dct
        for bpd in BPD.all_unreduced_bpds(pw0, n):
            cobpd = bpd.co_bpd()
            if cobpd.is_reduced:
                the_perm = cobpd.perm
                if the_perm.max_descent <= numvars:
                    dct[(the_perm, numvars)] = dct.get((the_perm, numvars), 0) + (S.NegativeOne) ** (perm.inv - the_perm.inv)
        return dct
        # if perm.inv == 0:
        #     return {(Permutation([]), numvars): S.One}
        # from schubmult import PipeDream, RCGraph, WCGraph

        # n = len(perm)
        # # pw0 = perm * Permutation.w0(n)
        # #pw0 = perm * Permutation.w0(n)

        # dct = {}
        # if perm.max_descent > numvars:
        #     return dct
        # seen = set()
        # for pd in WCGraph.all_wc_graphs(perm, n):
        #     copd = PipeDream.from_rc_graph(pd).co_pipe_dream().to_wc_graph()
        #     if copd.is_reduced:
        #         #print(f"copd={copd} {copd.perm=} {copd.length_vector=}")
        #         the_perm = copd.perm
        #         # rc = RCGraph(copd)
        #         # if rc in seen:
        #         #     continue
        #         # print(rc)
        #         # seen.add(rc)
        #         rewc = PipeDream.from_rc_graph(copd).co_pipe_dream().to_wc_graph()
        #         if rewc != pd:
        #             continue
        #         # if rewc.resize(len(pd)) != pd:
        #         #     print("NONCOMPLIANT PD", pd, rewc.resize(len(pd)))
        #         #     continue
        #         if the_perm.max_descent <= numvars:
        #             dct[(the_perm, numvars)] = dct.get((the_perm, numvars), 0) + (S.NegativeOne) ** (the_perm.inv - perm.inv)
        # return dct

    # @classmethod
    # def transition_word(cls, perm, numvars):
    #     """Transition a Groth basis key to the word basis.

    #     Returns ``{(word, numvars): coeff}`` where ``coeff`` is the coefficient
    #     of ``G_perm`` in ``w_{word;numvars}`` on the polynomial side, interpreted
    #     as the dual-basis coefficient.
    #     """
    #     if numvars < 0:
    #         raise ValueError(f"numvars must be nonnegative, got {numvars}")
    #     if numvars == 0:
    #         return {(): S.One} if perm == Permutation([]) else {}
    #     if numvars == 1:
    #         return {(perm.inv,): S.One}
    #     if perm.max_descent > numvars:
    #         return {}
    #     cd = perm.pad_code(numvars)

    #     if cd[0] == 0:
    #         return {(0, *k): v for k, v in cls.transition_word(uncode(cd[1:]), numvars - 1).items()}
    #     wetbag = cls.product((uncode([cd[0]]),1), (uncode(cd[1:]),numvars - 1))
    #     #wetbag1 = (cd[0],)
    #     donkeydict = {(cd[0], *k): v for k, v in cls.transition_word(uncode(cd[1:]), numvars - 1).items()}
    #     #donkeydict = {cd: S.One}
    #     for (pinkbat, _), coeff in wetbag.items():
    #         if pinkbat == perm:
    #             continue
    #         wtt = cls.transition_word(uncode(pinkbat), numvars)
    #         donkeydict = {k: v for k, v in add_perm_dict_with_coeff(donkeydict, wtt, coeff=coeff).items() if v != S.Zero}
    #     return donkeydict

    # @classmethod
    # def transition_grove(cls, key):
    #     """Transition a Groth basis key to the Grove basis.

    #     Returns ``{(grove_key): coeff}`` where ``coeff`` is the coefficient
    #     of ``G_perm`` in ``g_{grove_key}`` on the polynomial side, interpreted
    #     as the dual-basis coefficient.
    #     """
    #     from schubmult.combinatorics.wc_graph import WCGraph
    #     dct = {}
    #     for wc in WCGraph.all_wc_graphs(key[0], key[1]):
    #         if wc.forest_weight == wc.length_vector:
    #             dct[wc.forest_weight] = dct.get(wc.forest_weight, S.Zero) + S.One
    #     return dct

    @classmethod
    def transition(cls, other_basis):
        """Key -> ``{key: coeff}`` function into ``other_basis``: identity on Grothendieck subclasses,
        `transition_schubert` for `SchubertBasis`, and Schubert-then-onward for everything else.
        """
        # from .elementary_basis import ElementaryBasis
        from .schubert_basis import SchubertBasis
        #from .word_basis import WordBasis

        if isinstance(other_basis, type) and issubclass(other_basis, GrothendieckBasis):
            return lambda x: {x: S.One}
        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(*x)

        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(*x))

    @classmethod
    def printing_term(cls, k):
        """Display as ``AGx(perm, numvars)``."""
        perm, numvars = cls.as_key(k)
        return GenericPrintingTerm((perm, numvars), "AGx")

    @classmethod
    @cache
    def product(cls, key1, key2, coeff=S.One):
        """Multiply by expanding both keys in `SchubertBasis`, using its separated-descents product,
        and converting the result back.
        """
        from schubmult.utils._mul_utils import add_perm_dict

        from .schubert_basis import SchubertBasis

        left = cls.transition(SchubertBasis)(key1)
        right = cls.transition(SchubertBasis)(key2)
        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, FreeAlgebraBasis.compose_transition(SchubertBasis.transition(cls), SchubertBasis.product(key_schub_left, key_schub_right, v * v2 * coeff)))
        return ret

    def __hash__(self):
        return hash((self.__class__, "Hot potato"))

    @classmethod
    def dual_basis(cls):
        """``GrothendieckPolyBasis``: Grothendieck polynomials are the dual basis."""
        from ..polynomial_algebra.grothendieck_poly_basis import GrothendieckPolyBasis
        return GrothendieckPolyBasis
