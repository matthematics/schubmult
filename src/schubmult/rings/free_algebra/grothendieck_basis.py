from functools import cache

from schubmult.abc import x
from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.symbolic import S, Symbol
from schubmult.symbolic.poly.schub_poly import schub_elem_sym_to_groth_elem_sym_dict
from schubmult.symbolic.poly.variables import ZeroGeneratingSet
from schubmult.utils.perm_utils import add_perm_dict, add_perm_dict_with_coeff
from schubmult.utils.schub_lib import groth_pieri_mul

from ..printing import GrothendieckPoly
from ..schubert.grothendieck_ring import Gx
from ..schubert.separated_descents import SeparatedDescentsRing
from .free_algebra_basis import FreeAlgebraBasis

splugGx = SeparatedDescentsRing(Gx([]).ring)


class GrothendieckBasis(FreeAlgebraBasis):
    """Grothendieck basis for FreeAlgebra.

    Keys match the Schubert basis shape: ``(Permutation, numvars)``.
    The deformation parameter ``beta`` is a class variable so the class can be
    passed directly to ``FreeAlgebra`` without requiring a basis constructor.
    """

    beta = Symbol("\u03B2")
    zero_monom = (Permutation([]), 0)

    @classmethod
    def with_beta(cls, beta):
        """Return a configured subclass with a fixed class-level beta."""
        return type("GrothendieckBasis", (cls,), {"beta": beta, "zero_monom": cls.zero_monom})

    @classmethod
    def set_beta(cls, beta):
        """Set class-level beta in-place for this basis class."""
        cls.beta = beta

    @classmethod
    def is_key(cls, x):
        return (len(x) == 1 and isinstance(x[0], Permutation | list | tuple)) or (len(x) == 2 and isinstance(x[0], Permutation | list | tuple) and isinstance(x[1], int))

    @classmethod
    def as_key(cls, x):
        if len(x) == 1:
            perm = Permutation(x[0])
            return (perm, 0) if len(perm.descents()) == 0 else (perm, max(perm.descents()) + 1)
        return (Permutation(x[0]), x[1])

    @classmethod
    @cache
    def transition_from_elementary(cls, comp, numvars=None):
        """Transition an n-elementary monomial composition to Groth basis.

        If ``comp = (c1, ..., cN)`` and ``numvars=n``, interpret the factors as:
        - degree ``c_i`` in ``i`` vars for ``i <= n``
        - degree ``c_i`` in ``n`` vars for ``i > n``
        then multiply these Groth factors via ``groth_mul_full``.
        """
        comp = tuple(int(a) for a in comp)
        if any(a < 0 for a in comp):
            raise ValueError(f"Composition entries must be nonnegative, got {comp}")
        if numvars is None:
            numvars = len(comp)
        if numvars < 0:
            raise ValueError(f"numvars must be nonnegative, got {numvars}")

        if len(comp) > numvars and any(comp[i] < comp[i + 1] for i in range(numvars - 1, len(comp) - 1)):
            raise ValueError(f"Tail entries must be weakly decreasing after index {numvars}: {comp}")

        zz = ZeroGeneratingSet()
        coeff_dict = {uncode([]): S.One}
        for i, degree in enumerate(comp, start=1):
            if degree == 0:
                continue
            k = min(i, numvars)
            if degree > k:
                return {}

            elem_to_groth = schub_elem_sym_to_groth_elem_sym_dict(degree, k, x, zz, cls.beta)
            build = {}
            for (pp, kk), coeff in elem_to_groth.items():
                piece = groth_pieri_mul(coeff_dict, pp, kk, cls.beta)
                build = add_perm_dict(build, {perm: coeff * val for perm, val in piece.items()})
            coeff_dict = build

        return {(perm, numvars): coeff for perm, coeff in coeff_dict.items() if coeff != S.Zero}

    @classmethod
    def detransition(cls, degree):
        if degree < 0:
            raise ValueError(f"Degree must be nonnegative, got {degree}")
        if degree == 0:
            # G_{w_I} => beta**(|I|-1)G_w
            raise NotImplementedError("WIP")
        raise NotImplementedError("WIP")


    @classmethod
    @cache
    def transition_elementary(cls, perm, numvars):
        """Transition a Groth basis key to the dual elementary basis.

        Returns ``{(comp, numvars): coeff}`` where ``coeff`` is the coefficient
        of ``G_perm`` in ``e_{comp;numvars}`` on the polynomial side, interpreted
        as the dual-basis coefficient.

        This is computed dually by extracting the ``(perm, numvars)`` entry from
        ``transition_from_elementary(comp, numvars)`` over valid compositions
        ``comp`` with a finite search bound on ``sum(comp)`` given by the top
        homogeneous degree of ``G_perm``.
        """
        from schubmult.rings.polynomial_algebra import MonomialBasis, PolynomialAlgebra
        from schubmult.symbolic.poly.schub_poly import grothendieck_poly

        perm = Permutation(perm)
        if numvars < 0:
            raise ValueError(f"numvars must be nonnegative, got {numvars}")

        expr = grothendieck_poly(perm, x, ZeroGeneratingSet(), cls.beta)
        poly = PolynomialAlgebra(MonomialBasis(x))
        monoms = poly.from_expr(expr, length=max(numvars, len(perm) + numvars))
        top_degree = max((sum(monom) for monom, coeff in monoms.items() if coeff != S.Zero), default=0)

        ret = {}

        if numvars == 0:
            coeff0 = cls.transition_from_elementary((), 0).get((perm, 0), S.Zero)
            return ({((), 0): coeff0} if coeff0 != S.Zero else {})

        def _trim_trailing_zeros(comp):
            comp = list(comp)
            while comp and comp[-1] == 0:
                comp.pop()
            return tuple(comp)

        def _record_comp(comp):
            key_comp = _trim_trailing_zeros(comp)
            coeff = cls.transition_from_elementary(key_comp, numvars).get((perm, numvars), S.Zero)
            if coeff != S.Zero:
                ret[(key_comp, numvars)] = coeff

        prefix_len = max(numvars - 1, 0)

        def _prefixes_with_sum(total):
            # For i < numvars-1: comp[i] <= i+1, no monotonic constraint.
            if prefix_len == 0:
                if total == 0:
                    yield ()
                return

            def _rec(i, remaining, acc):
                if i == prefix_len:
                    if remaining == 0:
                        yield tuple(acc)
                    return
                bound = i + 1
                for v in range(min(bound, remaining) + 1):
                    acc.append(v)
                    yield from _rec(i + 1, remaining - v, acc)
                    acc.pop()

            yield from _rec(0, total, [])

        @cache
        def _tail_parts(total, prev):
            # For i >= numvars-1: entries are <= numvars and weakly decreasing.
            if total == 0:
                return ((),)
            out = []
            for v in range(min(prev, numvars, total), 0, -1):
                for tail in _tail_parts(total - v, v):
                    out.append((v, *tail))
            return tuple(out)

        for total in range(top_degree + 1):
            for prefix_sum in range(total + 1):
                for prefix in _prefixes_with_sum(prefix_sum):
                    rem = total - prefix_sum
                    for tail in _tail_parts(rem, numvars):
                        _record_comp((*prefix, *tail))

        return ret

    @classmethod
    def transition_word(cls, perm, numvars):
        """Transition a Groth basis key to the word basis.

        Returns ``{(word, numvars): coeff}`` where ``coeff`` is the coefficient
        of ``G_perm`` in ``w_{word;numvars}`` on the polynomial side, interpreted
        as the dual-basis coefficient.
        """
        if numvars < 0:
            raise ValueError(f"numvars must be nonnegative, got {numvars}")
        if numvars == 0:
            return {(): S.One} if perm == Permutation([]) else {}
        if numvars == 1:
            return {(perm.inv,): S.One}
        if perm.max_descent > numvars:
            return {}
        cd = perm.pad_code(numvars)

        if cd[0] == 0:
            return {(0, *k): v for k, v in cls.transition_word(uncode(cd[1:]), numvars - 1).items()}
        #wetbag = cls.product(uncode([cd[0]]), )
        wetbag = cls.product((uncode([cd[0]]), 1), (uncode(cd[1:]), numvars - 1))
        donkeydict = {cd: S.One}
        for (pinkbat, _), coeff in wetbag.items():
            if pinkbat == perm:
                continue
            wtt = cls.transition_word(uncode(pinkbat), numvars)
            donkeydict = {k: v for k, v in add_perm_dict_with_coeff(donkeydict, wtt, coeff=coeff).items() if v != S.Zero}
        return donkeydict

    @classmethod
    def transition(cls, other_basis):
        #from .elementary_basis import ElementaryBasis
        from .word_basis import WordBasis
        if other_basis == cls:
            return lambda x: {x: S.One}
        if other_basis == WordBasis:
            return lambda x: cls.transition_word(*x)
        return lambda x: FreeAlgebraBasis.compose_transition(WordBasis.transition(other_basis), cls.transition_word(*x))

    @classmethod
    def printing_term(cls, k):
        perm, numvars = cls.as_key(k)
        return GrothendieckPoly((perm, numvars), "x", prefix="A")

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        """Multiply two Schubert basis keys via the separated-descents ring."""
        return dict(coeff * splugGx(*cls.as_key(key1)) * splugGx(*cls.as_key(key2)))


    @classmethod
    def dual_basis(cls):
        from .schubert_basis import SchubertBasis

        # Placeholder until a dedicated Groth polynomial dual basis is added.
        return SchubertBasis.dual_basis()
