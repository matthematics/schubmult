from functools import cache
from typing import TYPE_CHECKING

from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict


class FreeAlgebraBasis:
    @classmethod
    def is_key(cls, x): ...

    @classmethod
    def from_rc_graph(cls, rc_graph): ...

    @classmethod
    def as_key(cls, x): ...

    @classmethod
    def transition(cls, other_basis): ...

    @classmethod
    def printing_term(cls, key): ...

    @classmethod
    def compose_transition(cls, tkeyfunc, output):
        ret = {}
        for key, v in output.items():
            ret = add_perm_dict(ret, {k: v * v0 for k, v0 in tkeyfunc(key).items()})
        return ret

    @classmethod
    def change_tensor_basis(cls, tensor_elem, basis1, basis2):
        from ..tensor_ring import TensorRing

        ring1 = tensor_elem.ring.rings[0]
        ring2 = tensor_elem.ring.rings[1]
        Tring2 = TensorRing(ring1.__class__(basis=basis1), ring2.__class__(basis=basis2))
        res = Tring2.zero
        for (key1, key2), v in tensor_elem.items():
            new_elem1 = ring1(*key1).change_basis(basis1)
            new_elem2 = ring2(*key2).change_basis(basis2)
            res += v * Tring2.ext_multiply(new_elem1, new_elem2)
        return res

    @classmethod
    @cache
    def coproduct(cls, key):
        from ...utils._mul_utils import _tensor_product_of_dicts_first
        from .word_basis import WordBasis

        return FreeAlgebraBasis.compose_transition(
            lambda x: _tensor_product_of_dicts_first(WordBasis.transition(cls)(x[0]), WordBasis.transition(cls)(x[1])),
            FreeAlgebraBasis.compose_transition(lambda y: WordBasis.coproduct(y), cls.transition(WordBasis)(key)),
        )

    @classmethod
    @cache
    def bcoproduct(cls, key):
        from ...utils._mul_utils import _tensor_product_of_dicts_first
        from .word_basis import WordBasis

        return FreeAlgebraBasis.compose_transition(
            lambda x: _tensor_product_of_dicts_first(WordBasis.transition(cls)(x[0]), WordBasis.transition(cls)(x[1])),
            FreeAlgebraBasis.compose_transition(lambda y: WordBasis.bcoproduct(y), cls.transition(WordBasis)(key)),
        )

    @classmethod
    @cache
    def product(cls, key1, key2, coeff=S.One):
        from .word_basis import WordBasis

        left = cls.transition(WordBasis)(key1)
        right = cls.transition(WordBasis)(key2)
        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, FreeAlgebraBasis.compose_transition(WordBasis.transition(cls), WordBasis.product(key_schub_left, key_schub_right, v * v2 * coeff)))
        return ret

    @classmethod
    def internal_product(cls, key1, key2, coeff=S.One):
        from .word_basis import WordBasis

        left = cls.transition(WordBasis)(key1)
        right = cls.transition(WordBasis)(key2)
        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, FreeAlgebraBasis.compose_transition(WordBasis.transition(cls), WordBasis.internal_product(key_schub_left, key_schub_right, v * v2 * coeff)))
        return ret


def __getattr__(name):
    if name == "WordBasis":
        from .word_basis import WordBasis

        return WordBasis
    if name == "JBasis":
        from .j_basis import JBasis

        return JBasis
    if name == "JTBasis":
        from .jt_basis import JTBasis

        return JTBasis
    if name == "SchubertBasis":
        from .schubert_basis import SchubertBasis

        return SchubertBasis
    if name == "SchubertSchurBasis":
        from .schubert_schur_basis import SchubertSchurBasis

        return SchubertSchurBasis
    if name == "ElementaryBasis":
        from .elementary_basis import ElementaryBasis

        return ElementaryBasis
    if name in {"SeparatedDescentsBasis", "_SeparatedDescentsBasis"}:
        from .separated_descents_basis import SeparatedDescentsBasis, _SeparatedDescentsBasis

        return {"SeparatedDescentsBasis": SeparatedDescentsBasis, "_SeparatedDescentsBasis": _SeparatedDescentsBasis}[name]
    if name == "NElementaryBasis":
        from .nelementary_basis import NElementaryBasis

        return NElementaryBasis
    if name == "ZBasis":
        from .z_basis import ZBasis

        return ZBasis
    raise AttributeError(f"module {__name__} has no attribute {name}")


if TYPE_CHECKING:
    from .elementary_basis import ElementaryBasis
    from .j_basis import JBasis
    from .jt_basis import JTBasis
    from .nelementary_basis import NElementaryBasis
    from .schubert_basis import SchubertBasis
    from .schubert_schur_basis import SchubertSchurBasis
    from .separated_descents_basis import SeparatedDescentsBasis, _SeparatedDescentsBasis
    from .word_basis import WordBasis
    from .z_basis import ZBasis

__all__ = [
    "ElementaryBasis",
    "FreeAlgebraBasis",
    "JBasis",
    "JTBasis",
    "NElementaryBasis",
    "SchubertBasis",
    "SchubertSchurBasis",
    "SeparatedDescentsBasis",
    "WordBasis",
    "ZBasis",
    "_SeparatedDescentsBasis",
]
