"""Type aliases for the vocabulary shared by the kernels and rings.

The aliases record the *intent* of the polymorphic entry points -- what counts as "a permutation",
"an alphabet", "a coefficient" -- since the runtime coercions (``Permutation(...)``, ``_genset``)
accept several representations.  Symbolic coefficients are SymEngine/SymPy expressions mixed
freely with Python numbers; SymEngine ships no stubs, so ``Expr`` is ``Any`` for the checker and
documentation for the reader.

The aliases are quoted and their targets imported only under ``TYPE_CHECKING``: ``from __future__
import annotations`` defers annotations, not the right-hand side of an alias assignment, and a real
import here would make this module unusable inside `schubmult.combinatorics.permutation` and
`schubmult.symbolic.poly.variables` (its own dependencies).  Modules that *use* the aliases import
them under ``TYPE_CHECKING`` too and write signatures unquoted.
"""

from __future__ import annotations

from collections.abc import Hashable, Sequence
from typing import TYPE_CHECKING, Any, TypeAlias, TypeVar

if TYPE_CHECKING:
    from schubmult.combinatorics.permutation import Permutation
    from schubmult.symbolic.poly.variables import GeneratingSet_base

__all__ = ["Alphabet", "Coeff", "CoeffDict", "Expr", "PermCoeffDict", "PermLike"]

K = TypeVar("K", bound=Hashable)

Expr: TypeAlias = Any
"""A SymEngine or SymPy expression."""

Coeff: TypeAlias = Any
"""A structure constant: an ``Expr`` or a Python ``int``/``Fraction``."""

PermLike: TypeAlias = "Permutation | Sequence[int]"
"""Anything ``Permutation(...)`` accepts: a `Permutation` or its one-line array form."""

CoeffDict: TypeAlias = dict[K, Coeff]
"""``{key: coeff}``, the expansion ``sum coeff * (basis element indexed by key)``; the key type is the
basis index of the ring in question (``Permutation``, a tuple of them, an ``RCGraph``, a word, ...)."""

PermCoeffDict: TypeAlias = "CoeffDict[Permutation]"
"""``{w: coeff_w}`` indexed by permutations: the currency of the ``schubmult.mult`` kernels."""

Alphabet: TypeAlias = "GeneratingSet_base | Sequence[Expr]"
"""A coefficient alphabet ``y_1, y_2, ...``: a generating set or a plain sequence of symbols."""
