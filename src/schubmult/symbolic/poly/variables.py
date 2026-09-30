"""Generating sets: indexed families of variables ``x_0, x_1, x_2, ...`` used as ring generators.

`GeneratingSet("x")` interns ``DEF_GENSET_SIZE`` symbols ``x_0..x_99``; ``gs[i]`` is the symbol
``x_i``, and the polynomial variables are ``x_1, x_2, ...`` (``x_0`` is unused), so an exponent
tuple ``(a_1, ..., a_n)`` means ``x_1^{a_1} ... x_n^{a_n}``. `MaskedGeneratingSet` hides a set of indices of a base set,
`CustomGeneratingSet` wraps an arbitrary sequence of expressions, and `ZeroGeneratingSet`
returns 0 for every index (used for single Schubert polynomials as a degenerate coefficient
set). `genset_dict_from_expr` converts a polynomial expression into ``{exponent_tuple: coeff}``.
"""

# class generators with base
# symbols cls argument!

from __future__ import annotations

import re
from bisect import bisect_left
from collections.abc import Iterable, Iterator, Sequence
from functools import cache
from typing import TYPE_CHECKING, Any, ClassVar, overload

from schubmult.symbolic import Add, Mul, Pow, S, SympifyError, expand, symbols, sympify
from schubmult.utils.logging import get_logger

if TYPE_CHECKING:
    from schubmult._typing import Expr

logger = get_logger(__name__)

DEF_GENSET_SIZE = 100


class GeneratingSet_base:
    """Interface for generating sets: indexing, length, ``index(symbol)`` (``-1`` if absent), and ``label``."""

    _args: tuple[Any, ...]

    def __new__(cls, *args: Any) -> GeneratingSet_base:
        obj = object.__new__(cls)
        obj._args = args
        return obj

    @property
    def args(self) -> tuple[Any, ...]:
        return self._args

    @overload
    def __getitem__(self, i: int) -> Expr: ...

    @overload
    def __getitem__(self, i: slice) -> Sequence[Expr]: ...

    def __getitem__(self, i):
        return NotImplemented

    def __len__(self) -> int:
        return NotImplemented

    def __iter__(self) -> Iterator[Expr]:
        yield from [self[i] for i in range(len(self))]

    def index(self, other: object) -> int:
        """Position of the symbol ``other`` in this set, or ``-1``."""
        raise NotImplementedError

    def __contains__(self, other: object) -> bool:
        return self.index(other) != -1

    @property
    def label(self) -> str | None:
        return None


class ZeroGeneratingSet(GeneratingSet_base):
    """A generating set every entry of which is ``0``; contains no symbols."""

    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.stop is None:
                return self
            start = index.start if index.start is not None else 0
            stop = index.stop
            return [S.Zero for i in range(start, stop)]
        return S.Zero

    def __contains__(self, item: object) -> bool:
        return False

    def index(self, _: object) -> int:
        return -1

    def __iter__(self) -> Iterator[Expr]:
        if False:
            yield


# variable registry
# TODO: ensure sympifies
# TODO: masked generating set
class GeneratingSet(GeneratingSet_base):
    """The interned family ``name_0, name_1, ...``; ``gs[i]`` is the symbol ``name_i`` and ``gs(i)`` is ``gs[i - 1]``."""

    _symbols_arr: tuple[Expr, ...]
    _index_lookup: dict[Any, int]
    _hash: int

    def __new__(cls, name: str) -> GeneratingSet:
        return GeneratingSet.__xnew_cached__(cls, name)  # type: ignore[arg-type]  # mypy: type[...] vs Hashable (python/mypy#11470)

    _registry: ClassVar[dict[str, GeneratingSet]] = {}

    _index_pattern = re.compile("^([^_]+)_([0-9]+)$")
    _sage_index_pattern = re.compile("^([^0-9]+)([0-9]+)$")

    # is_Atom = True
    # TODO: masked generating set
    @staticmethod
    @cache
    def __xnew_cached__(_class, name):
        return GeneratingSet.__xnew__(_class, str(name))

    @staticmethod
    def __xnew__(_class, name):
        obj = GeneratingSet_base.__new__(_class, name)
        obj._symbols_arr = tuple([symbols(f"{name}_{i}") for i in range(DEF_GENSET_SIZE)])
        obj._index_lookup = {obj._symbols_arr[i]: i for i in range(len(obj._symbols_arr))}
        obj._hash = hash(name)
        return obj

    # def shift(self, index):
    #     return CustomGeneratingSet(self._symbols_arr[])

    def __call__(self, index: int) -> Expr:
        """1-indexed"""
        return self[index - 1]

    @property
    def label(self) -> str:
        """The variable name, e.g. ``"x"``."""
        return str(self.args[0])

    # index of v in the genset
    def index(self, v: object) -> int:
        """Position of the symbol ``v`` in this set, or ``-1``."""
        try:
            return self._index_lookup.get(v, self._index_lookup.get(sympify(v), -1))
        except SympifyError:
            return -1
        except TypeError:
            return -1

    def __repr__(self):
        return f"GeneratingSet('{self.label}')"

    def __str__(self):
        return self.label

    def _latex(self, printer):
        return printer.doprint(self.label)

    def _sympystr(self, printer):
        return printer.doprint(self.label)

    @overload
    def __getitem__(self, i: int) -> Expr: ...

    @overload
    def __getitem__(self, i: slice) -> tuple[Expr, ...]: ...

    def __getitem__(self, i):
        return self._symbols_arr[i]

    def __len__(self) -> int:
        return len(self._symbols_arr)

    def __hash__(self) -> int:
        return self._hash

    def __iter__(self) -> Iterator[Expr]:
        yield from self._symbols_arr

    def __eq__(self, other: object) -> bool:
        return self is other or (isinstance(other, GeneratingSet) and self.label == other.label)


class MaskedGeneratingSet(GeneratingSet_base):
    """A base generating set with the (1-indexed) positions in ``index_mask`` removed and the rest
    renumbered consecutively; ``complement()`` gives the set of the masked variables instead.
    """

    _mask: dict[int, int]
    _index_lookup: dict[Any, int]
    _label: str
    _symbols_arr: tuple[Expr, ...]

    def __new__(cls, gset: GeneratingSet_base, index_mask: Iterable[int]) -> MaskedGeneratingSet:
        return MaskedGeneratingSet.__xnew_cached__(cls, gset, tuple(sorted(index_mask)))  # type: ignore[arg-type]  # mypy: type[...] vs Hashable (python/mypy#11470)

    @staticmethod
    @cache
    def __xnew_cached__(_class, gset, index_mask):
        return MaskedGeneratingSet.__xnew__(_class, gset, index_mask)

    @staticmethod
    def __xnew__(_class, gset, index_mask):
        obj = GeneratingSet_base.__new__(_class, gset, index_mask)
        # obj._symbols_arr = tuple([symbols(f"{name}_{i}") for i in range(100)])
        # obj._index_lookup = {obj._symbols_arr[i]: i for i in range(len(obj._symbols_arr))}
        mask_dict = {}
        mask_dict[0] = 0
        for i in range(1, len(gset._symbols_arr)):
            index = bisect_left(index_mask, i)
            # logger.debug(f"{index=}")
            if index >= len(index_mask) or index_mask[index] != i:
                # logger.debug(f"{i - index} mapsto {i} and {index_mask=}")
                mask_dict[i - index] = i
            # if index>=len(index_mask) or index_mask[index] != i:
            #     mask_dict[cur_index] = i
            #     cur_index += 1
        # print(f"{index_mask=} {mask_dict=}")
        obj._mask = mask_dict
        obj._index_lookup = {gset[mask_dict[i]]: i for i in range(len(gset) - len(index_mask))}
        obj._label = gset.label
        obj._symbols_arr = tuple(iter(obj))
        return obj

    @property
    def base_genset(self) -> GeneratingSet_base:
        """The underlying unmasked generating set."""
        return self.args[0]

    @property
    def label(self) -> str:
        return str(self._label)

    def set_label(self, label: str) -> None:
        self._label = label

    @property
    def index_mask(self) -> tuple[int, ...]:
        """Sorted tuple of the hidden 1-indexed positions."""
        return tuple(self.args[1])

    def complement(self) -> MaskedGeneratingSet:
        """The masked set on the complementary positions."""
        return MaskedGeneratingSet(self.base_genset, [i for i in range(1, len(self.base_genset)) if i not in set(self.index_mask)])

    def __call__(self, index: int) -> Expr:
        """1-indexed"""
        return self[index - 1]

    @overload
    def __getitem__(self, index: int) -> Expr: ...

    @overload
    def __getitem__(self, index: slice) -> list[Expr]: ...

    def __getitem__(self, index):
        if isinstance(index, slice):
            start = index.start if index.start is not None else 0
            stop = index.stop if index.stop is not None else len(self)
            return [self[ii] for ii in range(start, stop)]
        return self.base_genset[self._mask[index]]

    def __iter__(self) -> Iterator[Expr]:
        yield from [self[i] for i in range(len(self))]

    def index(self, v: object) -> int:
        try:
            return self._index_lookup.get(v, self._index_lookup.get(sympify(v), -1))
        except SympifyError:
            return -1
        except TypeError:
            return -1

    def __hash__(self) -> int:
        return hash((self.base_genset, self.index_mask))

    def __len__(self) -> int:
        return len(self.base_genset) - len(self.index_mask)

    def __eq__(self, other: object) -> bool:
        return type(self) is type(other) and other.base_genset == self.base_genset and other.index_mask == self.index_mask  # type: ignore[attr-defined]


class CustomGeneratingSet(GeneratingSet_base):
    """A generating set over an explicit sequence of expressions (sympified on construction)."""

    _symbols_arr: tuple[Expr, ...]
    _index_lookup: dict[Any, int]

    def __new__(cls, gens: Iterable[Expr]) -> CustomGeneratingSet:
        return CustomGeneratingSet.__xnew_cached__(cls, tuple(gens))  # type: ignore[arg-type]  # mypy: type[...] vs Hashable (python/mypy#11470)

    @staticmethod
    @cache
    def __xnew_cached__(_class, gens):
        return CustomGeneratingSet.__xnew__(_class, gens)

    @staticmethod
    def __xnew__(_class, gens):
        obj = GeneratingSet_base.__new__(_class, gens)
        obj._symbols_arr = tuple([sympify(gens[i]) for i in range(len(gens))])
        obj._index_lookup = {obj._symbols_arr[i]: i for i in range(len(obj._symbols_arr))}
        return obj

    @overload
    def __getitem__(self, index: int) -> Expr: ...

    @overload
    def __getitem__(self, index: slice) -> tuple[Expr, ...]: ...

    def __getitem__(self, index):
        return self._symbols_arr[index]

    def __iter__(self) -> Iterator[Expr]:
        yield from self._symbols_arr

    def index(self, v: object) -> int:
        try:
            return self._index_lookup.get(v, self._index_lookup.get(sympify(v), -1))
        except SympifyError:
            return -1
        except TypeError:
            return -1

    def __call__(self, index: int) -> Expr:
        """1-indexed"""
        return self[index - 1]

    def __len__(self) -> int:
        return len(self.args[0])

    def __hash__(self) -> int:
        return hash(self._symbols_arr)

    def __eq__(self, other: object) -> bool:
        return type(self) is type(other) and other._symbols_arr == self._symbols_arr  # type: ignore[attr-defined]


NoneVar = 1e10
ZeroVar = 0


class NotEnoughGeneratorsError(ValueError):
    """Raised when an operation needs more generators than a generating set provides."""


@cache
def poly_genset(v: str | int | float) -> GeneratingSet_base:
    """``GeneratingSet(v)``, or a `ZeroGeneratingSet` for the sentinels ``ZeroVar``/``NoneVar``."""
    if v == ZeroVar:
        return ZeroGeneratingSet(tuple([sympify(0) for i in range(DEF_GENSET_SIZE)]))
    if v == NoneVar:
        return ZeroGeneratingSet(tuple([sympify(0) for i in range(DEF_GENSET_SIZE)]))
    return GeneratingSet(str(v))


def genset_dict_from_expr(expr: Expr, genset: GeneratingSet_base, length: int | None = None) -> dict[tuple[int, ...], Expr]:
    """Write a polynomial in the generators of ``genset`` as ``{exponent_tuple: coeff}``.

    Exponent tuples are 0-indexed by generator position ``genset(i) -> tuple[i - 1]`` and have
    length ``length`` (default: the largest generator index present). Factors free of the
    generators go into the coefficient; a factor mixing generators with other symbols raises.
    """
    if length is not None:
        k = length
    else:
        try:
            k = max([genset.index(a) for a in expr.free_symbols])
        except Exception:
            return {(): expr}
    poly: dict[tuple[int, ...], Expr] = {}
    expr = expand(expr)
    for term in Add.make_args(expr):
        coeff, exps = [], [0] * k

        for factor in Mul.make_args(term):
            if factor.is_Number:
                coeff.append(factor)
            else:
                try:
                    if isinstance(factor, Pow):
                        base, exp = factor.args[0], int(factor.args[1])
                        if base not in genset:
                            raise IndexError
                        exps[genset.index(base) - 1] = exp
                    else:
                        if factor not in genset:
                            raise IndexError
                        exps[genset.index(factor) - 1] = 1
                except IndexError:
                    if not any(a in factor.free_symbols for a in genset[:k]):
                        coeff.append(factor)
                    else:
                        raise Exception(f"{factor} contains an element of the set of generators.")

        monom = tuple(exps)

        if monom in poly:
            poly[monom] += Mul(*coeff)
        else:
            poly[monom] = Mul(*coeff)

    return poly
