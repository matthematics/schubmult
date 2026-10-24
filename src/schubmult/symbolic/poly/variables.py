"""Generating sets: indexed families of variables ``x_0, x_1, x_2, ...`` used as ring generators.

`GeneratingSet("x")` starts with ``DEF_GENSET_SIZE`` symbols and doubles its cache as needed;
nonnegative indexing and bounded slices have no fixed upper limit. Length, iteration, negative
indices, and open-ended slices refer to the currently cached finite prefix. ``gs[i]`` is the symbol
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
from bisect import bisect_left, bisect_right
from collections.abc import Callable, Iterable, Iterator, Sequence
from functools import cache
from typing import TYPE_CHECKING, Any, ClassVar, TypeVar, overload

from schubmult.symbolic import Add, Mul, Pow, S, SympifyError, efficient_subs, expand, symbols, sympify
from schubmult.utils.logging import get_logger

if TYPE_CHECKING:
    from schubmult._typing import Expr

logger = get_logger(__name__)

DEF_GENSET_SIZE = 10

_GeneratingSetT = TypeVar("_GeneratingSetT", bound="GeneratingSet_base")


class GeneratingSet_base:
    """Interface for generating sets: indexing, length, ``index(symbol)`` (``-1`` if absent), and ``label``."""

    _args: tuple[Any, ...]
    is_infinite = False

    def __new__(cls: type[_GeneratingSetT], *args: Any) -> _GeneratingSetT:
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

    is_infinite = True

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


# TODO: ensure sympifies
class GeneratingSet(GeneratingSet_base):
    """The interned family ``name_0, name_1, ...``; ``gs[i]`` is the symbol ``name_i`` and ``gs(i)`` is ``gs[i - 1]``."""

    _symbols_arr: tuple[Expr, ...]
    _index_lookup: dict[Any, int]
    _hash: int
    is_infinite = True

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
            return self._index_lookup[v]
        except (KeyError, TypeError):
            pass
        try:
            v = sympify(v)
            index = self._index_lookup.get(v)
            if index is not None:
                return index
            if getattr(v, "is_Symbol", False):
                name, separator, suffix = str(v).rpartition("_")
                if separator and name == self.label and suffix.isascii() and suffix.isdecimal():
                    index = int(suffix)
                    if suffix == str(index) and self[index] == v:
                        return index
        except (SympifyError, TypeError):
            return -1
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
        if i.__class__ is slice:
            if (i.step is None or i.step > 0) and i.stop is not None and i.stop > len(self._symbols_arr):
                self._grow(i.stop)
            elif i.step is not None and i.step < 0 and i.start is not None and i.start >= len(self._symbols_arr):
                self._grow(i.start + 1)
            return self._symbols_arr[i]
        try:
            return self._symbols_arr[i]
        except IndexError:
            if i < 0:
                raise
            self._grow(i + 1)
            return self._symbols_arr[i]

    def _grow(self, size: int) -> None:
        old_size = len(self._symbols_arr)
        new_size = max(size, 2 * old_size)
        name = self.label
        new_symbols = tuple(symbols(f"{name}_{i}") for i in range(old_size, new_size))
        self._symbols_arr += new_symbols
        self._index_lookup.update((symbol, i) for i, symbol in enumerate(new_symbols, old_size))

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
    _mask_offsets: tuple[int, ...]

    def __new__(cls, gset: GeneratingSet_base, index_mask: Iterable[int]) -> MaskedGeneratingSet:
        return MaskedGeneratingSet.__xnew_cached__(cls, gset, tuple(sorted(set(index_mask))))  # type: ignore[arg-type]  # mypy: type[...] vs Hashable (python/mypy#11470)

    @staticmethod
    @cache
    def __xnew_cached__(_class, gset, index_mask):
        return MaskedGeneratingSet.__xnew__(_class, gset, index_mask)

    @staticmethod
    def __xnew__(_class, gset, index_mask):
        obj = GeneratingSet_base.__new__(_class, gset, index_mask)
        if any(i <= 0 for i in index_mask):
            raise ValueError("masked indices must be positive")
        if not gset.is_infinite and index_mask and index_mask[-1] >= len(gset):
            raise IndexError("masked index out of range")
        mask_dict = {0: 0} if len(gset) else {}
        for i in range(1, len(gset)):
            index = bisect_left(index_mask, i)
            if index >= len(index_mask) or index_mask[index] != i:
                mask_dict[i - index] = i
        obj._mask = mask_dict
        obj._mask_offsets = tuple(index - i for i, index in enumerate(index_mask))
        obj._index_lookup = {gset[base_index]: i for i, base_index in mask_dict.items()}
        obj._label = gset.label
        obj.is_infinite = gset.is_infinite
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
        if self.is_infinite:
            return MaskedGeneratingSet(CustomGeneratingSet([self.base_genset[0], *(self.base_genset[i] for i in self.index_mask)]), ())
        return MaskedGeneratingSet(self.base_genset, [i for i in range(1, len(self.base_genset)) if i not in set(self.index_mask)])

    def __call__(self, index: int) -> Expr:
        """1-indexed"""
        return self[index - 1]

    @overload
    def __getitem__(self, index: int) -> Expr: ...

    @overload
    def __getitem__(self, index: slice) -> list[Expr]: ...

    def __getitem__(self, index):
        if index.__class__ is slice:
            stop = index.stop
            size = len(self)
            if self.is_infinite and stop is not None and stop > size and (index.step is None or index.step > 0):
                size = stop
            elif self.is_infinite and index.step is not None and index.step < 0 and index.start is not None and index.start >= size:
                size = index.start + 1
            return [self[ii] for ii in range(*index.indices(size))]
        try:
            return self.base_genset[self._mask[index]]
        except KeyError:
            if index < 0:
                index += len(self)
            if index < 0 or (not self.is_infinite and index >= len(self)):
                raise IndexError("generating-set index out of range") from None
            base_index = index + bisect_right(self._mask_offsets, index)
            self._mask[index] = base_index
        return self.base_genset[base_index]

    def __iter__(self) -> Iterator[Expr]:
        yield from [self[i] for i in range(len(self))]

    def index(self, v: object) -> int:
        try:
            return self._index_lookup[v]
        except (KeyError, TypeError):
            pass
        index = self.base_genset.index(v)
        if index < 0:
            return -1
        offset = bisect_left(self.index_mask, index)
        if offset < len(self.index_mask) and self.index_mask[offset] == index:
            return -1
        position = index - offset
        self._index_lookup[self.base_genset[index]] = position
        return position

    def __hash__(self) -> int:
        return hash((self.base_genset, self.index_mask))

    def __len__(self) -> int:
        size = len(self.base_genset)
        return size - bisect_left(self.index_mask, size)

    def __eq__(self, other: object) -> bool:
        return type(self) is type(other) and other.base_genset == self.base_genset and other.index_mask == self.index_mask


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
        return len(self._symbols_arr)

    def __hash__(self) -> int:
        return hash(self._symbols_arr)

    def __eq__(self, other: object) -> bool:
        return type(self) is type(other) and other._symbols_arr == self._symbols_arr


class ChangeOfVariables(GeneratingSet_base):
    """An indexed change of variables from ``source`` to expressions in ``target``.

    ``rule(i, target)`` gives the image of the zero-indexed generator ``source[i]``.
    Use ``subs(expr)`` to apply the resulting substitution to an expression.
    """

    _source: GeneratingSet_base
    _target: GeneratingSet_base
    _rule: Callable[[int, GeneratingSet_base], Expr]
    _images: dict[int, Expr]

    def __new__(
        cls,
        source: GeneratingSet_base,
        target: GeneratingSet_base,
        rule: Callable[[int, GeneratingSet_base], Expr],
    ) -> ChangeOfVariables:
        obj = GeneratingSet_base.__new__(cls, source, target, rule)
        obj._source = source
        obj._target = target
        obj._rule = rule
        obj._images = {}
        obj.is_infinite = source.is_infinite and target.is_infinite
        return obj

    @property
    def source(self) -> GeneratingSet_base:
        """The generating set being replaced."""
        return self._source

    @property
    def target(self) -> GeneratingSet_base:
        """The generating set used to form replacement expressions."""
        return self._target

    @property
    def label(self) -> str | None:
        """The label of the target generating set."""
        return self._target.label

    def __getitem__(self, index):
        if isinstance(index, slice):
            size = len(self)
            if self.is_infinite and index.stop is not None and index.stop > size and (index.step is None or index.step > 0):
                size = index.stop
            elif self.is_infinite and index.step is not None and index.step < 0 and index.start is not None and index.start >= size:
                size = index.start + 1
            return tuple(self[i] for i in range(*index.indices(size)))
        if index < 0:
            index += len(self)
        if index < 0 or (not self.is_infinite and index >= len(self)):
            raise IndexError("generating-set index out of range")
        if index not in self._images:
            self._images[index] = sympify(self._rule(index, self._target))
        return self._images[index]

    def __call__(self, index: int) -> Expr:
        """Return the image of the 1-indexed source generator."""
        return self[index - 1]

    def __len__(self) -> int:
        if self._target.is_infinite:
            return len(self._source)
        if self._source.is_infinite:
            return len(self._target)
        return min(len(self._source), len(self._target))

    def index(self, value: object) -> int:
        """Position of a source generator, or ``-1`` if it is not in the source."""
        return self._source.index(value)

    def subs(self, expr: Expr) -> Expr:
        """Apply this change of variables to ``expr``."""
        expr = sympify(expr)
        substitutions = {}
        for symbol in expr.free_symbols:
            index = self._source.index(symbol)
            if index >= 0:
                substitutions[symbol] = self[index]
        return efficient_subs(expr, substitutions)

    def __repr__(self) -> str:
        return f"ChangeOfVariables({self._source!r}, {self._target!r})"


NoneVar = 1e10
ZeroVar = 0


class NotEnoughGeneratorsError(ValueError):
    """Raised when an operation needs more generators than a generating set provides."""


def _generator_prefix(gens: Iterable[Expr], length: int) -> tuple[Expr, ...]:
    """Materialize only the required prefix when given an indexed generating set."""
    if isinstance(gens, GeneratingSet_base):
        return tuple(gens[:max(0, length)])
    return tuple(gens)


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
