# class generators with base
# symbols cls argument!

import re
from bisect import bisect_left
from functools import cache
from typing import ClassVar

from symengine import SympifyError, symbols, sympify
from sympy import Basic, Tuple
from sympy.core.symbol import Str

from schubmult.utils.logging import get_logger

logger = get_logger(__name__)


class GeneratingSet_base(Basic):
    def __new__(cls, *args):
        return Basic.__new__(cls, *args)

    def __getitem__(self, i):
        return NotImplemented

    def __len__(self):
        return NotImplemented


# variable registry
# TODO: ensure sympifies
# TODO: masked generating set
class GeneratingSet(GeneratingSet_base):
    def __new__(cls, name):
        return GeneratingSet.__xnew_cached__(cls, name)

    _registry: ClassVar = {}

    _index_pattern = re.compile("^([^_]+)_([0-9]+)$")
    _sage_index_pattern = re.compile("^([^0-9]+)([0-9]+)$")

    # is_Atom = True
    # TODO: masked generating set
    @staticmethod
    @cache
    def __xnew_cached__(_class, name):
        return GeneratingSet.__xnew__(_class, Str(str(name)))

    @staticmethod
    def __xnew__(_class, name):
        obj = GeneratingSet_base.__new__(_class, name)
        obj._symbols_arr = tuple([symbols(f"{name}_{i}") for i in range(100)])
        obj._index_lookup = {obj._symbols_arr[i]: i for i in range(len(obj._symbols_arr))}
        return obj

    def __call__(self, index):
        """1-indexed"""
        return self[index - 1]

    @property
    def label(self):
        return str(self.args[0])

    # index of v in the genset
    def index(self, v):
        try:
            return self._index_lookup.get(v, self._index_lookup.get(sympify(v), -1))
        except SympifyError:
            return -1

    def __repr__(self):
        return f"GeneratingSet('{self.label}')"

    def __str__(self):
        return self.name

    def _latex(self, printer):
        return printer.doprint(self.label)

    def _sympystr(self, printer):
        return printer.doprint(self.label)

    def __getitem__(self, i):
        return self._symbols_arr[i]

    def __len__(self):
        return len(self._symbols_arr)

    def __hash__(self):
        return hash(self.label)

    def __iter__(self):
        yield from [self[i] for i in range(len(self))]

    def __eq__(self, other):
        return isinstance(other, GeneratingSet) and self.label == other.label


class MaskedGeneratingSet(GeneratingSet_base):
    def __new__(cls, gset, index_mask):
        return MaskedGeneratingSet.__xnew_cached__(cls, gset, tuple(sorted(index_mask)))

    @staticmethod
    @cache
    def __xnew_cached__(_class, gset, index_mask):
        return MaskedGeneratingSet.__xnew__(_class, gset, index_mask)

    @staticmethod
    def __xnew__(_class, gset, index_mask):
        obj = GeneratingSet_base.__new__(_class, gset, Tuple(*index_mask))
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
        return obj

    @property
    def label(self):
        return str(self._label)

    def set_label(self, label):
        self._label = label

    @property
    def index_mask(self):
        return tuple(self.args[1])

    def complement(self):
        return MaskedGeneratingSet(self.base_genset, [i for i in range(1, len(self.base_genset)) if i not in set(self.index_mask)])

    @property
    def base_genset(self):
        return self.args[0]

    def __call__(self, index):
        """1-indexed"""
        return self[index - 1]

    def __getitem__(self, index):
        if isinstance(index, slice):
            start = index.start if index.start is not None else 0
            stop = index.stop if index.stop is not None else len(self)
            return [self[ii] for ii in range(start, stop)]
        return self.base_genset[self._mask[index]]

    def __iter__(self):
        yield from [self[i] for i in range(len(self))]

    def index(self, v):
        try:
            return self._index_lookup.get(v, self._index_lookup.get(sympify(v), -1))
        except SympifyError:
            return -1

    def __len__(self):
        return len(self.base_genset) - len(self.index_mask)


class CustomGeneratingSet(GeneratingSet_base):
    def __new__(cls, gens):
        return CustomGeneratingSet.__xnew_cached__(cls, tuple(gens))

    @staticmethod
    @cache
    def __xnew_cached__(_class, gens):
        return CustomGeneratingSet.__xnew__(_class, gens)

    @staticmethod
    def __xnew__(_class, gens):
        obj = GeneratingSet_base.__new__(_class, Tuple(*gens))
        obj._symbols_arr = [sympify(gens[i]) for i in range(len(gens))]
        obj._index_lookup = {obj._symbols_arr[i]: i for i in range(len(obj._symbols_arr))}
        return obj

    def __getitem__(self, index):
        return self._symbols_arr[index]

    def __iter__(self):
        yield from [self[i] for i in range(len(self))]

    def index(self, v):
        try:
            return self._index_lookup.get(v, self._index_lookup.get(sympify(v), -1))
        except SympifyError:
            return -1

    def __call__(self, index):
        """1-indexed"""
        return self[index - 1]

    def __len__(self):
        return len(self.args[0])
