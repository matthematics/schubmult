import logging  # noqa: F401
from functools import cache, cached_property

from sympy.printing.defaults import Printable

import schubmult.schub_lib.schub_lib as schub_lib
from schubmult.perm_lib import Permutation, uncode
from schubmult.rings import ASx
from schubmult.symbolic import S, prod
from schubmult.utils.logging import get_logger, init_logging
from schubmult.utils.perm_utils import add_perm_dict

from .free_algebra import FreeAlgebra, FreeAlgebraElement
from .free_algebra_basis import WordBasis
from .nil_hecke import NilHeckeRing

init_logging(debug=False)
logger = get_logger(__name__)
# logging.basicConfig(level=logging.DEBUG,
#                      format='%(asctime)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s',)


def _is_row_root(row, root):
    return root[0] <= row and root[1] > row


FA = FreeAlgebra(WordBasis)
FAS = FA


def debug_print(*args, debug=False):
    if debug:
        print(*args)  # noqa: T201

class _SparseTupleList:

    def __new__(cls, rows=None, logical_length=None):
        obj = object.__new__(cls)
        obj._rows = {i: row for i, row in enumerate(rows) if len(row) > 0} if rows is not None else {}
        obj._length = logical_length if logical_length is not None else len(rows) if rows is not None else 0
        return obj

    def __init__(self, rows, logical_length=None):
        pass

    def _getitem(self, idx):
        if isinstance(idx, int):
            if idx >= self._length:
                raise IndexError("Index out of range")
            if idx < 0:
                idx += self._length
            return self._rows.get(idx, ())
        elif isinstance(idx, slice):
            return tuple(self[i] for i in range(*idx.indices(self._length)))
        else:
            raise TypeError("Invalid index type")

    def __getitem__(self, idx):
        raise NotImplementedError

    def __len__(self):
        return self._length

    def __iter__(self):
        for i in range(self._length):
            yield self._rows.get(i, ())


class RCGraph(Printable, _SparseTupleList):
    def __eq__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return self._rows == other._rows and self._length == other._length

    @property
    def is_valid(self):
        if self.perm.inv != len(self.perm_word):
            return False
        return True

    def shiftup(self, shift):
        return [tuple([a + shift for a in rrow]) for rrow in self]

    @cache
    def right_root_at(self, i, j, debug=False):  # noqa: ARG002
        if i <= 0 or j <= 0:
            raise IndexError("i and j must be positive")
        if len(self.perm_word) > 0:
            index = self.bisect_left_coords_index(i, j)
            if index < len(self.perm_word):
                if self.left_to_right_inversion_coords(index) == (i, j):
                    return self.perm.right_root_at(index, word=self.perm_word)
                word_piece = list(self.perm_word[index:])
            else:
                word_piece = []
            refl = ~Permutation.ref_product(*word_piece)
            result = refl.act_root(i + j - 1, i + j)

        else:
            result = (i + j - 1, i + j)
        return result

    @cache
    def left_root_at(self, i, j):
        start_root = (i + j - 1, i + j)
        for j2 in range(j + 1, self.cols):
            if self[i - 1, j2 - 1]:
                start_root = Permutation.ref_product(self[i - 1, j2 - 1]).act_root(*start_root)
        for i2 in range(i - 1, 0, -1):
            for j2 in range(1, self.cols + 1):
                if self[i2 - 1, j2 - 1]:
                    start_root = Permutation.ref_product(self[i2 - 1, j2 - 1]).act_root(*start_root)
        return start_root

    def reverse_kogan_kumar_insert(self, descent, reflection_path, return_rows=False, debug=False):
        # pair_sequence = sorted(pair_sequence, key=lambda x: x[0]
        pair_dict = {}
        for ref in reflection_path:
            a, b = ref
            if a not in pair_dict:
                pair_dict[a] = set()
            pair_dict[a].add(b)
        pair_dict_rev = {}
        # ref_by_index = {}
        for a, b_list in pair_dict.items():
            for b in b_list:
                assert a <= descent and descent < b  # noqa: PT018
                pair_dict_rev[b] = a

        def is_relevant_crossing(root, prm):  # noqa: ARG001
            # min_root = max(pair_dict.keys())

            # if root[0] > descent:
            #     if root[0] not in pair_dict_rev or pair_dict_rev.get(root[0], 0) != pair_dict_rev.get(root[1], 0):
            #         return False
            #     return True
            if root[0] not in pair_dict:
                if root[0] in pair_dict_rev and root[1] in pair_dict_rev and pair_dict_rev[root[0]] == pair_dict_rev[root[1]]:
                    return True
                return False
            return root[0] in pair_dict and root[1] in pair_dict[root[0]]

        # may have to add q, s or a_i, q
        def is_relevant_noncrossing(root):
            top, bottom = max(root[1], root[0]), min(root[0], root[1])
            return (bottom <= descent and descent < top and top not in pair_dict_rev) or (bottom in pair_dict_rev and top > descent and top not in pair_dict_rev)

        # Add this intersection. If we are in the first case, insert (s, q) into the sequence (ai, bi) in the rightmost position, such that ai’s remain nondecreasing in the # noqa: RUF003
        # sequence. ((s, q) are the rows where the two strands shown in Figure 3 originate.) If
        # we are in the second case, add (ai, q) just before where (a, bi) is in the sequence.

        working_rc = self
        if debug:
            debug_print("Starting with", debug=debug)
            debug_print(working_rc)
            debug_print(working_rc.perm)
        rows = []
        # read_rows = list(range(1, len(self) + 1))
        for row in range(1, len(self) + 1):
            for col in range(1, max(descent, self.cols) + 1):
                if working_rc.has_element(row, col):
                    a, b = working_rc.right_root_at(row, col)
                    if debug:
                        debug_print(f"{a,b=} at {row,col=}", debug=debug)
                        debug_print(working_rc, debug=debug)
                    # print(f"{reflection_path=}")
                    # print(f"{(a,b)=}")
                    if is_relevant_crossing((a, b), working_rc.perm):
                        working_rc = working_rc.toggle_ref_at(row, col)
                        if debug:
                            debug_print(f"Toggled {a,b=}", debug=debug)
                            debug_print(working_rc, debug=debug)
                        a2 = a
                        if a2 in pair_dict_rev:
                            a2 = pair_dict_rev[a2]
                        #     pair_dict[a2].remove(a)
                        #     del pair_dict_rev[a]
                        # else:
                        pair_dict[a2].remove(b)
                        del pair_dict_rev[b]
                        if len(pair_dict[a2]) == 0:
                            del pair_dict[a2]
                        # pair_dict_rev.pop(a, None)
                        rows.append(row)
                        # if len(pair_dict_rev) == 0:
                        #     break
                        for col2 in range(1, col):
                            if not working_rc.has_element(row, col):
                                a2, b2 = working_rc.right_root_at(row, col2)
                                # if a2 > b2:
                                #     continue
                                if is_relevant_noncrossing((a2, b2)):
                                    # print(f"Rect {a2, b2}")
                                    if a2 <= descent:
                                        assert b2 not in pair_dict
                                        if a2 not in pair_dict:
                                            pair_dict[a2] = set()
                                        pair_dict[a2].add(b2)
                                        pair_dict_rev[b2] = a2
                                        working_rc = working_rc.toggle_ref_at(row, col2)
                                        rows.pop()
                                    else:
                                        assert a2 not in pair_dict_rev, f"{pair_dict_rev=}"
                                        assert b2 in pair_dict_rev, f"{pair_dict_rev=}"
                                        a = pair_dict_rev[a2]
                                        pair_dict[a].add(b2)
                                        pair_dict_rev[b2] = a
                                        working_rc = working_rc.toggle_ref_at(row, col2)
                                        rows.pop()
                                    break

        assert len(pair_dict_rev) == 0, f"{pair_dict=}, {pair_dict_rev=}, {working_rc=}"
        assert working_rc.perm.bruhat_leq(self.perm)
        if return_rows:
            return working_rc, rows
        return working_rc

    @cache
    def inversion_label(self, i, j):
        if i >= j:
            raise ValueError("i must be less than j")
        if self.perm[i] < self.perm[j]:
            raise ValueError("Not an inversion")
        for index in len(self.perm_word):
            if self.left_to_right_inversion(index) == (i + 1, j + 1):
                return self.left_to_right_inversion_coords(index)[0]
        raise ValueError("Could not find inversion")

    @cache
    def lehmer_label(self, i, j):
        value = self.inversion_label(i, j)
        numeros = set(range(1, value + 1))
        for ip in range(i):
            try:
                numeros.remove(self.inversion_label(ip, j))
            except ValueError:
                pass
            except KeyError:
                pass
        return len(numeros)

    def __new__(cls, *args):
        new_args = tuple(tuple(arg) for arg in args)
        return RCGraph.__xnew_cached__(cls, *new_args)

    @staticmethod
    @cache
    def __xnew_cached__(_class, *args):
        return RCGraph.__xnew__(_class, *args)

    @staticmethod
    def __xnew__(_class, *args, **kwargs):
        obj = _SparseTupleList.__new__(_class, *args, **kwargs)
        obj._hashcode = hash(tuple(obj))
        return obj

    def __init__(self, *args, **kwargs):
        pass

    @cached_property
    def perm_word(self):
        ret = []
        for row in self:
            ret = [*ret, *row]
        return tuple(ret)

    # def edelman_greene(self):
    #     from schubmult.perm_lib import NilPlactic

    #     word1 = []
    #     word2 = []
    #     index = 0
    #     evil_self = list(reversed([list(reversed(row)) for row in self]))
    #     for i in range(len(evil_self)):
    #         for a in evil_self[i]:
    #             to_insert = len(self) - i
    #             word1, word2 = NilPlactic.ed_insert_rsk(word1, word2, a, to_insert)
    #             index += 1
    #     P = Tableau(word1)
    #     Q = Tableau(word2)
    #     # reg._rc_graph = self
    #     return (P, Q)

    # def __matmul__(self, other):
    #     if isinstance(other, RCGraph):
    #         return RCGraphModule({self: 1}, generic_key_type=self.__class__) @ RCGraphModule({other: 1}, generic_key_type=self.__class__)
    #     if isinstance(other, ModuleType):
    #         return RCGraphModule({self: 1}, generic_key_type=self.__class__) @ other
    #     return NotImplemented

    # def __rmul__(self, other):
    #     if isinstance(other, RCGraphModule):
    #         return NotImplemented
    #     return other * RCGraphModule({self: 1}, generic_key_type=self.__class__)

    # def __mul__(self, other):
    #     if isinstance(other, RCGraph):
    #         return self.prod_with_rc(other)
    #     return NotImplemented

    def asdtype(self, cls):
        return cls.dtype().ring.from_rc_graph(self)

    def as_nil_hecke(self, x, y=None):
        R = NilHeckeRing(x)
        return self.polyvalue(x, y) * R(self.perm)

    @cache
    def has_element(self, i, j):
        return i <= len(self) and j + i - 1 in self[i - 1]

    @cached_property
    def length_vector(self):
        return tuple([len(row) for row in self])

    @cache
    def lehmer_partial_leq(self, other):
        try:
            for i in range(self.perm.inv):
                a, b = self.perm.right_root_at(i)
                if self.lehmer_label(a - 1, b - 1) > other.lehmer_label(a - 1, b - 1):
                    return False
        except ValueError:
            return False
        return True

    def rowrange(self, start, end):
        if start == end:
            return RCGraph(())
        return RCGraph([tuple([a - start for a in row]) for row in self[start:end]])

    def polyvalue(self, x, y=None):
        ret = S.One
        for i, row in enumerate(self):
            if y is None:
                ret *= x[i + 1] ** len(row)
            else:
                ret *= prod([x[i + 1] - y[row[j] - i] for j in range(len(row))])
        return ret

    _graph_cache = {}  # noqa: RUF012

    _cache_by_weight = {}  # noqa: RUF012

    @classmethod
    def all_rc_graphs(cls, perm, length=-1, weight=None):
        if length > 0 and length < len(perm.trimcode):
            raise ValueError(f"Length must be at least the last descent of the permutation, permutation has {len(perm.trimcode)} rows and {perm=}, got {length=}")
        if length < 0:
            length = len(perm.trimcode)
        if weight and len(weight) != length:
            raise ValueError("Weight must have length equal to the number of rows")
        if weight:
            if (perm, tuple(weight)) in cls._cache_by_weight:
                return cls._cache_by_weight[(perm, tuple(weight))]
        elif (perm, length) in cls._graph_cache:
            return cls._graph_cache[(perm, length)]
        if perm.inv == 0:
            return {RCGraph([()] * length if length > 0 else [])}
        if len(perm.trimcode) == 1:
            nrc = RCGraph((tuple(range(perm.code[0], 0, -1)),))
            if len(nrc) < length:
                nrc = RCGraph((*nrc, *tuple([()] * (length - len(nrc)))))
                assert len(nrc) == length
            return {nrc}
        ret = set()
        pm = perm
        L = schub_lib.pull_out_var(1, pm)
        for _, new_perm in L:
            new_row = [new_perm[i] for i in range(max(len(pm), len(new_perm))) if new_perm[i] == pm[i + 1]]
            if weight and len(new_row) != weight[0]:
                continue
            new_row.sort(reverse=True)
            if weight:
                oldset = cls.all_rc_graphs(new_perm, length=length - 1, weight=weight[1:])
            else:
                oldset = cls.all_rc_graphs(new_perm, length=length - 1)
            for old_rc in oldset:
                nrc = RCGraph([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in old_rc]])
                assert nrc.perm == perm
                assert len(nrc) == length
                ret.add(nrc)
        if weight:
            cls._cache_by_weight[(perm, tuple(weight))] = ret
        else:
            cls._graph_cache[(perm, length)] = ret
        return ret

    def extend(self, extra_rows):
        return RCGraph([*self, *tuple([()] * extra_rows)])

    def _kogan_kumar_insert_row(self, row, descent, dict_by_a, dict_by_b, num_times, debug=False, start_index=-1, backwards=True):
        working_rc = self
        if row > descent:
            raise ValueError("All rows must be less than or equal to descent")
        if backwards:
            debug_print("WARNING BACKWARDS IS ON", debug=False)

        i = start_index

        if i == -1:
            if backwards:
                i = 0
            elif len(self[row - 1]) == 0:
                i = descent + 5
            else:
                i = max(self[row - 1]) + descent + 5
        num_done = 0
        flag = True
        while num_done < num_times:
            if i <= 1 and not backwards:
                debug_print("Resetting i", debug=debug)
                i = working_rc.cols + descent + 5
            if not backwards:
                i -= 1
            else:
                i += 1
            flag = False
            if debug:
                debug_print(f"Trying column {i=} {descent=} {row=} {num_done=} {num_times=} {dict_by_a=} {dict_by_b=}", debug=debug)
                debug_print(f"{working_rc=}", debug=debug)
            if not working_rc.has_element(row, i):
                a, b = working_rc.right_root_at(row, i)
                if a < b:
                    if debug:
                        debug_print(f"root is {a, b}", debug=debug)
                    flag = False
                    if debug:
                        debug_print(f"{dict_by_b=}", debug=debug)
                    if _is_row_root(descent, (a, b)) and b not in dict_by_b:
                        working_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[a] = dict_by_a.get(a, set())
                        dict_by_a[a].add(b)
                        dict_by_b[b] = a
                        flag = True

                        if debug:
                            debug_print("Toggled a", debug=debug)
                            debug_print(working_rc, debug=debug)
                    elif a in dict_by_b and b > descent and b not in dict_by_b:
                        working_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[dict_by_b[a]].add(b)
                        dict_by_b[b] = dict_by_b[a]
                        flag = True

                        if debug:
                            debug_print("Toggled b", debug=debug)
                            debug_print(working_rc, debug=debug)
                if flag:
                    num_done += 1
                    debug_print("did it", debug=debug)
                    debug_print(f"{working_rc=}", debug=debug)
                if row > 1 and not working_rc.is_valid:
                    working_rc = working_rc._kogan_kumar_rectify(row - 1, descent, dict_by_a, dict_by_b, backwards=backwards)  # minus one?
        return working_rc

    def _kogan_kumar_rectify(self, row_below, descent, dict_by_a, dict_by_b, debug=False, backwards=True):
        debug_print(f"In rectify {row_below=} {self.is_valid=} {self=}", debug=debug)
        working_rc = self
        if row_below == 0:
            assert working_rc.is_valid, f"{working_rc=}, {dict_by_a=}, {dict_by_b=}"
            return working_rc
        if working_rc.is_valid:
            return working_rc
        for j in range(1, working_rc.cols + descent + 5):
            flag = False
            if working_rc.is_valid:
                return working_rc
            debug_print(f"Considering column {j=} {row_below=} {descent=}", debug=debug)
            if working_rc.has_element(row_below, j):
                debug_print("Has element", debug=debug)
                a, b = working_rc.right_root_at(row_below, j)

                top, bottom = max(a, b), min(a, b)
                debug_print("root=", (a, b), debug=debug)
                debug_print("Entered", debug=debug)
                if a < b:
                    debug_print("Root is positive", debug=debug)
                    continue
                if debug:
                    debug_print(f"Considering bad at {row_below, j}")
                    debug_print(f"{dict_by_a=}, {dict_by_b=}")
                    debug_print(f"root = ({a, b})")
                if bottom in dict_by_a and top in dict_by_a[bottom]:
                    new_rc = working_rc.toggle_ref_at(row_below, j)
                    dict_by_a[bottom].remove(top)
                    if len(dict_by_a[bottom]) == 0:
                        del dict_by_a[bottom]
                    del dict_by_b[top]
                    working_rc = new_rc
                    flag = True
                    debug_print("Toggle bad a", debug=debug)
                    debug_print(working_rc)
                elif bottom in dict_by_b and top in dict_by_b and dict_by_b[top] == dict_by_b[bottom]:
                    new_rc = working_rc.toggle_ref_at(row_below, j)
                    dict_by_a[dict_by_b[bottom]].remove(top)
                    del dict_by_b[top]
                    if len(dict_by_a[dict_by_b[top]]) == 0:
                        del dict_by_a[dict_by_b[top]]
                    flag = True
                    working_rc = new_rc
                    debug_print("Toggle bad c", debug=debug)
                else:
                    debug_print(f"{working_rc=}, {dict_by_a=}, {dict_by_b=}, {a=}, {b=} {row_below=} {j=}")
                    raise ValueError(f"Could not rectify at {(row_below, j)} with root {(a, b)}")
                if flag:
                    working_rc = working_rc._kogan_kumar_insert_row(row_below, descent, dict_by_a, dict_by_b, num_times=1, debug=debug, backwards=backwards)
        return working_rc._kogan_kumar_rectify(row_below - 1, descent, dict_by_a, dict_by_b, backwards=backwards)

    def associative_kogan_kumar_insert(self, descent, rows, debug=False):
        """WIP"""
        if len(self.perm.trimcode) <= descent:
            return self.kogan_kumar_insert(descent, rows, debug=debug)
        max_desc = len(self.perm.trimcode)
        diff_rows_stack = []
        desc_stack = []
        interim = RCGraph([*self])
        while max_desc > descent:
            diff_rows = []

            while interim.perm.inv > 0 and len(interim.perm.trimcode) > max_desc - 1:
                prev_interim = interim
                interim = interim.exchange_property(len(interim.perm.trimcode))
                for i in range(len(interim)):
                    if len(interim[i]) < len(prev_interim[i]):
                        rw = i + 1
                        diff_rows.append(rw)
                        break

            diff_rows_stack.append(diff_rows)
            desc_stack.append(max_desc)
            max_desc = len(interim.perm.trimcode)
        interim = interim.kogan_kumar_insert(descent, rows, debug=debug)
        for i in range(len(diff_rows_stack) - 1, -1, -1):
            diff_rows = diff_rows_stack[i]
            prev_descent = desc_stack[i]
            interim = interim.kogan_kumar_insert(prev_descent, diff_rows, debug=debug)
        return interim

    def associative_product(self, other, debug=False):
        """WIP"""
        max_desc = len(self.perm.trimcode)
        diff_rows_stack = []
        desc_stack = []
        interim = RCGraph([*self])
        while max_desc > 0:
            diff_rows = []

            while interim.perm.inv > 0 and len(interim.perm.trimcode) > max_desc - 1:
                prev_interim = interim
                interim = interim.exchange_property(len(interim.perm.trimcode))
                for i in range(len(interim)):
                    if len(interim[i]) < len(prev_interim[i]):
                        rw = i + 1
                        diff_rows.append(rw)
                        break

            diff_rows_stack.append(diff_rows)
            desc_stack.append(max_desc)
            max_desc = len(interim.perm.trimcode)
        for i in range(len(diff_rows_stack) - 1, -1, -1):
            diff_rows = diff_rows_stack[i]
            prev_descent = desc_stack[i]
            other = other.associative_kogan_kumar_insert(prev_descent, diff_rows, debug=debug)
        return other

    # VERIFY
    def kogan_kumar_insert(self, descent, rows, debug=False, return_reflections=False, backwards=True):
        dict_by_a = {}
        dict_by_b = {}
        # row is descent
        # inserting times

        working_rc = RCGraph([*self])
        if len(rows) == 0:
            if return_reflections:
                return working_rc, ()
            return self
        rows_grouping = {}

        for r in rows:
            rows_grouping[r] = rows_grouping.get(r, 0) + 1
        if max(rows) > len(working_rc):
            working_rc = working_rc.extend(max(rows) - len(working_rc))
        rows = sorted(rows, reverse=True)
        if debug:
            debug_print(f"inserting {rows=} into")
            debug_print(self)
        for row in sorted(rows_grouping.keys(), reverse=True):
            num_times = rows_grouping[row]
            if debug:
                debug_print(f"Inserting {row=} {num_times=}")
                debug_print(working_rc)
                debug_print(f"{working_rc.perm.inv=}, {self.perm.inv=}")
            last_working_rc = working_rc
            working_rc = working_rc._kogan_kumar_insert_row(row, descent, dict_by_a, dict_by_b, num_times, debug=debug, backwards=backwards)

            if row > 1 and not working_rc.is_valid:
                working_rc = working_rc._kogan_kumar_rectify(row - 1, descent, dict_by_a, dict_by_b, backwards=backwards)  # minus one?

            try:
                assert len(working_rc[row - 1]) == len(last_working_rc[row - 1]) + num_times
            except AssertionError:
                # print("Assertion failed")
                # print(working_rc)
                # print(f"{working_rc.perm.inv=}, {self.perm.inv + total_num=}")
                # print(f"{dict_by_a=}, {dict_by_b=}")
                # print(f"{working_rc.perm=}, {self.perm=}")
                #  if debug:
                #      raise
                #  self.kogan_kumar_insert(descent, rows, debug=False)
                raise
        if return_reflections:
            reflections = []
            for a in sorted(dict_by_a.keys()):
                reflections = [*reflections, *[(a, b) for b in dict_by_a[a]]]
            return working_rc, tuple(reflections)
        return working_rc

    @cached_property
    def perm(self):
        perm = Permutation([])
        for row in self:
            for p in row:
                perm = perm.swap(p - 1, p)
        return perm

    def transpose(self):
        newrc = []
        trimself = [list(row) for row in self]
        i = 0
        while len(newrc) < len(self) or any(len(row) > 0 for row in trimself):
            new_row = []
            for index in range(len(trimself)):
                if len(trimself[index]) > 0 and trimself[index][-1] == index + i + 1:
                    new_row += [index + i + 1]
                    trimself[index].pop()
            new_row.reverse()
            newrc.append(tuple(new_row))
            i += 1
        new_rc = RCGraph(newrc)

        assert new_rc.perm == ~self.perm
        return new_rc

    @classmethod
    def one_row(cls, p):
        return RCGraph((tuple(range(p, 0, -1)),))

    def weak_order_leq(self, other):
        for i in range(self.perm.inv):
            a, b = self.perm.right_root_at(i)
            try:
                if self.inversion_label(a - 1, b - 1) > other.inversion_label(a - 1, b - 1):
                    return False
            except ValueError:
                return False
        return True

    w_key_cache = {}  # noqa: RUF012
    rc_cache = set()  # noqa: RUF012

    @classmethod
    def fa_hom(cls, *word):
        if len(word) == 0:
            return cls()
        if len(word) == 1:
            return FA(*word) * cls()
        return cls.fa_hom(*word[:-1]) * cls.fa_hom(word[-1])

    def toggle_ref_at(self, i, j, debug=False):
        if i <= 0 or j <= 0:
            raise IndexError()
        new_row = [*self[i - 1]]
        if i + j - 1 in new_row:
            index = new_row.index(i + j - 1)
            new_row = [*new_row[:index], *new_row[index + 1 :]]
        else:
            index = 0
            debug_print("Searching", debug=debug)
            if len(new_row) > 0:
                while index < len(new_row) and new_row[index] > i + j - 1:
                    index += 1
            debug_print(f"Index {index=} inserting at row {i} column {j} the refl {i+j-1=} int {new_row=}", debug=debug)
            new_row.insert(index, i + j - 1)
            debug_print(f"Result {new_row=}", debug=debug)
        return RCGraph([*self[: i - 1], tuple(new_row), *self[i:]])

    # # THIS IS KEY
    # # EXCHANGE PROPERTY GOES TO UNIQUE PERMUTATION
    # # KOGAN INSERT ENSURES WE GO UP PROPERLY

    _z_cache = {}  # noqa: RUF012

    @cache
    def zero_out_last_row(self, debug=False):
        # this is important!
        # transition formula
        if len(self[-1]) != 0:
            raise ValueError("Last row not empty")
        if self.perm.inv == 0:
            return self.rowrange(0, len(self) - 1)
        interim = RCGraph([*self])

        if debug:
            debug_print("Zeroing out last row", debug=debug)
            debug_print(f"{self=} {len(self)=}", debug=debug)
            debug_print("-------", debug=debug)
        diff_rows = []
        descs = []
        extend_amount = 1

        while len(interim.perm.trimcode) > len(self) - 1:
            debug_print(f"Exchanging {interim.perm=} {len(interim.perm.trimcode)=}", debug=debug)
            descs += [len(interim.perm.trimcode)]
            interim, row = interim.exchange_property(len(interim.perm.trimcode), return_row=True)
            diff_rows += [row]
        debug_print(f"Exchanged to {interim=}, {diff_rows=} {descs=}", debug=debug)

        interim2 = RCGraph([*interim[:-1], tuple(sorted(descs, reverse=True))])
        interim = interim2.kogan_kumar_insert(len(self.perm.trimcode) - extend_amount, diff_rows, debug=False)
        if debug:
            debug_print("Got", debug=debug)
            debug_print(interim, debug=debug)
            debug_print(f"{interim.perm=}", debug=debug)

        # if interim.perm[len(self) - 1] != len(interim.perm):
        #     #raise ValueError(f"Last element not fixed {interim=}, {self=}, {diff_rows=}, {descs=}, {bubbles=}")
        #     pass

        interim = interim.rowrange(0, len(self) - 1).extend(1)
        assert interim.length_vector[:-1] == self.length_vector[:-1]
        assert len(interim.perm.trimcode) <= len(self) - 1, f"{interim.perm.trimcode=} {self.perm.trimcode=} {interim.perm=} {self.perm=} {interim=} {self=}"
        # assert in pull out var
        assert (self.perm, len(self)) in (ASx(interim.perm, len(self) - 1) * ASx(uncode([]), 1)).keys(), f"{self=}, {interim=} not zero correct"
        return interim.rowrange(0, len(self) - 1)

    @staticmethod
    def complete_sym_perms_op(orig_perm, p, k):
        from schubmult.utils.perm_utils import has_bruhat_descent

        orig_perm = Permutation(orig_perm)
        total_dict = {orig_perm: []}
        up_perm_list = [(orig_perm, 1000000000)]
        for pp in range(p):
            perm_list = []
            new_total_dict = {}
            for up_perm, last in up_perm_list:
                pos_list = [j for j in range(len(orig_perm)) if up_perm[j] == orig_perm[j]]
                for i in range(k):
                    for j in pos_list:
                        if has_bruhat_descent(up_perm, i, j):
                            new_perm_add = up_perm.swap(i, j)
                            perm_list += [(new_perm_add, up_perm[j])]
                            new_total_dict[new_perm_add] = [(i + 1, j + 1)] + total_dict[up_perm]
            up_perm_list = perm_list
            total_dict = new_total_dict
        return total_dict

    @staticmethod
    def complete_sym_perms(orig_perm, p, k):
        from schubmult.utils.perm_utils import has_bruhat_ascent

        orig_perm = Permutation(orig_perm)
        total_dict = {orig_perm: []}
        up_perm_list = [(orig_perm, 1000000000)]
        for pp in range(p):
            perm_list = []
            new_total_dict = {}
            for up_perm, last in up_perm_list:
                pos_list = [j for j in range(len(orig_perm) + 10) if up_perm[j] == orig_perm[j]]
                for i in range(k):
                    for j in pos_list:
                        if has_bruhat_ascent(up_perm, i, j):
                            new_perm_add = up_perm.swap(i, j)
                            perm_list += [(new_perm_add, up_perm[j])]
                            new_total_dict[new_perm_add] = total_dict[up_perm] + [(i + 1, j + 1)]
            up_perm_list = perm_list
            total_dict = new_total_dict
        return total_dict

    def vertical_cut(self, row):
        if row <= 0 or row > len(self):
            raise ValueError("Row out of range")
        front = RCGraph([*self[:row]]).extend(len(self) - row)
        for _ in range(len(self) - row):
            front = front.zero_out_last_row()
        if row == len(self):
            back = RCGraph()
        else:
            back = self.rowrange(row, len(self))
        return (front, back)

    # @staticmethod
    # def cap_perm_length(perm, spot):
    #     index = spot - 1
    #     spot_perm = perm
    #     top_spot = len(perm)
    #     while spot_perm[index] != len(spot_perm):
    #         spot_perm = spot_perm.swap(top_spot, top_spot + 1)
    #         top_spot -= 1
    #     return

    def right_zero_act(self, debug=False):
        # NOTE THAT THIS IS STILL USING THE OLD METHOD
        if self.perm.inv == 0:
            return {RCGraph([*self, ()])}
    
        if self in RCGraph._z_cache:
            return RCGraph._z_cache[self]
        # assert len(self.perm.trimcode) <= len(self), f"{self=}, {self.perm=}"
        up_perms = ASx(self.perm, len(self)) * ASx(uncode([0]), 1)

        rc_set = set()

        # if len(self.perm.trimcode) < len(self):
        #     rc_set0 = self.rowrange(0, len(self) - 1).right_zero_act(debug=debug)
        #     return {rc.extend(1) for rc in rc_set0}

        debug_print(f"{self=} {self.perm=}", debug=debug)

        # should be able to move in descents and reverse insert
        interim = RCGraph([*self])
        #while interim.perm[len(self)] != len(interim.perm) or len(interim.perm) != len(self.perm) + 1:
        # insert anything and remove descents
        # put up descents
        # remove n-1 roots
        # remove remaining descents

        if False:
            for (perm, _) in up_perms.keys():
                interim = RCGraph([*self, tuple(range(len(self.perm), len(self), -1))])
                # print(f"{interim.perm=}")
                # print(f"{lower_perm=}")
                print(f"{interim=}, {perm=}")
                ref_dict = RCGraph.complete_sym_perms(perm, interim.perm.inv - perm.inv, len(self) + 1)
                #print(f"{ref_dict=}")
                ref_path = ref_dict[interim.perm]
                #new_ref_path = [(root[0], root[1] - 1) for root in ref_path if root[0] <= len(self) + 1]
                # THIS SHOULD WORK BUT DOESN'T
                rc, diff_rows = interim.reverse_kogan_kumar_insert(len(self) + 1, ref_path, return_rows=True)
                print(f"Got {rc=}, {diff_rows=}, {ref_path=}")
                if len(rc[-1]) == 0:
                    rc_set.add(rc)
                else:
                    descs = []
                    assert len(rc) == len(self) + 1
                    diff_rows2 = []
                    while len(rc.perm.trimcode) > len(self):
                        descs += [len(rc.perm.trimcode)]
                        rc, row = rc.exchange_property(len(rc.perm.trimcode), return_row=True)
                        #diff_rows2 += [row]
                    #int_rc = RCGraph([*rc[:-1], tuple(sorted(descs, reverse=True))])
                    print(f"Took out exchange {rc=}")
                    # diff from self
                    diff_rows2 = []
                    for i in range(len(self)):
                        if len(rc[i]) != len(self[i]):
                            diff_rows2 += [i + 1]*(len(self[i]) - len(rc[i]))
                    rc2 = rc.kogan_kumar_insert(len(perm.trimcode), diff_rows2, debug=debug, backwards=True)
                    # remove the rows from self
                    # self_rc = RCGraph([*self, tuple(sorted(descs, reverse=True))])
                    # ref_dict2 = RCGraph.complete_sym_perms(int_rc.perm, self_rc.perm.inv - int_rc.perm.inv, len(self) + 1)
                    # print(f"{ref_dict2=}")
                    # print(f"{self_rc=}")
                    # print(f"{int_rc=}")
                    # print(f"{self_rc.perm=}")
                    # ref_path2 = ref_dict2[self_rc.perm]
                    
                    # rc2 = self_rc.reverse_kogan_kumar_insert(len(self) + 1, ref_path2, debug=debug).rowrange(0,len(self)).extend(1)
                    print(f"Got {rc2=}")
                    print(f"{rc2.zero_out_last_row()=} {self=}")
                    assert rc2.zero_out_last_row() == self
                    rc_set.add(rc2)
                


        # rc_set.add(self.extend(1))

        # # insert monks < len(self)
        # # pull out the r root
        # for row in range(1, len(self.perm.trimcode) + 1):
        #     test_rc, (ref,) = self.kogan_kumar_insert(len(self), [row], debug=debug, return_reflections=True)
        #     # find if it has a root that makes tn length one longer
        #     # root that goes from r to ref[1]
        #     debug_print("Inserted at row ", row)
        #     debug_print(test_rc)
        #     debug_print("Reflection is ", ref)
        #     found = False
        #     for row_below in range(1, row):
        #         if found:
        #             break
        #         for j in range(test_rc.cols + len(self) + 5):
        #             if not test_rc[row_below - 1, j]:
        #                 a, b = test_rc.right_root_at(row, j + 1)
        #                 if b == ref[1] and a == len(self.perm.trimcode) + 1:
        #                     debug_print(f"Toggling at {row, j + 1} with root {a, b}")
        #                     new_rc = test_rc.toggle_ref_at(row, j + 1)

        #                     debug_print(new_rc)
        #                     found = True
        #                     break
        #                     #else:
        #                     #
        #     if found:
        #         new_rc = new_rc.reverse_kogan_kumar_insert(len(self), ).extend(1)
        #         rc_set.add(new_rc)
        #         debug_print("Added", debug=debug)
        #     if not found:
        #         debug_print("Did not find root to toggle", debug=debug)
        #         debug_print(test_rc)
        #         debug_print(f"{self=}, {row=}, {ref=}")
        #         #raise ValueError("Did not find root to toggle")


        for perm, _ in up_perms.keys():
            # rc = RCGraph([*self, tuple(range(len(self.perm), len(self), -1))])

            for rc in RCGraph.all_rc_graphs(perm, len(self) + 1, weight=(*self.length_vector, 0)):
                if rc.zero_out_last_row() == self:
                    rc_set.add(rc)
                    debug_print("Added", debug=debug)
                    debug_print(rc, debug=debug)
                else:
                    debug_print("Skipped", debug=debug)
                    debug_print(rc, debug=debug)
                    debug_print("because", debug=debug)
                    debug_print(rc.zero_out_last_row(), debug=debug)
                    debug_print("Is not", debug=debug)
                    debug_print(self, debug=debug)

            # r = len(perm.trimcode)
            # s = max([i + 1 for i in range(r, len(perm) + 1) if perm[i] < perm[r - 1]])

            # #if r < len(self):

            # down_perm = perm.swap(r - 1, s - 1)
            # (r2, s2) = sorted([i + 1 for i in range(len(perm)) if down_perm[i] != self.perm[i]])
            # print(f"{r2, s2=} {r=} {perm=} {self.perm=} {down_perm=}")
            # withdrawn_rc, (row,) = self.reverse_kogan_kumar_insert(r - 1, [(r2, s2)], return_rows=True)
            # rc = withdrawn_rc.extend(1).kogan_kumar_insert(r, [row])
            # print(f"Got {rc=}, {rc.perm=}, {perm=}, {row=}, {withdrawn_rc=}")
            # #if rc.zero_out_last_row() == self:
            # assert rc.perm == perm, f"{rc=}, {rc.perm=}, {perm=}, {row=}, {withdrawn_rc=}"
            # rc_set.add(rc)
            # print("Added")
            # print(rc)

            # # top_perm = self.perm * Permutation.ref_product(*tuple(range(len(self.perm), len(self), -1)))
            # top_perm = rc.perm
            # ref_path = RCGraph.complete_sym_perms(perm, top_perm.inv - perm.inv, len(self)+1)[top_perm]
            # ref_path2 = [(ref[0], ref[1] - 1) for ref in ref_path if ref[0] <= len(self)]
            # # #print(f"{ref_path=}")
            # # new_ref_path = [(root[0], root[1] - 1) for root in ref_path if root[0] <= len(self)]
            # # #print(f"{new_ref_path=}")
            # lower_rc, diff_rows = self.reverse_kogan_kumar_insert(len(self), ref_path2, return_rows=True)

            # #diff_rows_low = [r for r in diff_rows if r <= len(self)]
            # rc = lower_rc.kogan_kumar_insert(len(self) + 1, diff_rows, debug=debug).rowrange(0,len(self)).extend(1)
            # print(rc)
            # assert rc.perm == perm
            # rc_set.add(rc)
            # lower_rc = RCGraph([*lower_rc.rowrange(0,len(self)),tuple(range(len(self.perm), len(self), -1))])
            # lower_rc = lower_rc.kogan_kumar_insert(len(self) + 1, diff_rows_low, debug=debug).rowrange(0,len(self)).extend(1)
            # print(f"After pullout {lower_rc=}, {diff_rows=}, {ref_path=}, {up_perms=}, {diff_rows_low=}")
            # assert (lower_rc.perm, len(lower_rc)) in up_perms.keys(), f"{lower_rc=}, {self=}, {diff_rows=}, {ref_path=}, {up_perms=}"
            # # rc_set.add(lower_rc)
            # assert lower_rc.length_vector[:-1] == self.length_vector, f"{lower_rc=}, {self=}, {diff_rows=}, {ref_path=}, {up_perms=}"
            # rc = lower_rc.kogan_kumar_insert(len(self) + 1, diff_rows, debug=debug).extend(1)
            # rc_set.add(rc)
            # nperm = perm
            # descs = []
            # while (len(nperm.trimcode)) > len(self):
            #     descs.append(len(nperm.trimcode))
            #     nperm = nperm.swap(len(nperm.trimcode) - 1, len(nperm.trimcode))

            # descs.sort()
            # rc = self.kogan_kumar_insert(len(self) + 1, [len(self) + 1]* (len(self.perm)-len(self)), debug=debug)

            # assert rc.perm[len(self)] == len(rc.perm)
            # print(f"Afterpullout {rc=}")
            # ref_dict = RCGraph.complete_sym_perms(perm, rc.perm.inv - perm.inv, len(self) + 1)
            # if rc.perm not in ref_dict:
            #     debug_print(f"{perm=}")
            #     debug_print(f"{ref_dict=}")
            #     debug_print("rc.perm=", rc.perm)
            #     continue
            # refl = ref_dict[rc.perm]
            # refl2 = []
            # for r in refl:
            #     assert r[1] > len(self)
            #     if r[1] > len(self):
            #         refl2.append((r[0],r[1]-1))
            # #rc = RCGraph([*self, ()])
            # rc, diff_rows = rc.reverse_kogan_kumar_insert(len(self) + 1, refl, return_rows=True)
            # print(f"{rc=}")
            # # print(f"{interim=}, {diff_rows=}, {refl=}, {up_perms=} {len(self)=}")
            # # # except AssertionError:
            # # #     continue
            # # print(f"{interim=}")
            # # # descs.sort(reverse=True)
            # # # desc_word = [*descs]

            # # #rc = interim.rowrange(0,len(self)).extend(1)
            # # diff_rows.sort(reverse=True)

            # # #for i, r in enumerate(diff_rows):
            # # rc = interim.rowrange(0,len(self)).extend(1)
            # rc = rc.rowrange(0, len(self)).extend(1)
            # rc = rc.kogan_kumar_insert(len(self), diff_rows, debug=debug)
            # print(f"{rc=}")
            # # # rc = rc.kogan_kumar_insert(len(self) + 1, diff_rows)
            # # # print(f"{self.perm=} {rc.perm=}")
            # # rc = rc.rowrange(0, len(self)).extend(1)
            # # if len(rc.perm.trimcode) > len(self) + 1:
            # #     debug_print(f"{rc=}")
            # #     debug_print(f"{rc.perm=}, {self.perm=}, {diff_rows=}, {refl=}, {up_perms=}")
            # #     continue
            # assert len(rc.perm.trimcode) <= len(self) + 1, f"{rc=}, {rc.perm=}, {self=}, {self.perm=}, {diff_rows=}, {refl=}, {up_perms=}, {self=}, {len(self)=}, {rc.perm.trimcode=}, {len(rc.perm.trimcode)=}"  # noqa: E501
            # assert rc.is_valid, f"{rc=}, {rc.perm=}, {self=}, {self.perm=}, {diff_rows=}, {refl=}, {up_perms=}"
            # # print("Finished with valid rc")
            # # print(rc)
            # #try:
            # assert rc.perm.inv == perm.inv
            #     #assert rc.perm == perm, f"{rc.perm=} {perm=}, {self.perm=}, {diff_rows=}"
            # assert len(rc) == len(self) + 1, f"{len(rc)=}, {len(self)=}, {rc=}, {self=}"
            # assert rc.zero_out_last_row() == self, f"{rc=} {rc.zero_out_last_row()=} {self=}, {descs=}"
            # assert len(rc) == len(self) + 1, f"{len(rc)=}, {len(self)=}, {rc=}, {self=}"
            # rc_set.add(rc)
            # assert rc.zero_out_last_row() == self, f"{rc=} {rc.perm=} {self=}, {self.perm=}, {diff_rows=}, {refl=}, {up_perms=} {rc.zero_out_last_row()=}"
            # except AssertionError:
            #     for rc in RCGraph.all_rc_graphs(perm, len(self) + 1, weight=tuple([*self.length_vector, 0])):
            #         if rc.length_vector[:-1] == self.length_vector and rc.zero_out_last_row() == self:
            #             assert rc in rc_set, f"{rc=} {rc.perm=} {self=}, {self.perm=}, {rc_set=}"
            #         # print(f"Added {rc=}, {rc.perm=}, {rc.length_vector=}")
        # CAN PUT THIS BACK IN TO CHECK
        # try:
        #     assert len(rc_set) == len(up_perms), f"{len(rc_set)=}, {len(up_perms)=}, {rc_set=}, {up_perms=}"
        #     # assert rc.zero_out_last_row() == self, f"{rc=} {rc.perm=} {self=}, {self.perm=}, {diff_rows=}, {refl=}, {up_perms=} {rc.zero_out_last_row()=}"
        # except AssertionError:
        #     assert len(rc_set) < len(up_perms), f"{rc_set=}, {up_perms=}, {self=}, {self.perm=}"
        #     for perm, lent in up_perms.keys():
        #         assert lent == len(self) + 1
        #         for rc0 in RCGraph.all_rc_graphs(perm, len(self) + 1, weight=(*self.length_vector, 0)):
        #             if rc0.length_vector[:-1] == self.length_vector:
        #                 assert rc0 in rc_set, f"{rc0=} {rc0.perm=} {self=}, {self.perm=}, {rc_set=}"
        #                 # print(f"Added {rc=}, {rc.perm=}, {rc.length_vector=}")
        RCGraph._z_cache[self] = rc_set
        return rc_set

    def __hash__(self):
        return self._hashcode

    @cache
    def bisect_left_coords_index(self, row, col, lo=0, hi=None):
        from bisect import bisect_left, bisect_right  # noqa: F401

        if hi is None:
            hi = len(self.perm_word)

        while lo < hi:
            mid = (lo + hi) // 2
            i, j = self.left_to_right_inversion_coords(mid)
            if i < row or (i == row and j > col):
                lo = mid + 1
            else:
                hi = mid
        return lo

    @cache
    def bisect_right_coords_index(self, row, col):
        raise NotImplementedError("Right bisect not implemented")

    def exchange_property(self, descent, left=False, return_row=False):
        debug_print(f"{descent=} in exchange property {self=} {self.perm=}")
        for i in range(len(self.perm_word)):
            if not left:
                a, b = self.left_to_right_inversion(i)
            else:
                a, b = self.left_to_right_left_inversion(i)
            if a == descent and b == descent + 1:
                row, col = self.left_to_right_inversion_coords(i)
                if return_row:
                    return self.toggle_ref_at(row, col), row
                return self.toggle_ref_at(row, col)
        raise ValueError("No such descent")

    @classmethod
    def generate_permuted_rcs(cls, perm, length, ordering):
        from schubmult.schub_lib.schub_lib import pull_out_var

        if length == 0:
            return {cls()}
        rc_set = set()

        for L, new_perm in pull_out_var(ordering[0], perm):
            extra = [perm[a] for a in range(ordering[0] - 1) if perm[a] != new_perm[a]]
            add_word = tuple(sorted([*L, *extra], reverse=True))
            rc_set_old = cls.generate_permuted_rcs(new_perm, length - 1, uncode(ordering.code[1:]))
            for rc in rc_set_old:
                new_rc = cls((tuple(add_word), *rc.shiftup(1)))
                rc_set.add(new_rc)
                assert new_rc.perm == perm, f"{new_rc=} {new_rc.perm=}, {new_perm=} {perm=}, {new_rc}"
        return rc_set

    @cache
    def _inversion_to_coord(self):
        roots_in = {}
        roots_out = {}
        for i in range(len(self)):
            for j in range(max(self[i]) - i - 1, -1, -1):
                a, b = self.right_root_at(i + 1, j + 1)
                if self.has_element(i + 1, j + 1):
                    roots_in[(a, b)] = (i + 1, j + 1)
                else:
                    roots_out[(a, b)] = (i + 1, j + 1)
        return roots_in, roots_out

    @cache
    def left_to_right_inversion(self, index):
        return self.right_root_at(*self.left_to_right_inversion_coords(index))

    @cache
    def left_to_right_left_inversion(self, index):
        # return self.left_root_at(*self.left_to_right_inversion_coords(index))
        raise NotImplementedError("Left inversions not implemented")

    @cache
    def left_to_right_inversion_coords(self, index, debug=False):
        if index < 0 or index >= len(self.perm_word):
            raise ValueError(f"Index {index} out of range {self.perm.inv}")
        index_find = 0
        for i in range(len(self)):
            index_find2 = index_find + len(self[i])
            if index_find2 > index:
                break
            index_find = index_find2
            # INDEX=0 WAH?

        debug_print(f"{index=} {index_find=} {i-1=} {self[i]=} {len(self.perm_word)} {(index - index_find)=} {self[i][(index - index_find)]}", debug=debug)
        debug_print(f"{(i + 1, self[i][(index - index_find)] - i)=}", debug=debug)
        return (i + 1, self[i][(index - index_find)] - i)

    @cache
    def max_of_row(self, r):
        if len(self[r - 1]) == 0:
            the_max = 1
            while self.right_root_at(r, the_max + 1) != (r + the_max, r + the_max + 1):
                the_max += 1
        else:
            the_max = max(self[r - 1]) - r + 1
        return the_max

    @classmethod
    def principal_rc(cls, perm, length):
        cd = perm.trimcode
        graph = []
        for i in range(len(cd)):
            row = tuple(range(i + cd[i], i, -1))
            graph.append(row)
        graph = [*graph, *[()] * (length - len(graph))]
        return cls(graph)

    # THE ZERO MAKES SCHUB PROD

    def prod_with_rc(self, other):
        if self.perm.inv == 0:
            return {RCGraph([*self, *other.shiftup(len(self))]): 1}
        num_zeros = max(len(other), len(other.perm))
        debug_print(f"Multiplying {self=} {other=} with {num_zeros=} zeros")
        assert len(self.perm.trimcode) <= len(self), f"{self=}, {self.perm=}"
        base_rc = self
        buildup_module = {base_rc: 1}

        for _ in range(num_zeros):
            new_buildup_module = {}
            for rc, coeff in buildup_module.items():
                new_buildup_module = add_perm_dict(new_buildup_module, dict.fromkeys(rc.right_zero_act(), coeff))
                # print("Ding")
            buildup_module = new_buildup_module
        # print("Barg")
        ret_module = {}

        for rc, coeff in buildup_module.items():
            new_rc = RCGraph([*rc[: len(self)], *other.shiftup(len(self))])
            assert len(new_rc) == len(self) + len(other)
            # if len(new_rc.perm.trimcode) > len(new_rc):
            #     new_rc = new_rc.extend(len(new_rc.perm.trimcode) - len(new_rc))
            if new_rc.is_valid and len(new_rc.perm.trimcode) <= len(new_rc):
                ret_module = add_perm_dict(ret_module, {new_rc: coeff})

        return ret_module

    def ring_act(self, elem):
        if isinstance(elem, FreeAlgebraElement):
            wd_dict = elem.change_basis(WordBasis)
            ret = {}
            for k, v in wd_dict.items():
                acted_element = {self: v}
                for a in reversed(k):
                    acted_element2 = {}
                    for k2, v2 in acted_element.items():
                        acted_element2 = add_perm_dict(acted_element2, dict.fromkeys(k2.act(a), v2))
                    acted_element = acted_element2
                ret = add_perm_dict(ret, acted_element)
            return ret
        raise ValueError(f"Cannot act by {type(elem)} {elem=}")

    def right_ring_act(self, elem):
        if isinstance(elem, FreeAlgebraElement):
            wd_dict = elem.change_basis(WordBasis)
            ret = {}
            for k, v in wd_dict.items():
                # print(f"{k=} {v=}")
                acted_element = {self: v}
                for a in reversed(k):
                    acted_element2 = {}
                    for k2, v2 in acted_element.items():
                        # print(f"{dict.fromkeys(k2.act(a), v2)=}")
                        acted_element2 = add_perm_dict(acted_element2, dict.fromkeys(k2.act(a), v2))
                    acted_element = acted_element2
                # print(f"{acted_element=}")
                ret = add_perm_dict(ret, acted_element)
            # print(f"{ret=}")
            return ret
        raise ValueError(f"Cannot act by {type(elem)} {elem=}")

    def right_act(self, p):
        st = self.transpose().act(p)
        return {rc.transpose() for rc in st}

    def act(self, p):
        pm = self.perm
        elem = FAS(pm, len(self))
        bumpup = FAS(uncode([p]), 1) * elem
        ret = set()
        for k, v in bumpup.items():
            perm2 = k[0]
            new_row = [pm[i] for i in range(max(len(pm), len(perm2))) if pm[i] == perm2[i + 1]]
            new_row.sort(reverse=True)
            nrc = RCGraph([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in self]])
            assert nrc.perm == perm2
            ret.add(nrc)
        assert ret == self.iterative_act(p), f"{ret=}\n{self.iterative_act(p)=}"
        return ret

    def iterative_act(self, p, insert=True):
        if p == 0:
            if insert:
                return {RCGraph([(), *[tuple([row[i] + 1 for i in range(len(row))]) for row in self]])}
            return {RCGraph(self)}
        last = self.iterative_act(p - 1, insert=insert)
        ret = set()
        for rc in last:
            # top row has som
            last_desc = 0
            old_perm = ~rc.perm
            if len(rc[0]) > 0:
                last_desc = rc[0][0]
            mx = len(rc)
            for row in rc:
                if len(row) > 0:
                    mx = max(mx, max(row) + 2)
            for i in range(last_desc + 1, mx + 1):
                if old_perm[i - 1] < old_perm[i]:
                    new_perm = old_perm.swap(i - 1, i)
                    if max((~new_perm).descents()) + 1 > len(rc):
                        continue
                    new_top_row = [i, *rc[0]]
                    new_rc = RCGraph([tuple(new_top_row), *rc[1:]])
                    ret.add(new_rc)
        return ret

    def __ge__(self, other):
        return not (self < other)

    def __gt__(self, other):
        return not (self <= other)

    @property
    def is_principal(self):
        return self.perm == uncode(self.length_vector)

    @property
    def inv(self):
        return self.perm.inv

    def _sympystr(self, printer=None):
        from sympy.printing.str import StrPrinter

        if not printer:
            printer = StrPrinter()
        return printer._print_MatrixBase(self)

    def _pretty(self, printer):
        return printer._print_MatrixBase(self)

    @property
    def rows(self):
        return len(self)

    @property
    def width(self):
        return self.cols

    @property
    def height(self):
        return self.rows

    @cached_property
    def cols(self):
        return max([self[i][0] - i if len(self[i]) > 0 else 0 for i in range(len(self))]) if len(self) > 0 else 0

    def __getitem__(self, key):
        # FLIPPED FOR PRINTING
        if isinstance(key, int):
            return self._getitem(key)
        if isinstance(key, tuple):
            i, j = key
            if not self.has_element(i + 1, self.cols - j):
                return " "
            return i + self.cols - j
        is_slice = isinstance(key, slice)

        if is_slice:
            return self._getitem(key)

        raise ValueError(f"Bad indexing {key=}")

    def _format_str(self, printer=None) -> str:
        from sympy.printing import StrPrinter

        if not printer:
            printer = StrPrinter()
        # Handle zero dimensions:
        if self.rows == 1:
            return "RCGraph([%s])" % self.table(printer, rowsep=",\n")  # noqa: UP031
        return "RCGraph([\n%s])" % self.table(printer, rowsep=",\n")  # noqa: UP031

    def table(self, printer, rowstart="|", rowend="|", rowsep="\n", colsep=" ", align="right"):
        table: list[list[str]] = []
        # Track per-column max lengths for pretty alignment
        maxlen = [0] * self.cols
        for i in range(self.rows):
            table.append([])
            for j in range(self.cols):
                s = printer._print(self[i, j]) if self[i, j] is not None else " "
                table[-1].append(s)
                maxlen[j] = max(len(s), maxlen[j])
        # Patch strings together
        align = {
            "left": "ljust",
            "right": "rjust",
            "center": "center",
            "<": "ljust",
            ">": "rjust",
            "^": "center",
        }[align]
        res = [""] * len(table)
        for i, row in enumerate(table):
            for j, elem in enumerate(row):
                row[j] = getattr(elem, align)(maxlen[j])
            res[i] = rowstart + colsep.join(row) + rowend
        return rowsep.join(res)

    def __lt__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        if self.perm.trimcode < other.perm.trimcode:
            return True
        if self.perm.trimcode > other.perm.trimcode:
            return False
        for i in range(self.perm.inv):
            a, b = self.perm.right_root_at(i)

            if self.inversion_label(a - 1, b - 1) < other.inversion_label(a - 1, b - 1):
                return True
            if self.inversion_label(a - 1, b - 1) > other.inversion_label(a - 1, b - 1):
                return False
        return False

    def __le__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return self < other or self == other
