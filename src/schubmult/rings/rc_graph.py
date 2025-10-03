import logging
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

init_logging(debug=True)
logger = get_logger(__name__)
# logging.basicConfig(level=logging.DEBUG,
#                      format='%(asctime)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s',)


def _is_row_root(row, root):
    return root[0] <= row and root[1] > row


FA = FreeAlgebra(WordBasis)
FAS = FA


class RCGraph(Printable, tuple):
    # def __str__(self):
    #     return super(Printable).__str__()

    # def _pretty(self, printer=None):
    #     return printer._print_DiagramGrid(self)

    def __eq__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return tuple(self) == tuple(other)

    @property
    def is_valid(self):
        if self.perm.inv != len(self.perm_word):
            return False
        # if max(self.perm.descents()) + 1 > len(self):
        #     return False
        return True

    def shiftup(self, shift):
        return [tuple([a + shift for a in rrow]) for rrow in self]

    @cache
    def right_root_at(self, i, j):
        from bisect import bisect_left

        start_root = (i + j - 1, i + j)
        if i > len(self):
            return start_root
        row = self[i - 1]
        revved = [*row]
        revved.reverse()

        index = bisect_left(revved, i + j - 1)
        assert index >= len(revved) or revved[index] != i + j - 1 or self.has_element(i, j)
        perm = Permutation.ref_product(*revved[:index])
        start_root = (perm[start_root[0] - 1], perm[start_root[1] - 1])
        lower_perm = Permutation([])

        for rrow in self[i:]:
            lower_perm *= Permutation.ref_product(*rrow)

        return ((~lower_perm)[start_root[0] - 1], (~lower_perm)[start_root[1] - 1])

    @cache
    def left_root_at(self, i, j):
        from bisect import bisect_left

        start_root = (i + j - 1, i + j)
        if i > len(self):
            return start_root
        row = self[i - 1]
        revved = [*row]
        # revved.reverse()

        index = bisect_left(revved, i + j - 1)
        assert index >= len(revved) or revved[index] != i + j - 1 or self.has_element(i, j)
        perm = Permutation.ref_product(*revved[:index])
        start_root = (perm[start_root[0] - 1], perm[start_root[1] - 1])
        lower_perm = Permutation([])
        if i > 1:
            for rrow in reversed(self[: i - 1]):
                lower_perm *= Permutation.ref_product(*rrow)

        return ((lower_perm)[start_root[0] - 1], (lower_perm)[start_root[1] - 1])

    def reverse_kogan_kumar_insert(self, descent, reflection_path, return_rows = False, debug=True):
        from schubmult.utils.perm_utils import has_bruhat_descent

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

        def is_relevant_crossing(root, prm):
            # min_root = max(pair_dict.keys())

            if root[0] not in pair_dict:
                if root[0] not in pair_dict_rev or pair_dict_rev.get(root[0], 0) != pair_dict_rev.get(root[1], 0):
                    return False
                return True
            return root[1] in pair_dict[root[0]] and has_bruhat_descent(prm, root[0] - 1, root[1] - 1)

        # may have to add q, s or a_i, q
        def is_relevant_noncrossing(root):
            # top, bottom = max(root[1], root[0]), min(root[0], root[1])
            return (root[0] <= descent and descent < root[1] and root[1] not in pair_dict_rev) or (root[0] in pair_dict_rev and root[1] > descent and root[1] not in pair_dict_rev)

        # Add this intersection. If we are in the first case, insert (s, q) into the sequence (ai, bi) in the rightmost position, such that ai’s remain nondecreasing in the # noqa: RUF003
        # sequence. ((s, q) are the rows where the two strands shown in Figure 3 originate.) If
        # we are in the second case, add (ai, q) just before where (a, bi) is in the sequence.

        working_rc = self
      #  if debug:
            # print("Starting with")
            # print(working_rc)
            # print(working_rc.perm)
        rows = []
        for row in range(1, len(self) + 1):
            for col in range(self.cols, 0, -1):
                if working_rc.has_element(row, col):
                    a, b = working_rc.right_root_at(row, col)
                    if is_relevant_crossing((a, b), working_rc.perm):
                        working_rc = working_rc.toggle_ref_at(row, col)
                        a2 = a
                        if a2 in pair_dict_rev:
                            a2 = pair_dict_rev[a2]
                        pair_dict[a2].remove(b)
                        if len(pair_dict[a2]) == 0:
                            del pair_dict[a2]
                        del pair_dict_rev[b]
                        rows.append(row)
                        for col2 in range(1, col):
                            if not working_rc.has_element(row, col2):
                                a2, b2 = working_rc.right_root_at(row, col2)
                                if a2 > b2:
                                    continue
                                if is_relevant_noncrossing((a2, b2)):
                                    if a2 <= descent:
                                        assert b2 not in pair_dict
                                        if a2 not in pair_dict:
                                            pair_dict[a2] = set()
                                        pair_dict[a2].add(b2)
                                        pair_dict_rev[b2] = a2
                                        working_rc = working_rc.toggle_ref_at(row, col2)
                                        rows.pop()
                                    else:
                                        assert a2 in pair_dict_rev
                                        assert b2 not in pair_dict_rev
                                        a = pair_dict_rev[a2]
                                        pair_dict[a].add(b2)
                                        pair_dict_rev[b2] = a
                                        working_rc = working_rc.toggle_ref_at(row, col2)
                                        rows.pop()
                                    break

        assert len(pair_dict_rev) == 0, f"{pair_dict=}, {pair_dict_rev=}, {working_rc=}"
    
        if return_rows:
            return working_rc, rows
        return working_rc

    # def insert_ref_rows(self, descent, rows, reflection_path):
    #     dict_by_a = {}
    #     dict_by_b = {}
    #     # row is descent
    #     # inserting times
    #     pair_dict = {}
    #     for ref in reflection_path:
    #         a, b = ref
    #         if a not in pair_dict:
    #             pair_dict[a] = set()
    #         pair_dict[a].add(b)
    #     pair_dict_rev = {}
    #     for a, b_list in pair_dict.items():
    #         for b in b_list:
    #             assert a <= descent and descent < b  # noqa: PT018
    #             pair_dict_rev[b] = a
    #     working_rc = RCGraph([*self])
    #     if len(rows) == 0:
    #         return self
    #     rows_grouping = {}
    #     for r in rows:
    #         rows_grouping[r] = rows_grouping.get(r, 0) + 1
    #     if max(rows) > len(working_rc):
    #         working_rc = working_rc.extend(max(rows) - len(working_rc))
    #     rows = sorted(rows, reverse=True)
    #   #  if debug:
    #         # print(f"inserting {rows=}")
    #     for row in sorted(rows_grouping.keys(), reverse=True):
    #         num_times = rows_grouping[row]
    #       #  if debug:
    #             # print(f"Inserting {row=} {num_times=}")
    #             # print(working_rc)
    #             # print(f"{working_rc.perm.inv=}, {self.perm.inv=}")
    #         last_working_rc = working_rc
    #         working_rc = working_rc._kogan_kumar_insert_row(row, descent, dict_by_a, dict_by_b, num_times)

    #         if not working_rc.is_valid:
    #             working_rc = working_rc._kogan_kumar_rectify(row, descent, dict_by_a, dict_by_b)  # minus one?
    #       #  if debug:
    #             # print("Next iteration")
    #             # print(working_rc)
    #             # print(f"{working_rc.perm.inv=}, {self.perm.inv + index + 1=}")

    #         try:
    #             assert len(working_rc[row - 1]) == len(last_working_rc[row - 1]) + num_times
    #         except AssertionError:
    #             # print("Assertion failed")
    #             # print(working_rc)
    #             # print(f"{working_rc.perm.inv=}, {self.perm.inv + total_num=}")
    #             # print(f"{dict_by_a=}, {dict_by_b=}")
    #             # print(f"{working_rc.perm=}, {self.perm=}")
    #           #  if debug:
    #           #      raise
    #           #  self.kogan_kumar_insert(descent, rows, debug=True)
    #             raise
    #     if return_reflections:
    #         reflections = []
    #         for a in sorted(dict_by_a.keys()):
    #             reflections = [*reflections, *[(a, b) for b in dict_by_a[a]]]
    #         # print(f"Returning {working_rc=}, {reflections=}")
    #         return working_rc, tuple(reflections)
    #     return working_rc  # , tuple(reflections)

    # def find_ref_spots(self, descent, reflection_path, rows, debug=True):
    #     # from schubmult.utils.perm_utils import has_bruhat_descent

    #     # pair_sequence = sorted(pair_sequence, key=lambda x: x[0]
    #     pair_dict = {}
    #     for ref in reflection_path:
    #         a, b = ref
    #         if a not in pair_dict:
    #             pair_dict[a] = set()
    #         pair_dict[a].add(b)
    #     pair_dict_rev = {}
    #     # ref_by_index = {}
    #     for a, b_list in pair_dict.items():
    #         for b in b_list:
    #             assert a <= descent and descent < b  # noqa: PT018
    #             pair_dict_rev[b] = a

    #     def is_relevant_noncrossing(root):
    #         # min_root = max(pair_dict.keys())
    #         # return root[1] in pair_dict_rev and pair_dict_rev[root[1]] == root[0]
    #         if root[0] not in pair_dict:
    #             if root[0] not in pair_dict_rev or pair_dict_rev.get(root[0], 0) != pair_dict_rev.get(root[1], 0):
    #                 return False
    #             return True
    #         return root[1] in pair_dict[root[0]]
    #         # if root[0] not in pair_dict:
    #         #     if root[1] in pair_dict or pair_dict_rev.get(root[0], 0) != pair_dict_rev.get(root[1], 0):
    #         #         return False
    #         #     return True
    #         # return root[1] in pair_dict[root[0]] and has_bruhat_descent(prm, root[0] - 1, root[1] - 1)

    #     # may have to add q, s or a_i, q
    #     def is_relevant_crossing(root):
    #         # top, bottom = max(root[1], root[0]), min(root[0], root[1])
    #         return ((root[0] > descent and descent >= root[1]) and (root[0] not in pair_dict_rev or (root[1] in pair_dict_rev and root[0] > descent and root[0] not in pair_dict_rev)))

    #     # Add this intersection. If we are in the first case, insert (s, q) into the sequence (ai, bi) in the rightmost position, such that ai’s remain nondecreasing in the # noqa: RUF003
    #     # sequence. ((s, q) are the rows where the two strands shown in Figure 3 originate.) If
    #     # we are in the second case, add (ai, q) just before where (a, bi) is in the sequence.

    #     working_rc = self
    #   #  if debug:
    #         # print("Starting with")
    #         # print(working_rc)
    #         # print(working_rc.perm)
    #     times = 0
    #     while len(rows) > 0 and times < 10:
    #         for row in sorted(rows):
    #             print(row)
    #             for col in range(self.cols + descent + 1, -1, -1):
    #                 if not working_rc.has_element(row, col):
    #                     a, b = working_rc.right_root_at(row, col)
    #                     print(f"{a,b}")
    #                     if is_relevant_noncrossing((a, b)):
    #                         working_rc = working_rc.toggle_ref_at(row, col)
    #                         print(working_rc)
    #                         a2 = a
    #                         if a2 in pair_dict_rev:
    #                             a2 = pair_dict_rev[a2]
    #                         pair_dict[a2].remove(b)
    #                         if len(pair_dict[a2]) == 0:
    #                             del pair_dict[a2]
    #                         del pair_dict_rev[b]
    #                         rows.remove(row)
    #                         for col2 in range(1, col):
    #                             if working_rc.has_element(row, col2):
    #                                 a2, b2 = working_rc.right_root_at(row, col2)
    #                                 if a2 < b2:
    #                                     continue
    #                                 if is_relevant_crossing((a2, b2)):
    #                                     if b2 <= descent:
    #                                         assert b2 not in pair_dict
    #                                         if b2 not in pair_dict:
    #                                             pair_dict[a2] = set()
    #                                         pair_dict[a2].add(b2)
    #                                         pair_dict_rev[b2] = a2
    #                                         working_rc = working_rc.toggle_ref_at(row, col2)
    #                                         rows.append(row)
    #                                     else:
    #                                         assert b2 in pair_dict_rev
    #                                         assert a2 not in pair_dict_rev
    #                                         a = pair_dict_rev[b2]
    #                                         pair_dict[a].add(a2)
    #                                         pair_dict_rev[a2] = a
    #                                         working_rc = working_rc.toggle_ref_at(row, col2)
    #                                         rows.append(row)
                                        
    #                         break
                    
    #             print(rows)
    #         times += 1

    #     assert len(pair_dict_rev) == 0, f"{pair_dict=}, {pair_dict_rev=}, {working_rc=}"
    
        
    #     return working_rc


    @cache
    def inversion_label(self, i, j):
        if i >= j:
            raise ValueError("i must be less than j")
        if self.perm[i] < self.perm[j]:
            raise ValueError("Not an inversion")
        for index in len(self.perm_word):
            if self.left_to_right_inversion(index) == (i + 1, j + 1):
                return self.left_to_right_inversion_coord(index)[0]
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

    # def __len__(self):
    #     return len(self.P)

    def __new__(cls, *args):
        new_args = tuple(tuple(arg) for arg in args)
        return RCGraph.__xnew_cached__(cls, *new_args)

    @staticmethod
    @cache
    def __xnew_cached__(_class, *args):
        return RCGraph.__xnew__(_class, *args)

    @staticmethod
    def __xnew__(_class, *args):
        return tuple.__new__(_class, *args)

    def __init__(self, *args):
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

    def weight_word(self):
        perm = self.perm
        nz = len([a for a in perm.trimcode if a != 0])
        root_dict = {perm.right_root_at(index): index for index in range(perm.inv)}
        result_word = [9] * perm.inv
        index = 0
        perm_word = self.perm_word
        for i, row in enumerate(self):
            for _ in range(len(row)):
                root = perm.right_root_at(index, word=perm_word)
                result_word[root_dict[root]] = nz - perm.code_index_of_index(index)
                index += 1
        result_word.reverse()
        return tuple(result_word)

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

    # def __new__(cls, *args, **kwargs):
    #     return tuple.__new__(cls, *args)

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

    _cache_by_weight = {}

    @classmethod
    @cache
    def all_rc_graphs(cls, perm, length=-1, weight=None):
        if length > 0 and length < len(perm.trimcode):
            raise ValueError("Length must be at least the last descent of the permutation")
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

    # this doesn't work
    # strong exchange property
    def extract_row(self, row):
        debug = False
        if row <= 0 or row > len(self):
            raise ValueError("Row out of range")
        interim = RCGraph([*self[:row]])

        diff_rows = []

      #  if debug:
            # print("Zeroing out last row")
            # print(self)
            # print("-------")

        diff_row = None
        descs = []

        # while len(interim.perm.trimcode) > row - 1:
        #     prev_interim = interim
        #     descs += [max(interim.perm.descents()) + 1]
        #     interim = interim.exchange_property(max(interim.perm.descents()) + 1)
        #     for i in range(len(interim)):
        #         if len(interim[i]) < len(prev_interim[i]):
        #             rw = i + 1
        #             diff_rows.append(rw)
        #             break
        if len(interim.perm.trimcode) == len(self):
            interim, diff_row = interim.exchange_property(len(interim.perm.trimcode), return_row=True)
        
        if diff_row:
            interim = interim.kogan_kumar_insert(row, diff_rows, debug=debug)

      #  if debug:
            # print("Got")
            # print(interim)
        # assert interim.length_vector[:-1] == self.length_vector[:-1]
        if row != 1:
            return RCGraph([*interim[: row - 1], *RCGraph(self[row:]).shiftup(-1)])
        return RCGraph([*RCGraph(self[row:]).shiftup(-1)])

    # def __le__(self, other):
    #     if not isinstance(other, RCGraph):
    #         return NotImplemented
    #     if not isinstance(other, RCGraph):
    #         return NotImplemented
    #     if self.length_vector < other.length_vector:
    #         return True
    #     if self.length_vector == other.length_vector and self.perm.bruhat_leq(other.perm) and self.perm != other.perm:
    #         return True
    #     if self.length_vector == other.length_vector and self.perm == other.perm and self.perm_word < other.perm_word:
    #         return True
    #     return self.length_vector == other.length_vector and self.perm == other.perm and self.perm_word == other.perm_word

    # def __lt__(self, other):
    #     if not isinstance(other, RCGraph):
    #         return NotImplemented
    #     if self.length_vector < other.length_vector:
    #         return True
    #     if self.length_vector == other.length_vector and self.perm.bruhat_leq(other.perm) and self.perm != other.perm:
    #         return True
    #     if self.length_vector == other.length_vector and self.perm == other.perm and self.perm_word < other.perm_word:
    #         return True
    #     return False

    def _kogan_kumar_insert_row(self, row, descent, dict_by_a, dict_by_b, num_times, pair_dict_rev = None, pair_dict=None, debug=True, start_index=0):
        working_rc = self
        if row > descent:
            raise ValueError("All rows must be less than or equal to descent")

        i = start_index
        num_done = 0

        while num_done < num_times:
            i += 1
            flag = False
          #  if debug:
                # print(f"Trying column {i=} {descent=} {row=} {num_done=} {num_times=}")
            if not working_rc.has_element(row, i):
                a, b = working_rc.right_root_at(row, i)
                if a < b:
                  #  if debug:
                        # print(f"root is {a, b}")
                    flag = False
                  #  if debug:
                        # print(f"{dict_by_b=}")
                    if _is_row_root(descent, (a, b)) and b not in dict_by_b and (pair_dict_rev is None or b in pair_dict_rev and a == pair_dict_rev[b]):
                        new_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[a] = dict_by_a.get(a, set())
                        dict_by_a[a].add(b)
                        dict_by_b[b] = a
                        if pair_dict_rev:
                            pair_dict[a].remove(b)
                            if len(pair_dict[a]) == 0:
                                del pair_dict[a]
                            del pair_dict_rev[b]
                        flag = True
                        working_rc = new_rc
                      #  if debug:
                            # print("Toggled a")
                            # print(working_rc)
                    elif a in dict_by_b and b not in dict_by_b:
                        new_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[dict_by_b[a]].add(b)
                        dict_by_b[b] = dict_by_b[a]
                        flag = True
                        working_rc = new_rc
                      #  if debug:
                            # print("Toggled b")
                            # print(working_rc)
                    elif pair_dict_rev and a in pair_dict_rev and b not in dict_by_b and b not in pair_dict_rev:
                        new_rc = working_rc.toggle_ref_at(row, i)
                        if pair_dict_rev[a] not in dict_by_a:
                            dict_by_a[pair_dict_rev[a]] = set()
                        dict_by_a[pair_dict_rev[a]].add(b)
                        dict_by_b[b] = pair_dict_rev[a]
                        pair_dict[pair_dict_rev[a]].remove(a)
                        if len(pair_dict[pair_dict_rev[a]]) == 0:
                            del pair_dict[pair_dict_rev[a]]
                        del pair_dict_rev[a]
                        flag = True
                        working_rc = new_rc
                      #  if debug:
                            # print("Toggled b")
                            # print(working_rc)
                    elif b in dict_by_b and a not in dict_by_b and a > descent:
                        new_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[dict_by_b[b]].add(a)
                        dict_by_b[a] = dict_by_b[b]
                        flag = True
                        working_rc = new_rc
                    elif pair_dict_rev and a>descent and a not in pair_dict_rev and a not in dict_by_b and b in pair_dict_rev:
                        new_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[pair_dict_rev[b]].add(a)
                        dict_by_b[a] = pair_dict_rev[b]
                        pair_dict[pair_dict_rev[b]].remove(b)
                        if len(pair_dict[pair_dict_rev[b]]) == 0:
                            del pair_dict[pair_dict_rev[b]]
                        del pair_dict_rev[b]
                        flag = True
                        working_rc = new_rc
                      #  if debug:
                            # print("Toggled c")
                            # print(working_rc)
                if flag:
                    num_done += 1
                    # assert last_rc.perm.inv + 1 == working_rc.perm.inv
                  #  if debug:
                        # print("Inserted")
                        # print(working_rc)
                # if not working_rc.is_valid:
                #     working_rc = working_rc._kogan_kumar_rectify(row - 1, descent, dict_by_a, dict_by_b)
        return working_rc

    def _kogan_kumar_rectify(self, row_below, descent, dict_by_a, dict_by_b):
        # print("In rectify")
        working_rc = self
        debug = True
        if row_below == 0:
            assert working_rc.is_valid
            return working_rc
        if working_rc.is_valid:
            return working_rc
        for j in range(working_rc.max_of_row(row_below) + 1, 0, -1):
            flag = False
            if working_rc.is_valid:
                return working_rc
            if working_rc.has_element(row_below, j):
                # print("Has element")
                a, b = working_rc.right_root_at(row_below, j)
                # print("root=", (a, b))
                if a > b:
                    # print("Entered")
                  #  if debug:
                        # print(f"Considering bad at {row_below, j}")
                        # print(f"{dict_by_a=}, {dict_by_b=}")
                        # print(f"root = ({a, b})")
                    if b in dict_by_a and a in dict_by_a[b]:
                        new_rc = working_rc.toggle_ref_at(row_below, j)
                        dict_by_a[b].remove(a)
                        if len(dict_by_a[b]) == 0:
                            del dict_by_a[b]
                        del dict_by_b[a]
                        working_rc = new_rc
                        flag = True
                        # print("Toggle bad a")
                        # print(working_rc)
                    elif a in dict_by_b and b in dict_by_b and dict_by_b[a] == dict_by_b[b]:
                        new_rc = working_rc.toggle_ref_at(row_below, j)
                        if new_rc.perm[dict_by_b[a] - 1] < new_rc.perm[a - 1]:
                            dict_by_a[dict_by_b[a]].remove(a)
                            del dict_by_b[a]
                            if len(dict_by_a[dict_by_b[b]]) == 0:
                                del dict_by_a[dict_by_b[b]]
                            # print("Toggle bad b")
                            flag = True
                        dict_by_a[dict_by_b[b]].remove(b)
                        del dict_by_b[b]
                        if len(dict_by_a[dict_by_b[a]]) == 0:
                            del dict_by_a[dict_by_b[a]]
                        # print("Toggle bad c")
                if flag:
                    working_rc = working_rc._kogan_kumar_insert_row(row_below, descent, dict_by_a, dict_by_b, num_times=1, debug=debug)
        return working_rc._kogan_kumar_rectify(row_below - 1, descent, dict_by_a, dict_by_b)

    def associative_kogan_kumar_insert(self, descent, rows, debug=True):
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

    def associative_product(self, other, debug=True):
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
        # interim = interim.kogan_kumar_insert(descent, rows, debug=debug)
        for i in range(len(diff_rows_stack) - 1, -1, -1):
            diff_rows = diff_rows_stack[i]
            prev_descent = desc_stack[i]
            other = other.associative_kogan_kumar_insert(prev_descent, diff_rows, debug=debug)
        return other

    # VERIFY
    def kogan_kumar_insert(self, descent, rows, debug=True, reflections = None, return_reflections=False):
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
        pair_dict = None
        pair_dict_rev = None
        if reflections:
            pair_dict = {}
            for ref in reflections:
                a, b = ref
                if a not in pair_dict:
                    pair_dict[a] = set()
                pair_dict[a].add(b)
            pair_dict_rev = {}
            for a, b_list in pair_dict.items():
                for b in b_list:
                    assert a <= descent and descent < b  # noqa: PT018
                    pair_dict_rev[b] = a

        for r in rows:
            rows_grouping[r] = rows_grouping.get(r, 0) + 1
        if max(rows) > len(working_rc):
            working_rc = working_rc.extend(max(rows) - len(working_rc))
        rows = sorted(rows, reverse=True)
      #  if debug:
            # print(f"inserting {rows=}")
        for row in sorted(rows_grouping.keys(), reverse=True):
            num_times = rows_grouping[row]
          #  if debug:
                # print(f"Inserting {row=} {num_times=}")
                # print(working_rc)
                # print(f"{working_rc.perm.inv=}, {self.perm.inv=}")
            last_working_rc = working_rc
            working_rc = working_rc._kogan_kumar_insert_row(row, descent, dict_by_a, dict_by_b, num_times, pair_dict=pair_dict, pair_dict_rev=pair_dict_rev, debug=debug)

            if not working_rc.is_valid:
                working_rc = working_rc._kogan_kumar_rectify(row, descent, dict_by_a, dict_by_b)  # minus one?
          #  if debug:
                # print("Next iteration")
                # print(working_rc)
                # print(f"{working_rc.perm.inv=}, {self.perm.inv + index + 1=}")

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
              #  self.kogan_kumar_insert(descent, rows, debug=True)
                raise
        if return_reflections:
            reflections = []
            for a in sorted(dict_by_a.keys()):
                reflections = [*reflections, *[(a, b) for b in dict_by_a[a]]]
            # print(f"Returning {working_rc=}, {reflections=}")
            return working_rc, tuple(reflections)
        return working_rc  # , tuple(reflections)

    @property
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

    def toggle_ref_at(self, i, j):
        from bisect import bisect_left

        assert i > 0 and j > 0  # noqa: PT018
        a, b = self.right_root_at(i, j)
        row = self[i - 1]
        rev_row = [*row]
        rev_row.reverse()
        index = bisect_left(rev_row, i + j - 1)
        if index >= len(rev_row):
            new_row = [i + j - 1, *row]
        else:
            if rev_row[index] == i + j - 1:
                new_row = [*row[: len(row) - index - 1], *row[len(row) - index :]]
            else:
                new_row = [*row[: len(row) - index], i + j - 1, *row[len(row) - index :]]
        return RCGraph([*self[: i - 1], tuple(new_row), *self[i:]])

    # # THIS IS KEY
    # # EXCHANGE PROPERTY GOES TO UNIQUE PERMUTATION
    # # KOGAN INSERT ENSURES WE GO UP PROPERLY

    def zero_out_last_row(self, debug=True):
        from schubmult.schub_lib.schub_lib import pull_out_var
        # this is important!
        # transition formula
        if len(self[-1]) != 0:
            raise ValueError("Last row not empty")
        if self.perm.inv == 0 or len(self.perm.trimcode) < len(self) - 1:
            return self.rowrange(0, len(self) - 1)

        # Let r
        # be the largest value such that w(r)>w(r+1)
        # and s
        # be the maximum value so that w(s)<w(r)
        # . We define

        interim = RCGraph([*self])
        diff_rows = []
        print("self")
        print(self)
        while(len(interim.perm.trimcode)) > len(self) - 1:
            interim, diff_row = interim.exchange_property(len(interim.perm.trimcode), return_row=True)
            print("Removed descent ", len(interim.perm.trimcode) + 1, " at row", diff_row)
            print(interim)
            diff_rows.append(diff_row)
        print(f"{diff_rows=}")
        # r = len(interim.perm.trimcode)
        # if r < len(self):
        #     return interim.rowrange(0, len(self) - 1)
        # s = max([i+1 for i in range(r, len(self.perm)) if self.perm[i] < self.perm[r-1]])

        # diff_row = None
        # for i in range(self.inv):
        #     # print(f"{i=} {self.left_to_right_inversion(i)=}")
        #     if self.left_to_right_inversion(i) == (r, s):
        #         interim = self.toggle_ref_at(*self.left_to_right_inversion_coord(i))
        #         diff_row = self.left_to_right_inversion_coord(i)[0]
        #         break
        # assert diff_row is not None, f"Could not find inversion {(r,s)} in {self.perm}, {self.perm_word=}, {self=}"
        # refl = RCGraph.complete_sym_perms_op(self.perm, len(diff_rows), len(self))[interim.perm]
        # interim, diff_rows = interim.reverse_kogan_kumar_insert(len(self), refl, return_rows=True)

        #perm_down = interim.perm
        #reflections = RCGraph.complete_sym_perms_op(self.perm, self.perm.inv - perm_down.inv, len(self))[perm_down.inv]
        interim, reflections = interim.kogan_kumar_insert(len(self) - 1, diff_rows, return_reflections=True)

        #interim = interim.kogan_kumar_insert(len(self) - 1, diff_rows, debug=debug)
        assert len(interim.perm.trimcode) < len(self)
        print("Final")
        print(interim)
        return interim.rowrange(0, len(self) - 1)

    # put the descents back
    # kogan reverse the non-descents
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

    def right_zero_act(self, debug=True):
        if self.perm.inv == 0:
            return {RCGraph([*self, ()])}

        up_perms = ASx(self.perm, len(self)) * ASx(uncode([0]), 1)

        rc_set = set()

        # print(f"{self=} {self.perm=}")
        rc_set.add(self.extend(1))
        for (perm, _) in up_perms.keys():
            if perm == self.perm:
                continue
            nperm = perm
            descs = []
            while(len(nperm.trimcode)) > len(self):
                descs.append(len(nperm.trimcode))
                nperm = nperm.swap(len(nperm.trimcode) - 1,len(nperm.trimcode))
            refl = RCGraph.complete_sym_perms_op(self.perm, self.perm.inv - nperm.inv, len(self))[nperm]
            interim, diff_rows = self.reverse_kogan_kumar_insert(len(self), refl, return_rows=True)
            descs.reverse()
            desc_word = [*descs]
            word_perm = Permutation.ref_product(*desc_word)
            refs = [word_perm.right_root_at(i) for i in range(word_perm.inv)]
            
            print(f"{descs=}")
            print(f"{refs=}")
            print(f"{diff_rows=}")
            print(interim)
            interim = interim.kogan_kumar_insert(len(self) + 1, diff_rows, reflections=refs)
            #??
            # print(f"{refs=}")
            # print(F"{diff_rows=}")
            # for i in range(len(diff_rows) - 1, -1, -1):
            #     row = diff_rows[i]
            #     for j in range(self.cols - 1, -1, -1):
            #         print(f"Trying {row, j+1}")
            #         if not interim[row - 1, j]:
            #             print(f"{interim=}")
            #             a, b = interim.right_root_at(row, j + 1)
            #             if (a, b) in refs:
            #                 interim = interim.toggle_ref_at(row, j + 1)
            #                 refs.remove((a, b))
            #                 print(f"Inserted {(a,b)}")
            #                 print(interim)
            #                 interim = interim._kogan_kumar_rectify(row - 1, len(self) + 1, {a: [b]}, {b: a})
            #                 print(f"Rectified")
            #                 print(interim)
            #                 break
            # the_graph = interim.kogan_kumar_insert(len(self) + 1, diff_rows, debug=debug).extend(1)
            # now we need to get the refs back in
            #rc_set.add(the_graph)
            # num = 1
            # refs = [(len(self) + 1, a + 1) for a in reversed(descs)]
            # while len(refs) > 0:
            #     for r in sorted(diff_rows):
            #         for j in range(self.cols):
            #             if not interim[r - 1, j] and interim.right_root_at(r, j + 1) in refs:
            #                 interim2, (ref,) = interim.kogan_kumar_insert(len(self) + 1, [r], return_reflections=True)
            #                 print(f"{ref=}")
            #                 if ref in refs:
            #                     interim = interim2
            #                     ind = refs.index(ref)
            #                     a, b = ref
            #                     swp = Permutation([]).swap(a - 1, b - 1)
            #                     for i, ref2 in enumerate(refs):
            #                         if i > ind:
            #                             refs[i] = (swp[ref2[0]-1], swp[ref2[1]-1])
            #                     refs.remove(ref)
            #                     diff_rows.remove(r)
            #                     break
            #                 continue                            
                            # print(f"{interim=}")
                            # #interim = interim._kogan_kumar_rectify(r - 1, descs[0] if len(descs) > 0 else 0, {}, {})
                            # num += 1
                            # break
                # interim = interim.kogan_kumar_insert(len(self) + 1, [r], debug=debug)
            
            rc = interim.extend(1)
            assert rc.is_valid
            rc_set.add(rc)
            #assert rc.zero_out_last_row() == self, f"{rc=} {rc.perm=} {self=}, {self.perm=}, {diff_rows=}, {refl=}, {up_perms=} {rc.zero_out_last_row()=}"
        for (perm, _) in up_perms.keys():
            for rc in RCGraph.all_rc_graphs(perm, len(self) + 1, weight=tuple([*self.length_vector, 0])):
                if rc.length_vector[:-1]  == self.length_vector and rc.zero_out_last_row() == self:
                    assert rc in rc_set, f"{rc=} {rc.perm=} {self=}, {self.perm=}, {rc_set=}"
                        # print(f"Added {rc=}, {rc.perm=}, {rc.length_vector=}")

        return rc_set

    @cache
    def bisect_left_coords_index(self, row, col):
        from bisect import bisect_left

        inversions = [self.left_to_right_inversion_coord(i) for i in range(len(self.perm_word))]
        inversions.sort(key=lambda x: (x[0], -x[1]))
        # print(f"{inversions=}")

        return bisect_left(inversions, (row, col), key=lambda x: (x[0], -x[1]))

    def exchange_property(self, descent, left=False, return_row = False):
        for i in range(len(self.perm_word) + 1):
            if not left:
                a, b = self.left_to_right_inversion(i)
            else:
                a, b = self.left_to_right_left_inversion(i)
            if a == descent and b == descent + 1:
                row, col = self.left_to_right_inversion_coord(i)
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
        return self.right_root_at(*self.left_to_right_inversion_coord(index))

    @cache
    def left_to_right_left_inversion(self, index):
        return self.left_root_at(*self.left_to_right_inversion_coord(index))

    @cache
    def left_to_right_inversion_coord(self, index):
        index2 = 0
        for i in range(self.rows):
            for j in range(self.cols - 1, -1, -1):
                if self[i, j]:
                    if index2 == index:
                        return (i + 1, j + 1)
                    index2 += 1
        assert False, "Index out of range"

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
        num_zeros = len(other) if other.perm.inv == 0 else max(len(other), len(other.perm))
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
        return printer._print_DiagramGrid(self)

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
        return max(max((row[i] + i for i in range(len(row))) if len(row) > 0 else [0]) for row in self) if len(self) > 0 else 0

    def __getitem__(self, key):
        if isinstance(key, int):
            return tuple(self)[key]
        if isinstance(key, tuple):
            i, j = key
            if not self.has_element(i + 1, j + 1):
                return None
            return i + j + 1
        is_slice = isinstance(key, slice)

        if is_slice:
            return tuple(tuple(self)[n] for n in range(len(self))[key])

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

    def __hash__(self):
        return hash(tuple(self))

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
