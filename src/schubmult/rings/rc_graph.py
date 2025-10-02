import logging
from functools import cache, cached_property

from sympy.printing.defaults import Printable

import schubmult.schub_lib.schub_lib as schub_lib
from schubmult.perm_lib import Permutation, uncode
from schubmult.rings import ASx
from schubmult.symbolic import S, prod
from schubmult.utils.perm_utils import add_perm_dict

from .free_algebra import FreeAlgebra, FreeAlgebraElement
from .free_algebra_basis import WordBasis
from .nil_hecke import NilHeckeRing

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


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

    def reverse_kogan_kumar_insert(self, descent, reflection_path, debug=False):
        # logger.debug("Here")
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
                if root[0] not in pair_dict or pair_dict_rev.get(root[0], 0) != pair_dict_rev.get(root[1], 0):
                    return False
                return True
            return root[1] in pair_dict[root[0]] and has_bruhat_descent(prm, root[0] - 1, root[1] - 1)

        # may have to add q, s or a_i, q
        def is_relevant_noncrossing(root):
            top, bottom = root[1], root[0]
            return (root[0] <= descent and descent < root[1] and root[1] not in pair_dict_rev) or (
                bottom in pair_dict_rev and top not in pair_dict_rev
            )

        # Add this intersection. If we are in the first case, insert (s, q) into the sequence (ai, bi) in the rightmost position, such that ai’s remain nondecreasing in the # noqa: RUF003
        # sequence. ((s, q) are the rows where the two strands shown in Figure 3 originate.) If
        # we are in the second case, add (ai, q) just before where (a, bi) is in the sequence.

        #new_rc = self
        #index = new_rc.max_of_row(i) - 2 if len(new_rc[i]) > 0 else -1
        # logger.debug(pair_dict)

        #while len(pair_dict_rev) > 0
        working_rc = self
        if debug:
            print("Starting with")
            print(working_rc)
            print(working_rc.perm)
        for row in range(1, len(self) + 1):
            for col in range(1, self.max_of_row(row) + 1):
                if working_rc.has_element(row, col):
                    a, b = working_rc.right_root_at(row, col)
                    if is_relevant_crossing((a,b), working_rc.perm):
                        working_rc = working_rc.toggle_ref_at(row, col)
                        pair_dict[a].remove(b)
                        if len(pair_dict[a]) == 0:
                            del pair_dict[a]
                        del pair_dict_rev[b]
                        for col2 in range(col + 1, self.max_of_row(row) + 1):
                            if working_rc.has_element(row, col2):
                                a2, b2 = working_rc.right_root_at(row, col2)
                                if a2 > b2:
                                    continue
                                if is_relevant_noncrossing((a2,b2)):
                                    if a2 <= descent:
                                        assert b2 not in pair_dict
                                        if a2 not in pair_dict:
                                            pair_dict[a2] = set()
                                        pair_dict[a2].add(b2)
                                        pair_dict_rev[b2] = a2
                                        working_rc = working_rc.toggle_ref_at(row, col2)
                                    else:
                                        assert a2 in pair_dict_rev
                                        assert b2 not in pair_dict_rev
                                        a = pair_dict_rev[a2]
                                        pair_dict[a].add(b2)
                                        pair_dict[b2] = a
                                    break

            # build_graph = [*new_rc]
            # did_any = False
            # if debug:
            #     logger.debug("starting with")
            #     logger.debug(new_rc)
            # if len(new_rc[i]) > 0:
            #     # logger.debug("Row", i)
            #     while len(pair_dict) > 0 and index >= 0:
            #         # col_index = max(new_rc[i]) - i - index - 1
            #         # logger.debug(pair_dict)
            #         col_index = index
            #         refl = col_index + i + 1
            #         # index = col_index + 1 + 1
            #         # assert index != 0 or new_rc.has_element(i + 1, col_index + 1)
            #         if debug:
            #             logger.debug("Starting")
            #             logger.debug(new_rc)
            #         if new_rc.has_element(i + 1, col_index + 1):
            #             # logger.debug(f"Found at {col_index + 1}")
            #             root = new_rc.right_root_at(i + 1, col_index + 1)
            #             # logger.debug(f"{root=}")
            #             if is_relevant_crossing(root, new_rc.perm):
            #                 # logger.debug("Relevant")
            #                 did_any = True
            #                 if root[0] in pair_dict:
            #                     pair_dict[root[0]].remove(root[1])
            #                     del pair_dict_rev[root[1]]
            #                     if len(pair_dict[root[0]]) == 0:
            #                         del pair_dict[root[0]]
            #                 # may have to remember this root
            #                 build_graph[i] = tuple(a for a in build_graph[i] if a != refl)
            #                 new_rc = RCGraph(build_graph)
            #                 if debug:
            #                     logger.debug(f"removed")
            #                     logger.debug(build_graph[i])
            #                     logger.debug(new_rc)
            #                 if new_rc.is_valid:
            #                     index -= 1
            #                     continue

            #                 for index2 in range(new_rc.max_of_row(i) - 2, index, -1):
            #                     # col_index2 = max(new_rc[i]) - i - index2 - 1
            #                     col_index2 = index2
            #                     refl = col_index2 + i + 1
            #                     if not new_rc.has_element(i + 1, col_index2 + 1):
            #                         a, b = new_rc.right_root_at(i + 1, col_index2 + 1)
            #                         root = (b, a)
            #                         if is_relevant_noncrossing(root):
            #                             if debug:
            #                                 logger.debug(f"Putting it back")
            #                             if a not in pair_dict:
            #                                 pair_dict[a] = set()
            #                             pair_dict[a].add(b)
            #                             pair_dict_rev[b] = a
            #                             val_to_insert = refl
            #                             new_row = new_rc[i]
            #                             if index == 0:
            #                                 new_row = [*new_row, val_to_insert]
            #                             else:
            #                                 for j in range(len(new_row)):
            #                                     if new_row[j] > val_to_insert:
            #                                         continue
            #                                     new_row = [*new_row[:j], val_to_insert, *new_row[j:]]
            #                                     break
            #                             build_graph[i] = tuple(new_row)
            #                             new_rc = RCGraph(build_graph)
            #                             if debug:
            #                                 logger.debug("added back")
            #                                 logger.debug(new_rc)
            #                             did_any = False
            #         index -= 1
            #         if did_any:
            #             break

            # # logger.debug(f"{did_any=} {index=}")
            # if not did_any:
            #     i -= 1
            # if i > 0:
            #     # logger.debug(f"{i=}")
            #     index = new_rc.max_of_row(i) - 1
            #     # logger.debug(f"{i=}")
            #     # logger.debug(f"{index=}")
        assert len(pair_dict_rev) == 0, f"{pair_dict=}, {pair_dict_rev=}, {working_rc=}"
        # logger.debug("Got")
        # logger.debug(new_rc)
        if debug:
            print("Finished with")
            print(working_rc)
            print(working_rc.perm)
        return working_rc

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

        if debug:
            logger.debug("Zeroing out last row")
            logger.debug(self)
            logger.debug("-------")

        diff_rows = []
        descs = []

        while len(interim.perm.trimcode) > row - 1:
            prev_interim = interim
            descs += [max(interim.perm.descents()) + 1]
            interim = interim.exchange_property(max(interim.perm.descents()) + 1)
            for i in range(len(interim)):
                if len(interim[i]) < len(prev_interim[i]):
                    rw = i + 1
                    diff_rows.append(rw)
                    break

        # go up to len interim - 1
        if debug:
            logger.debug("Transformed")
            logger.debug(interim)
            logger.debug("=========")

        interim = interim.kogan_kumar_insert(row, diff_rows, debug=debug)

        if debug:
            logger.debug("Got")
            logger.debug(interim)
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

    def _kogan_kumar_insert_row(self, row, descent, dict_by_a, dict_by_b, num_times, debug=False, start_index=0):
        working_rc = self
        if row > descent:
            raise ValueError("All rows must be less than or equal to descent")

        i = start_index

        num_done = 0

        while num_done < num_times:
            i += 1
            flag = False
            if debug:
                logger.debug(f"Trying column {i=} {descent=} {row=} {num_done=} {num_times=}")
            if not working_rc.has_element(row, i):
                a, b = working_rc.right_root_at(row, i)
                if a < b:
                    if debug:
                        logger.debug(f"root is {a, b}")
                    flag = False
                    if debug:
                        logger.debug(f"{dict_by_b=}")
                    if _is_row_root(descent, (a, b)) and b not in dict_by_b:
                        new_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[a] = dict_by_a.get(a, set())
                        dict_by_a[a].add(b)
                        dict_by_b[b] = a
                        flag = True
                        working_rc = new_rc
                        if debug:
                            logger.debug("Toggled a")
                            logger.debug(working_rc)
                    elif a in dict_by_b and b not in dict_by_b:
                        new_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[dict_by_b[a]].add(b)
                        dict_by_b[b] = dict_by_b[a]
                        flag = True
                        working_rc = new_rc
                        if debug:
                            logger.debug("Toggled b")
                            logger.debug(working_rc)
                    elif b in dict_by_b and a not in dict_by_b and a > descent:
                        new_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[dict_by_b[b]].add(a)
                        dict_by_b[a] = dict_by_b[b]
                        flag = True
                        working_rc = new_rc
                        if debug:
                            logger.debug("Toggled c")
                            logger.debug(working_rc)
                if flag:
                    num_done += 1
                    # assert last_rc.perm.inv + 1 == working_rc.perm.inv
                    if debug:
                        logger.debug("Inserted")
                        logger.debug(working_rc)
                if not working_rc.is_valid:
                    working_rc = working_rc._kogan_kumar_rectify(row - 1, descent, dict_by_a, dict_by_b)
        return working_rc

    def _kogan_kumar_rectify(self, row_below, descent, dict_by_a, dict_by_b):
        # logger.debug("In rectify")
        working_rc = self
        debug = False
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
                # logger.debug("Has element")
                a, b = working_rc.right_root_at(row_below, j)
                # logger.debug("root=", (a, b))
                if a > b:
                    # logger.debug("Entered")
                    if debug:
                        logger.debug(f"Considering bad at {row_below, j}")
                        logger.debug(f"{dict_by_a=}, {dict_by_b=}")
                        logger.debug(f"root = ({a, b})")
                    if b in dict_by_a and a in dict_by_a[b]:
                        new_rc = working_rc.toggle_ref_at(row_below, j)
                        dict_by_a[b].remove(a)
                        if len(dict_by_a[b]) == 0:
                            del dict_by_a[b]
                        del dict_by_b[a]
                        working_rc = new_rc
                        flag = True
                        # logger.debug("Toggle bad a")
                        # logger.debug(working_rc)
                    elif a in dict_by_b and b in dict_by_b and dict_by_b[a] == dict_by_b[b]:
                        new_rc = working_rc.toggle_ref_at(row_below, j)
                        if new_rc.perm[dict_by_b[a] - 1] < new_rc.perm[a - 1]:
                            dict_by_a[dict_by_b[a]].remove(a)
                            del dict_by_b[a]
                            if len(dict_by_a[dict_by_b[b]]) == 0:
                                del dict_by_a[dict_by_b[b]]
                            # logger.debug("Toggle bad b")
                            flag = True
                        dict_by_a[dict_by_b[b]].remove(b)
                        del dict_by_b[b]
                        if len(dict_by_a[dict_by_b[a]]) == 0:
                            del dict_by_a[dict_by_b[a]]
                        # logger.debug("Toggle bad c")
                if flag:
                    working_rc = working_rc._kogan_kumar_insert_row(row_below, descent, dict_by_a, dict_by_b, num_times=1, debug=debug)
        return working_rc._kogan_kumar_rectify(row_below - 1, descent, dict_by_a, dict_by_b)

    def associative_kogan_kumar_insert(self, descent, rows, debug=False):
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

    def kogan_kumar_insert(self, descent, rows, debug=False):
        dict_by_a = {}
        dict_by_b = {}
        # row is descent
        # inserting times

        working_rc = RCGraph([*self])
        if len(rows) == 0:
            return self
        rows_grouping = {}
        total_num = 0
        for r in rows:
            rows_grouping[r] = rows_grouping.get(r, 0) + 1
        if max(rows) > len(working_rc):
            working_rc = working_rc.extend(max(rows) - len(working_rc))
        rows = sorted(rows, reverse=True)
        if debug:
            logger.debug(f"inserting {rows=}")
        for row in sorted(rows_grouping.keys(), reverse=True):
            num_times = rows_grouping[row]
            if debug:
                logger.debug(f"Inserting {row=} {num_times=}")
                logger.debug(working_rc)
                logger.debug(f"{working_rc.perm.inv=}, {self.perm.inv=}")
            last_working_rc = working_rc
            working_rc = working_rc._kogan_kumar_insert_row(row, descent, dict_by_a, dict_by_b, num_times, debug=debug)

            if not working_rc.is_valid:
                working_rc = working_rc._kogan_kumar_rectify(row, descent, dict_by_a, dict_by_b)  # minus one?
            if debug:
                logger.debug("Next iteration")
                logger.debug(working_rc)
                # logger.debug(f"{working_rc.perm.inv=}, {self.perm.inv + index + 1=}")

            try:
                assert len(working_rc[row - 1]) == len(last_working_rc[row - 1]) + num_times
            except AssertionError:
                logger.debug("Assertion failed")
                logger.debug(working_rc)
                logger.debug(f"{working_rc.perm.inv=}, {self.perm.inv + total_num=}")
                logger.debug(f"{dict_by_a=}, {dict_by_b=}")
                logger.debug(f"{working_rc.perm=}, {self.perm=}")
                if debug:
                    raise
                self.kogan_kumar_insert(descent, rows, debug=False)
                raise
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
    
    def zero_out_last_row(self, debug=False):
        # this is important!
        # transition formula
        if len(self[-1]) != 0:
            raise ValueError("Last row not empty")
        if self.perm.inv == 0:
            return self.rowrange(0, len(self) - 1)

        interim = RCGraph([*self])
        diff_rows = []

        if debug:
            logger.debug("Zeroing out last row")
            logger.debug(self)
            logger.debug("-------")

        diff_rows = []
        descs = []
        while len(interim.perm.trimcode) > len(self) - 1:
            prev_interim = interim
            descs += [max(interim.perm.descents()) + 1]
            interim = interim.exchange_property(max(interim.perm.descents()) + 1)
            for i in range(len(interim)):
                if len(interim[i]) < len(prev_interim[i]):
                    rw = i + 1
                    diff_rows.append(rw)
                    break

        # go up to len interim2 - 1
        interim2 = RCGraph([*interim[:-1], tuple(sorted(descs, reverse=True))])
        if debug:
            logger.debug("Transformed")
            logger.debug(interim2)
            logger.debug("=========")

        interim = interim2.kogan_kumar_insert(len(self), diff_rows, debug=debug)

        if debug:
            logger.debug("Got")
            logger.debug(interim)
        assert interim.length_vector[:-1] == self.length_vector[:-1]
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
                pos_list = [i for i in range(k) if up_perm[i] < last]
                for j in range(k, max(k + 2, len(up_perm) + 1)):
                    if up_perm[j] >= last:
                        continue
                    for i in pos_list:
                        if has_bruhat_descent(up_perm, i, j):
                            new_perm_add = up_perm.swap(i, j)
                            perm_list += [(new_perm_add, up_perm[j])]
                            new_total_dict[new_perm_add] = [(i + 1, j + 1)] + total_dict[up_perm]
            up_perm_list = perm_list
            total_dict = new_total_dict
        return total_dict

    def vertical_cut(self,row):
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

    def right_zero_act(self, debug=False):
        if self.perm.inv == 0:
            return {RCGraph([*self, ()])}
        if debug:
            logger.debug("Num rows:")
            logger.debug(len(self))
            logger.debug(self)



        up_perms = ASx(self.perm, len(self)) * ASx(uncode([0]), 1)

        rc_set = set()
        
        # cannot_do = 0
        # padded_weight = (*self.length_vector, 0)
        print("Self")
        print(self)
        for (perm, _), v in up_perms.items():
            assert v == 1
        #     # done_any = False
            print(f"Looking for {perm}")
            pull_graph = self.extend(1)
            #num_d = len(self.perm) + 1 - len(self)
            num_d = 0
            nperm = perm
            descs = []
            while nperm.inv > 0 and len(nperm.trimcode) > len(self):
                desc = len(nperm.trimcode)
                nperm = nperm.swap(desc - 1, desc)
                descs.append(desc)
                num_d += 1
            print(f"{descs=}")
            # print(f"{num_d=}")
            pull_graph = RCGraph([*pull_graph[:-1], range(len(self) + num_d, len(self), -1)])
            # print(pull_graph)
            ref_dict = RCGraph.complete_sym_perms_op(pull_graph.perm, num_d, len(self) + 1)
            if perm not in ref_dict:
                continue
            # print(f"ref_dict[perm]={ref_dict[perm]}")
            the_graph = pull_graph.reverse_kogan_kumar_insert(len(self) + 1, ref_dict[perm])#.rowrange(0,len(self)).extend(1)
            for d in descs:
                the_graph = the_graph.exchange_property(d)
            # now a different to put the descents back
            # exchange
            print("Afterward")
            print(the_graph)
            diff_rows = []
            for rw in range(len(self)):
                if len(the_graph[rw]) != len(self[rw]):
                    diff_rows = diff_rows + ([rw + 1] * (len(self[rw]) - len(the_graph[rw])))
            # print(f"{diff_rows=}")
            for d in reversed(descs):
                if debug:
                    logger.debug(f"Putting back descent at {d=}")
                for i, row in enumerate(the_graph[:-1]):
                    if len(row) == len(self[i]):
                        continue
                    found = False
                    for j in range(max(d + 1,self.cols)):
                        # print(f"{a=}, {j=}, {the_graph[a-1,j]=}")
                        if not the_graph[i, j]:
                            #the_graph = the_graph.toggle_ref_at(a, j + 1)
                            aa, b = the_graph.right_root_at(i + 1, j + 1)
                            # print(f"root {aa,b}")
                            if (aa, b) == (d, d + 1):
                                # print(f"Putting it in row {a} {j+1}")
                                # print("Got")
                                # print(the_graph)
                                found = True
                                the_graph = the_graph.toggle_ref_at(i + 1, j + 1)
                                break

                if not found:
                    raise ValueError("Could not put back descent")
                    #diff_rows.remove(a)
            # try:
            #     assert len(diff_rows) == 0, f"{diff_rows=}, {the_graph=}, {self=}, {perm=}"
            # except AssertionError:
            #     logger.debug("Could not put back descents")
            #     logger.debug(the_graph)
            #     logger.debug(f"{diff_rows=}, {the_graph=}, {self=}, {perm=}")
            #     continue
            # # take out the exchange
            # refs = [(d, len(self) + 2) for d in range(len(self) + 2 - num_d, len(self) + 2)]
            # print(f"refs={refs}, diff_rows={diff_rows}")
            # the_graph = the_graph.reverse_kogan_kumar_insert(len(self) + 1, refs)#.rowrange(0,len(self)).extend(1)
            # for boggle in diff_rows:
            #     the_graph = the_graph.kogan_kumar_insert(len(self) + 1, [boggle])
            # print("Final")
            # print(the_graph)
            assert the_graph.perm == perm, f"{the_graph.perm=}, {perm=}, {the_graph=}, {self=}"
            #assert the_graph.rowrange(0, len(the_graph) - 1).perm == self.perm
            #assert the_graph.perm == perm
            assert the_graph.zero_out_last_row() == self, f"{the_graph.zero_out_last_row()=}, {self=}, {the_graph=}, {perm=}"
            rc_set.add(the_graph)

            # for rc in RCGraph.all_rc_graphs(perm, len(self) + 1, weight=padded_weight):
            #     if debug:
            #         logger.debug("rc num rows:")
            #         logger.debug(len(rc))
            #     assert len(rc) == len(self) + 1
            #     z = rc.zero_out_last_row()
            #     assert len(z) == len(self)
            #     if debug:
            #         logger.debug("Zeroed out last row:")
            #         logger.debug(z)
            #     if z == self:
            #         if debug:
            #             logger.debug("Added")
            #         rc_set.add(rc)
        #     if not done_any:
        #         cannot_do += 1

        # try:
        #     assert len(rc_set) + cannot_do == len(up_perms), f"{rc_set=}, {len(up_perms)=}, {self=} {up_perms=} {cannot_do=}"
        # except AssertionError:
        #     logger.debug("Could not achieve some")
        #     lst = [perm[0] for perm in up_perms.keys() if perm[0] not in [rc.perm for rc in rc_set]]
        #     logger.debug(f"Missing {lst=}")
        #     logger.debug("Found:")
        #     for rc in rc_set:
        #         logger.debug(f"Perm: {rc.perm}")
        #         logger.debug(rc)
        #         logger.debug("===========")

        #     logger.debug(f"All rc's missing {len(self)=}")
        #     for perm in lst:
        #         for rc in RCGraph.all_rc_graphs(perm, len(self) + 1):
        #             if rc.length_vector[:-1] == self.length_vector and rc not in rc_set:
        #                 logger.debug(rc)
        #                 logger.debug("----------")
        #                 logger.debug("Zeroed")
        #                 logger.debug(rc.zero_out_last_row())
        #                 logger.debug("===========")
        #     raise
        return rc_set

    @cache
    def bisect_left_coords_index(self, row, col):
        from bisect import bisect_left

        inversions = [self.left_to_right_inversion_coord(i) for i in range(len(self.perm_word))]
        inversions.sort(key=lambda x: (x[0], -x[1]))
        # logger.debug(f"{inversions=}")

        return bisect_left(inversions, (row, col), key=lambda x: (x[0], -x[1]))

    def exchange_property(self, descent, left=False):
        for i in range(len(self.perm_word) + 1):
            if not left:
                a, b = self.left_to_right_inversion(i)
            else:
                a, b = self.left_to_right_left_inversion(i)
            if a == descent and b == descent + 1:
                return self.toggle_ref_at(*self.left_to_right_inversion_coord(i))
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
        row = 0
        sm = 0
        for i in range(len(self)):
            sm += len(self[i])
            if sm >= index and len(self[i]) > 0:
                row = i + 1
                break
        along = sm - index
        col_index = len(self[row - 1]) - 1 - along
        # logger.debug(f"{index=} {sm=} {along=} {len(self[row - 1])=}, {col_index=}")
        col = self[row - 1][col_index] - row + 1
        return (row, col)

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
            buildup_module = new_buildup_module

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
                # logger.debug(f"{k=} {v=}")
                acted_element = {self: v}
                for a in reversed(k):
                    acted_element2 = {}
                    for k2, v2 in acted_element.items():
                        # logger.debug(f"{dict.fromkeys(k2.act(a), v2)=}")
                        acted_element2 = add_perm_dict(acted_element2, dict.fromkeys(k2.act(a), v2))
                    acted_element = acted_element2
                # logger.debug(f"{acted_element=}")
                ret = add_perm_dict(ret, acted_element)
            # logger.debug(f"{ret=}")
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
