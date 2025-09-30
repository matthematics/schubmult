import logging
from functools import cache, cached_property
from itertools import zip_longest

from symengine import SympifyError

import schubmult.schub_lib.schub_lib as schub_lib
from schubmult.perm_lib import Permutation, uncode
from schubmult.rings import ASx
from schubmult.symbolic import DefaultPrinting, S, prod, sympify
from schubmult.utils.perm_utils import add_perm_dict

from .free_algebra import FreeAlgebra, FreeAlgebraElement
from .free_algebra_basis import SchubertBasis, WordBasis
from .nil_hecke import NilHeckeRing
from .tensor_ring import TensorRing, TensorRingElement

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def _is_row_root(row, root):
    return root[0] <= row and root[1] > row


ASx = FreeAlgebra(SchubertBasis)
FA = FreeAlgebra(WordBasis)
FAS = FA


class UnderlyingGraph(tuple):
    def __new__(cls, *args, value=1):
        new_args = tuple(tuple(arg) for arg in args)
        return UnderlyingGraph.__xnew_cached__(cls, *new_args)

    @staticmethod
    @cache
    def __xnew_cached__(_class, *args):
        return UnderlyingGraph.__xnew__(_class, *args)

    @staticmethod
    def __xnew__(_class, *args):
        return tuple.__new__(_class, *args)


class RCGraph(UnderlyingGraph, DefaultPrinting):
    # def __str__(self):
    #     return self._sympystr()

    def __eq__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return tuple(self) == tuple(other)

    def _toggle_root(self, row, col, roots, used_bs, root_set):
        root = self.right_root_at(row, col)
        # ret = None
        if root in root_set:
            roots[root[0]] = roots.get(root[0], [])
            roots[root[0]].append(root[1])
            used_bs[root[1]] = root[0]
            root_set.remove(root)
            working_rc = self.toggle_ref_at(row, col)
            return working_rc, not working_rc.is_valid
        if root[0] in used_bs and (used_bs[root[0]], root[1]) in root_set:
            a = used_bs[root[0]]
            root_set.remove((used_bs[root[0]], root[1]))
            roots[a].insert(roots[a].index(root[0]), root[1])
            used_bs[root[1]] = a
            working_rc = self.toggle_ref_at(row, col)
            return working_rc, not working_rc.is_valid
        if root[1] in used_bs and (used_bs[root[1]], root[0]) in root_set:
            a = used_bs[root[1]]
            root_set.remove((used_bs[root[1]], root[0]))
            roots[a].insert(roots[a].index(root[1]), root[0])
            used_bs[root[0]] = a
            working_rc = self.toggle_ref_at(row, col)
            return working_rc, not working_rc.is_valid
        return self, False

    def _rectify(self, row, row_below, roots, used_bs, root_set):
        working_rc = self
        if working_rc.is_valid:
            return working_rc
        if row_below <= 0:
            return self
        for j in range(self.max_of_row(row_below) + 1, 0, -1):
            reinsert = False
            if self.has_element(row_below, j):
                b, a = self.right_root_at(row_below, j)
                # logger.debug(f"Considering {a, b} at {row_below, j}")
                # logger.debug(f"{roots=}, {used_bs=}")
                if a in roots and b in roots[a]:
                    roots[a].remove(b)
                    if len(roots[a]) == 0:
                        del roots[a]
                    del used_bs[b]
                    root_set.add((a, b))
                    working_rc = self.toggle_ref_at(row_below, j)
                elif a in used_bs and b in used_bs and used_bs[a] == used_bs[b]:
                    working_rc = self.toggle_ref_at(row_below, j)
                    root_set.add((a, b))
                    if working_rc.perm[used_bs[a] - 1] < working_rc.perm[a - 1]:
                        roots[used_bs[a]].remove(a)
                        del used_bs[a]
                        if len(roots[used_bs[b]]) == 0:
                            del roots[used_bs[b]]
                    else:
                        roots[used_bs[b]].remove(b)
                        del used_bs[b]
                        if len(roots[used_bs[a]]) == 0:
                            del roots[used_bs[a]]
            if not working_rc.is_valid:
                for jp in range(working_rc.max_of_row(row_below) + 1, 0, -1):
                    if not working_rc.has_element(row_below, jp):
                        b2, a2 = working_rc.right_root_at(row_below, jp)
                        if (a2, b2) in root_set and b2 not in used_bs:
                            working_rc = working_rc.toggle_ref_at(row_below, jp)
                            root_set.remove(a2, b2)
                            roots[a2] = roots.get(a2, [])
                            roots[a2].append(b2)
                            used_bs[b2] = a2
            else:
                return working_rc
        # logger.debug("After rectify row ", row_below)
        # logger.debug(working_rc)
        return working_rc._rectify(row, row_below - 1, roots, used_bs)

    @property
    def is_valid(self):
        if self.perm.inv != len(self.perm_word()):
            return False
        # if max(self.perm.descents()) + 1 > len(self):
        #     return False
        return True

    # # kogan
    # def fill_row_with_boxes(self, row, roots_seq):
    #     roots = {}
    #     used_bs = {}
    #     for root in roots_seq:
    #         roots[root[0]] = roots.get(root[0], [])
    #         roots[root[0]].append(root[1])
    #         used_bs[root[1]] = root[0]

    #     working_rc = self
    #     while len(working_rc.perm_word()) < len(self.perm_word()) + num_boxes:
    #         for col in range(working_rc.max_of_row(row) + 1, 0, -1):
    #             if not working_rc.has_element(row, col):
    #                 working_rc, did = working_rc._toggle_root(row, col, roots, used_bs)
    #                 if did:
    #                     working_rc = working_rc._rectify(row, row - 1, roots, used_bs)
    #                     assert working_rc.is_valid
    #                     break
    #     return working_rc, roots

    # def insert_reflections(self, row0, roots_seq):
    #     roots = {}
    #     used_bs = {}
    #     roots = {}
    #     used_bs = {}
    #     # for root in roots_seq:
    #     #     roots[root[0]] = roots.get(root[0], [])
    #     #     roots[root[0]].append(root[1])
    #     #     used_bs[root[1]] = root[0]
    #     roots_set = set(roots_seq)
    #     working_rc = self
    #     while len(roots_set) > 0:
    #         for row in range(row0, 0, -1):
    #             for col in range(working_rc.max_of_row(row) + 1, 0, -1):
    #                 if not working_rc.has_element(row, col):
    #                     working_rc, did = working_rc._toggle_root(row, col, roots, used_bs, roots_set)
    #                     # logger.debug(working_rc)
    #                     if did:
    #                         working_rc = working_rc._rectify(row0, row, roots, used_bs, roots_set)
    #                         assert working_rc.is_valid
    #                         break
    #         # logger.debug(roots_set)
    #     return working_rc

    # def _rectify_remove(self, row, row_below, roots, used_bs):
    #     working_rc = self
    #     if row_below <= 0:
    #         return self
    #     working_rc = self
    #     for j in range(working_rc.max_of_row(row_below), 0, -1):
    #         if working_rc.has_element(row_below, j):
    #             a, b = working_rc.right_root_at(row_below, j)
    #             if b < a:
    #                 working_rc = working_rc.toggle_ref_at(row_below, j)
    #                 for j2 in range(self.max_of_row(row_below), 0, -1):
    #                     if not working_rc.has_element(row_below, j2) and _is_row_root(row, working_rc.right_root_at(row_below, j2)):
    #                         working_rc = working_rc.toggle_ref_at(row_below, j2)
    #                         break
    #     return working_rc

    # # kogan
    # def remove_boxes_from_row(self, row, indexes):
    #     indexes = sorted(indexes, reverse=True)
    #     working_rc = self
    #     for col in range(working_rc.max_of_row(row) + 1, 0, -1):
    #         if working_rc.has_element(row, col) and col not in indexes:
    #             working_rc = working_rc.toggle_ref_at(row, col)
    #     return working_rc

    def shiftup(self, shift):
        return [tuple([a + shift for a in rrow]) for rrow in self]

    def insert_boxes_into_row(self, row, indexes):
        """
        Insert reflections at the indexes in the given row
        """
        # check if valid
        # logger.debug(f"{indexes=}")
        indexes = sorted(indexes, reverse=True)
        working_rc = self
        while len(working_rc) >= row:
            working_rc = working_rc.normalize_and_remove_last_row()
        return RCGraph([*working_rc, tuple([a + row for a in indexes]), *self.rowrange(row - 1, len(self)).shiftup(row)])

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

    def reverse_kogan_kumar_insert(self, descent, reflection_path):
        # logger.debug("Here")
        from schubmult.utils.perm_utils import has_bruhat_descent

        # pair_sequence = sorted(pair_sequence, key=lambda x: x[0]
        pair_dict = reflection_path
        pair_dict_rev = {}
        # ref_by_index = {}
        working_perm = self.perm
        for a, b_list in pair_dict.items():
            for b in b_list:
                assert a <= descent and descent < b
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
            top, bottom = max(root), min(root)
            return (root[0] <= descent and descent < root[1] and root[1] not in pair_dict_rev) or (
                (bottom in pair_dict_rev or (bottom not in pair_dict_rev and top in pair_dict_rev)) and bottom > descent
            )

        # Add this intersection. If we are in the first case, insert (s, q) into the sequence (ai, bi) in the rightmost position, such that ai’s remain nondecreasing in the sequence. ((s, q) are the rows where the two strands shown in Figure 3 originate.) If
        # we are in the second case, add (ai, q) just before where (a, bi) is in the sequence.

        new_rc = self
        i = descent - 1
        index = new_rc.max_of_row(i) - 2 if len(new_rc[i]) > 0 else -1
        # logger.debug(pair_dict)

        while i >= 0:
            build_graph = [*new_rc]
            did_any = False
            # logger.debug("starting with")
            # logger.debug(new_rc)
            if len(new_rc[i]) > 0:
                # logger.debug("Row", i)
                while len(pair_dict) > 0 and index >= 0:
                    # col_index = max(new_rc[i]) - i - index - 1
                    # logger.debug(pair_dict)
                    col_index = index
                    refl = col_index + i + 1
                    # index = col_index + 1 + 1
                    # assert index != 0 or new_rc.has_element(i + 1, col_index + 1)
                    # logger.debug("Starting")
                    # logger.debug(new_rc)
                    if new_rc.has_element(i + 1, col_index + 1):
                        # logger.debug(f"Found at {col_index + 1}")
                        root = new_rc.right_root_at(i + 1, col_index + 1)
                        # logger.debug(f"{root=}")
                        if is_relevant_crossing(root, new_rc.perm):
                            # logger.debug("Relevant")
                            did_any = True
                            if root[0] in pair_dict:
                                pair_dict[root[0]].remove(root[1])
                                del pair_dict_rev[root[1]]
                                if len(pair_dict[root[0]]) == 0:
                                    del pair_dict[root[0]]
                            # may have to remember this root
                            build_graph[i] = tuple(a for a in build_graph[i] if a != refl)
                            new_rc = RCGraph(build_graph)
                            # logger.debug(f"removed")
                            # logger.debug(build_graph[i])
                            # logger.debug(new_rc)
                            if new_rc.is_valid:
                                index -= 1
                                continue

                            for index2 in range(new_rc.max_of_row(i) - 2, index, -1):
                                # col_index2 = max(new_rc[i]) - i - index2 - 1
                                col_index2 = index2
                                refl = col_index2 + i + 1
                                if not new_rc.has_element(i + 1, col_index2 + 1):
                                    a, b = new_rc.right_root_at(i + 1, col_index2 + 1)
                                    root = (b, a)
                                    if is_relevant_noncrossing(root):
                                        # logger.debug(f"Putting it back")
                                        if a not in pair_dict:
                                            pair_dict[a] = set()
                                        pair_dict[a].add(b)
                                        pair_dict_rev[b] = a
                                        val_to_insert = refl
                                        new_row = new_rc[i]
                                        if index == 0:
                                            new_row = [*new_row, val_to_insert]
                                        else:
                                            for j in range(len(new_row)):
                                                if new_row[j] > val_to_insert:
                                                    continue
                                                new_row = [*new_row[:j], val_to_insert, *new_row[j:]]
                                                break
                                        build_graph[i] = tuple(new_row)
                                        new_rc = RCGraph(build_graph)
                                        # logger.debug("added back")
                                        # logger.debug(new_rc)
                                        did_any = False
                    index -= 1
                    if did_any:
                        break

            # logger.debug(f"{did_any=} {index=}")
            if not did_any:
                i -= 1
            if i > 0:
                # logger.debug(f"{i=}")
                index = new_rc.max_of_row(i) - 1
                # logger.debug(f"{i=}")
                # logger.debug(f"{index=}")
        assert len(pair_dict_rev) == 0, f"{pair_dict=}, {pair_dict_rev=}, {new_rc=}"
        # logger.debug("Got")
        # logger.debug(new_rc)
        return new_rc

    def right_p_act(self, p):
        if len(self) == 0:
            return set((FA(p) * self).value_dict.keys())

        up_perms = ASx(self.perm, len(self)) * ASx(uncode([p]), 1)

        old_rc_set = self.rowrange(1, len(self)).right_p_act(p)
        rc_set = set()

        for rc in old_rc_set:
            new_rc_set = FA(len(self[0])) * rc
            for rc2 in new_rc_set.value_dict.keys():
                if (rc2.perm, len(self) + 1) in up_perms.keys():
                    rc_set.add(rc2)
        return rc_set

    @cache
    def inversion_label(self, i, j):
        if i >= j:
            raise ValueError("i must be less than j")
        if self.perm[i] < self.perm[j]:
            raise ValueError("Not an inversion")
        index = 0
        for i0, row in enumerate(self):
            for j0, a in enumerate(row):
                if self.perm.right_root_at(index, word=self.perm_word()) == (i + 1, j + 1):
                    return i0 + 1
                index += 1
        raise ValueError("Could not find inversion")

    @cache
    def lehmer_label(self, i, j):
        value = self.inversion_label(i, j)
        numeros = set(list(range(1, value + 1)))
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

    def __new__(cls, *args, value=1):
        rcc = UnderlyingGraph.__new__(cls, *args)
        rcc.value = value
        return rcc

    def __init__(self, *args, value=1):
        pass

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
        perm_word = self.perm_word()
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

    def has_element(self, i, j):
        return i <= len(self) and j + i - 1 in self[i - 1]

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

    _graph_cache = {}

    @classmethod
    def all_rc_graphs(cls, perm, length=-1):
        if length == -1:
            length = len(perm.trimcode)
        if (perm, length) in cls._graph_cache:
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
        for index_list, new_perm in L:
            new_row = [new_perm[i] for i in range(max(len(pm), len(new_perm))) if new_perm[i] == pm[i + 1]]
            new_row.sort(reverse=True)
            oldset = cls.all_rc_graphs(new_perm)
            for old_rc in oldset:
                nrc = RCGraph([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in old_rc]])
                if len(nrc) < length:
                    nrc = RCGraph((*nrc, *tuple([()] * (length - len(nrc)))))
                assert nrc.perm == perm
                assert len(nrc) == length
                ret.add(nrc)
        cls._graph_cache[(perm, length)] = ret
        return ret

    def extend(self, extra_rows):
        return RCGraph([*self, *tuple([()] * extra_rows)])

    # def __le__(self, other):
    #     if not isinstance(other, RCGraph):
    #         return NotImplemented
    #     if not isinstance(other, RCGraph):
    #         return NotImplemented
    #     if self.length_vector() < other.length_vector():
    #         return True
    #     if self.length_vector() == other.length_vector() and self.perm.bruhat_leq(other.perm) and self.perm != other.perm:
    #         return True
    #     if self.length_vector() == other.length_vector() and self.perm == other.perm and self.perm_word() < other.perm_word():
    #         return True
    #     return self.length_vector() == other.length_vector() and self.perm == other.perm and self.perm_word() == other.perm_word()

    # def __lt__(self, other):
    #     if not isinstance(other, RCGraph):
    #         return NotImplemented
    #     if self.length_vector() < other.length_vector():
    #         return True
    #     if self.length_vector() == other.length_vector() and self.perm.bruhat_leq(other.perm) and self.perm != other.perm:
    #         return True
    #     if self.length_vector() == other.length_vector() and self.perm == other.perm and self.perm_word() < other.perm_word():
    #         return True
    #     return False

    def _kogan_kumar_insert_row(self, row, descent, dict_by_a, dict_by_b, num_times, debug=False, start_index=0):
        # working_rc = self
        ARBITRARY_MAX_VALUE = 10
        working_rc = self
        if row > descent:
            raise ValueError("All rows must be less than or equal to descent")
        # for i in range(working_rc.max_of_row(row) + descent, 0, -1):
        i = start_index

        num_done = 0

        while num_done < num_times:
            last_rc = working_rc
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
                        logger.debug("_is_row_root:", _is_row_root(descent, (a, b)))
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
                self.kogan_kumar_insert(descent, rows, debug=True)
                raise
        # reflections = []
        # start_perm = self.perm
        # a_keys = sorted(dict_by_a.keys())
        # # logger.debug(dict_by_a)
        # # logger.debug(start_perm)
        # # input()
        # for a in a_keys:
        #     while a in dict_by_a.keys():
        #         for b in sorted(dict_by_b.keys(), reverse=True):
        #             if has_bruhat_ascent(start_perm, a - 1, b - 1):
        #                 reflections.append((a, b))
        #                 start_perm = start_perm.swap(a - 1, b - 1)
        #                 dict_by_a[a].remove(b)
        #                 if len(dict_by_a[a]) == 0:
        #                     del dict_by_a[a]
        #                 del dict_by_b[b]
        #                 break
        # logger.debug(f"Did not find {a,b} start_perm={start_perm} {dict_by_a=}")
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

    # def dualact(self, p):
    #     pm = ~self.perm
    #     elem = FAS(pm, len(self))
    #     bumpup = FAS(uncode([p]), 1) * elem
    #     ret = set()
    #     for k, v in bumpup.items():
    #         perm2 = k[0]
    #         new_row = [pm[i] for i in range(max(len(pm), len(perm2))) if pm[i] == perm2[i + 1]]
    #         new_row.sort()
    #         lst = [tuple([a + 1 for a in row]) for row in self]

    #         for index in range(max(len(self) + 1, *new_row)):
    #             if index < len(lst):
    #                 if index + 1 in new_row:
    #                     lst[index] = (*lst[index], index + 1)
    #             else:
    #                 if index + 1 in new_row:
    #                     lst += [(index + 1,)]
    #                 else:
    #                     lst += [()]
    #         nrc = RCGraph(lst)
    #         assert nrc.perm == ~perm2
    #         ret.add(nrc)
    #     return ret

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

    w_key_cache = {}
    rc_cache = set()

    @classmethod
    def fa_hom(cls, *word):
        if len(word) == 0:
            return cls()
        if len(word) == 1:
            return FA(*word) * cls()
        return cls.fa_hom(*word[:-1]) * cls.fa_hom(word[-1])

    @classmethod
    def fa_coprod(cls, *word):
        if len(word) == 0:
            return cls() @ cls()
        if all(a == 0 for a in word):
            return cls([() * len(word)]) @ cls([() * len(word)])
        cop = FA(*word).coproduct()
        result = 0
        for (word1, word2), coeff in cop.items():
            result += RCGraph.fa_hom(*word1) @ RCGraph.fa_hom(*word2)
        return result

    def toggle_ref_at(self, i, j):
        from bisect import bisect_left

        assert i > 0 and j > 0
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

    def monk_insert(self, descent, row):
        assert row <= descent
        result = None
        for i in range(max(len(self[row - 1]) - row + 1, descent), 0, -1):
            if not self.has_element(row, i) and _is_row_root(descent, self.right_root_at(row, i)):
                result = self.toggle_ref_at(row, i)
                break
        if result.is_valid:
            return result
        # rectify
        # logger.debug("Result")
        # logger.debug(result)

        return result._monk_rectify(descent, row)

    def _monk_rectify(self, descent, row_below):
        working_rc = self
        if working_rc.is_valid:
            return working_rc
        if row_below <= 0:
            # logger.debug("End result")
            # logger.debug(self)
            raise ValueError("Could not rectify")
        for j in range(working_rc.max_of_row(row_below) + 1, 0, -1):
            if working_rc.has_element(row_below, j):
                b, a = working_rc.right_root_at(row_below, j)
                if _is_row_root(descent, (a, b)):
                    working_rc = working_rc.toggle_ref_at(row_below, j)
                    for jp in range(working_rc.max_of_row(row_below) + 1, 0, -1):
                        if not working_rc.has_element(row_below, jp) and _is_row_root(descent, working_rc.right_root_at(row_below, jp)):
                            working_rc = working_rc.toggle_ref_at(row_below, jp)
                            if working_rc.is_valid:
                                return working_rc
                    return working_rc._monk_rectify(descent, row_below)
        return working_rc._monk_rectify(descent, row_below - 1)

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
        # if max(self.perm.descents()) + 1 <= len(self) - 1:
        #     return self.rowrange(0, len(self) - 1)
        # exchange property div diff sn
        interim = RCGraph([*self])
        diff_rows = []
        # we need to bump up additionally?
        # example
        #
        if debug:
            logger.debug("Zeroing out last row")
            logger.debug(self)
            logger.debug("-------")
        total = 0
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
        assert interim.length_vector()[:-1] == self.length_vector()[:-1]
        return interim.rowrange(0, len(self) - 1)

    def right_zero_act(self, debug=False):
        # logger.debug("Right zeroing")
        # logger.debug(self)
        if self.perm.inv == 0:
            return {RCGraph([*self, ()])}
        if debug:
            logger.debug("Num rows:")
            logger.debug(len(self))
            logger.debug(self)
        # if len(self) >= len(self.perm) - 1:
        #     return {RCGraph([*self, ()])}
        up_perms = ASx(self.perm, len(self)) * ASx(uncode([0]), 1)

        rc_set = set()

        cannot_do = 0
        for (perm, _), v in up_perms.items():
            assert v == 1
            done_any = False
            for rc in RCGraph.all_rc_graphs(perm, len(self) + 1):
                if debug:
                    logger.debug("rc num rows:")
                    logger.debug(len(rc))
                    assert len(rc) == len(self) + 1
                if rc.length_vector()[:-1] == self.length_vector():
                    done_any = True
                    z = rc.zero_out_last_row()
                    assert len(z) == len(self)
                    if debug:
                        logger.debug("Zeroed out last row:")
                        logger.debug(z)
                    if z == self:
                        if debug:
                            logger.debug("Added")
                        rc_set.add(rc)
            if not done_any:
                cannot_do += 1

        try:
            assert len(rc_set) + cannot_do == len(up_perms), f"{rc_set=}, {len(up_perms)=}, {self=} {up_perms=} {cannot_do=}"
        except AssertionError:
            logger.debug("Could not achieve some")
            lst = [perm[0] for perm in up_perms.keys() if perm[0] not in [rc.perm for rc in rc_set]]
            logger.debug(f"Missing {lst=}")
            logger.debug("Found:")
            for rc in rc_set:
                logger.debug(f"Perm: {rc.perm}")
                logger.debug(rc)
                logger.debug("===========")

            logger.debug(f"All rc's missing {len(self)=}")
            for perm in lst:
                for rc in RCGraph.all_rc_graphs(perm, len(self) + 1):
                    if rc.length_vector()[:-1] == self.length_vector() and rc not in rc_set:
                        logger.debug(rc)
                        logger.debug("----------")
                        logger.debug("Zeroed")
                        logger.debug(rc.zero_out_last_row())
                        logger.debug("===========")
            raise
        return rc_set

    @staticmethod
    def right_zero_act_pair(rc1, rc2):
        # logger.debug("Right zeroing")
        # logger.debug(self)
        if len(rc1) == 0:
            return {(RCGraph([()]), RCGraph([()]))}
        up_perms = (ASx @ ASx)(((rc1.perm, len(rc1)), (rc2.perm, len(rc2)))) * ASx(uncode([0]), 1).coproduct()

        rc_set = set()

        for ((perm1, _), (perm2, _)), v in up_perms.items():
            for rc01 in RCGraph.all_rc_graphs(perm1, len(rc1) + 1):
                for rc02 in RCGraph.all_rc_graphs(perm2, len(rc2) + 1):
                    # rc = RCGraph((tuple(sorted([*rc01[0], *rc02[0]])), *[tuple(sorted([*rc01[i], *[a + 1 for a in rc02[i]]])) for i in range(1, max(len(rc01), len(rc02)))]))
                    if rc02.length_vector()[:1] == rc2.length_vector() and rc01.length_vector()[:-1] == rc1.length_vector() and rc01.zero_out_last_row() == rc1 and rc02.zero_out_last_row() == rc2:
                        rc_set.add((rc01, rc02))
        return rc_set

    def normalize_and_remove_last_row(self):
        if len(self) == 0:
            raise ValueError("No rows to remove")
        if self.perm.inv == 0:
            return RCGraph(self[:-1])
        working_rc = RCGraph(*self)
        # logger.debug("working_rc")
        # logger.debug(working_rc)
        # logger.debug("---------")
        # logger.debug(f"{len(self)=}")
        # logger.debug(f"{max(self.perm.descents())=}")
        # input()
        while working_rc.perm[len(self) - 1] != len(working_rc.perm):
            # logger.debug("Down")
            working_rc = working_rc.fill_row_with_boxes(len(working_rc), 1)[0]
            # logger.debug("working_rc")
            # logger.debug(working_rc)
            # logger.debug("---------")
            # logger.debug("cut")
            # logger.debug(working_rc.rowrange(0,len(self)-1))
            # logger.debug("---------")
            # logger.debug("self")
            # logger.debug(self)
            # logger.debug(f"{len(self) - 1=}")
            # logger.debug(f"{max(working_rc.rowrange(0,len(self)-1).perm.descents())=} >= {max(self.perm.descents())=}")
            # input()
        reflections = {}

        reflections[len(working_rc)] = [len(working_rc) + a for a in range(len(working_rc[-1]) - 1, 0, -1)]
        # logger.debug(f"{reflections=}")
        working_rc = RCGraph(working_rc.reverse_kogan_kumar_insert(len(working_rc), reflections)[:-1])
        return RCGraph(working_rc)
        # logger.debug(working_rc)

    # certificate to insert the row

    def bisect_left_coords_index(self, row, col):
        from bisect import bisect_left

        inversions = [self.left_to_right_inversion_coord(i) for i in range(len(self.perm_word()))]
        inversions.sort(key=lambda x: (x[0], -x[1]))
        # logger.debug(f"{inversions=}")

        return bisect_left(inversions, (row, col), key=lambda x: (x[0], -x[1]))

    # def find_bruhat_embedding(self, perm, skip_row):
    #     working_perm = perm
    #     working_self_perm = self.perm
    #     word = self.perm_word()
    #     build_word = []
    #     for a in reversed(word):
    #         row, col = self.left_to_right_inversion_coord(len(word) - len(build_word) - 1)
    #         if working_perm[a-1] > working_perm[a]:
    #             build_word.append(a)
    #             working_perm = working_perm.swap(a - 1, a)
    #             working_self_perm = working_self_perm.swap(a - 1, a)
    #         else:
    #             if row != skip_row:
    #                 # want to stay between
    #             else:
    #                 desc = next(iter([b for b in working_perm.descents() if b not in working_self_perm.descents()]))
    #                 build_word.append(desc + 1)
    #                 working_perm = working_perm.swap(desc, desc + 1)

    def exchange_property(self, descent):
        for i in range(len(self.perm_word()) + 1):
            a, b = self.left_to_right_inversion(i)
            # logger.debug(f"{descent=}, {a=}, {b=}")
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

    def extract_row2(self, row):
        from schubmult.schub_lib.schub_lib import pull_out_var

        working_rc = self
        to_match = tuple(sorted([a - row + 1 for a in working_rc[row - 1]]))
        lower_perm = None
        # logger.debug(to_match)
        for L, new_perm in pull_out_var(row, working_rc.perm):
            # logger.debug("Trying   ", L)

            if tuple(sorted(L)) == to_match:
                lower_perm = new_perm
                # logger.debug("Good")
                break
        assert lower_perm is not None
        build_rc = RCGraph(working_rc[:-1])
        # SPOTS EQUAL
        # just those refs
        # just strings
        # reflections at the spots
        for r in range(row - 1, 0, -1):
            to_match = tuple(sorted([a - r + 1 for a in working_rc[r - 1]]))
            lower_perm2 = None
            # logger.debug(f"{lower_perm=}")
            # logger.debug(f"row {r=}")
            # logger.debug(to_match)
            for L, new_perm in pull_out_var(r, lower_perm):
                # logger.debug("Trying   ", L)
                if tuple(sorted(L)) == to_match:
                    lower_perm2 = new_perm
                    # logger.debug("Good")
                    break
            if lower_perm2 is None:
                # logger.debug("Failed on row ", r)
                break
            else:
                lower_perm = lower_perm2
                # logger.debug("Succeeded on row ", r)

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

    def left_to_right_inversion(self, index):
        return self.right_root_at(*self.left_to_right_inversion_coord(index))

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

    def max_of_row(self, r):
        if len(self[r - 1]) == 0:
            the_max = 1
            while self.right_root_at(r, the_max + 1) != (r + the_max, r + the_max + 1):
                the_max += 1
        else:
            the_max = max(self[r - 1]) - r + 1
        return the_max

    @staticmethod
    def find_reflection_path(bottom_perm, top_perm, row, last_spot=0, last_value=-1000):
        from schubmult.utils.perm_utils import has_bruhat_ascent

        working_perm = bottom_perm
        # logger.debug(working_perm)
        # logger.debug(f"{top_perm=}")
        i = row - 1
        for j in range(row, len(top_perm)):
            if has_bruhat_ascent(working_perm, i, j) and working_perm[j] > last_value:
                new_perm = working_perm.swap(i, j)
                if new_perm == top_perm:
                    reflection_path = {}
                    reflection_path[i + 1] = reflection_path.get(i + 1, [])
                    reflection_path[i + 1].append(j + 1)
                    return reflection_path
                if new_perm.bruhat_leq(top_perm):
                    reflection_path = RCGraph.find_reflection_path(new_perm, top_perm, row, last_spot=i, last_value=working_perm[j])
                    if reflection_path is not None:
                        reflection_path = dict(reflection_path)
                        reflection_path[i + 1] = reflection_path.get(i + 1, [])
                        reflection_path[i + 1].insert(0, j + 1)
                        last_value = working_perm[i]
                        working_perm = new_perm
                        return reflection_path
        return None

    @classmethod
    def schub_hom(cls, perm, length):
        word_elem = ASx(perm, length).change_basis(WordBasis)
        result = 0
        for word, coeff in word_elem.items():
            result += coeff * RCGraph.fa_hom(*word)
        return result

    @classmethod
    def schub_coproduct(cls, perm, length):
        word_elem = ASx(perm, length).change_basis(WordBasis)
        result = 0
        for word, coeff in word_elem.items():
            cop = coeff * FA(*word).coproduct()
            for (word1, word2), cf2 in cop.items():
                result += cf2 * RCGraph.fa_hom(*word1) @ RCGraph.fa_hom(*word2)
        return result

    # @cache
    # def coproduct(self):
    #     def trans_to_h_list(arr):
    #         h_list = []
    #         buildup_code = []
    #         cnt = 0
    #         arr = list(reversed(arr))
    #         # logger.debug(f"{arr=}")
    #         first_saw = -1
    #         for i in range(len(arr)):
    #             cnt += 1
    #             if first_saw == -1:
    #                 first_saw = arr[i]
    #             if i > 0 and arr[i - 1] != arr[i] - 1:
    #                 buildup_code = (first_saw - 1) * [0]
    #                 buildup_code.append(cnt)
    #                 h_list.append(uncode(buildup_code))
    #                 first_saw = arr[i]
    #                 # logger.debug(f"{buildup_code=}")
    #                 buildup_code = []
    #                 cnt = 0

    #         buildup_code = (first_saw - 1) * [0]
    #         buildup_code.append(cnt)
    #         h_list.append(uncode(buildup_code))
    #         lst = [cd.trimcode for cd in h_list]
    #         # logger.debug(f"{lst=}")
    #         # logger.debug(f"{arr=}")
    #         assert sum(cd.inv for cd in h_list) == len(arr), f"{h_list=}, {arr=}"
    #         return h_list

    #     # commuting h's
    #     if len(self) == 0:
    #         return RCGraphModule({RCGraph(): 1}) @ RCGraphModule({RCGraph(): 1})
    #     if len(self) == 1 and len(self[0]) == 0:
    #         return RCGraphModule({RCGraph([()]): 1}) @ RCGraphModule({RCGraph([()]): 1})
    #     if self.perm.inv == 0:
    #         return RCGraphModule({RCGraph([()] * len(self)): 1}) @ RCGraphModule({RCGraph([()] * len(self)): 1})
    #     buildup_module = RCGraphModule({RCGraph([]): 1}) @ RCGraphModule({RCGraph([]): 1})

    #     for row in range(len(self) - 1, -1, -1):
    #         ret_elem = 0
    #         # buildup_module = RCGraph.pad_tensor_module_with_zeros(buildup_module, 1)
    #         the_row = self.rowrange(row, row + 1)
    #         h_list = trans_to_h_list(the_row[0])
    #         # if len(self) == 1:

    #         # logger.debug("Buildup is")
    #         # logger.debug(buildup_module)

    #         for (rc1, rc2), coeff in buildup_module.items():
    #             # logger.debug("Multiplying")
    #             # logger.debug(rc1)
    #             # logger.debug("and")
    #             # logger.debug(rc2)

    #             # logger.debug("Product is")
    #             # logger.debug(prod_module)
    #             rc01 = RCGraph([(), *rc1.shiftup(1)])
    #             rc02 = RCGraph([(), *rc2.shiftup(1)])

    #             # logger.debug(f"{rc12=}")
    #             # logger.debug(f"{rc11}")
    #             for perm in h_list:
    #                 perm_set = ASx(perm).coproduct()
    #                 for ((p1, _), (p2, _)), _ in perm_set.items():
    #                     rc011 = rc01
    #                     rc012 = rc02
    #                     if p1.inv > 0:
    #                         rc011 = rc011.kogan_kumar_insert(len(p1.trimcode), [1] * p1.inv)
    #                         rc011 = rc011.rowrange(0, 1) * rc1
    #                         new_rc011 = 0
    #                         for rc, c in rc011.items():
    #                             new_rc011 += c * (RCGraph([(), *rc[1:]]).kogan_kumar_insert(len(p1.trimcode), [1] * p1.inv))
    #                         rc011 = new_rc011
    #                     if p2.inv > 0:
    #                         rc012 = rc012.kogan_kumar_insert(len(p2.trimcode), [1] * p2.inv)
    #                         rc012 = rc012.rowrange(0, 1) * rc2
    #                         new_rc012 = 0
    #                         for rc, c in rc012.items():
    #                             new_rc012 += c * (RCGraph([(), *rc[1:]]).kogan_kumar_insert(len(p2.trimcode), [1] * p2.inv))
    #                         rc012 = new_rc012
    #                     rc011 = 1 * rc011
    #                     rc012 = 1 * rc012
    #                     # rc012 = RCGraph([(),*rc012[1:]]).kogan_kumar_insert(len(p2.trimcode),[1]*p2.inv)
    #                     # if rc011.perm.bruhat_leq(self.rowrange(0,row+1).perm) and rc012.perm.bruhat_leq(self.rowrange(0,row+1).perm):
    #                     # logger.debug("Matched")
    #                     # logger.debug(rc011)
    #                     # logger.debug(rc012)
    #                     # logger.debug("With coeff", coeff)
    #                     # logger.debug("and perm", perm)
    #                     ret_elem += coeff * (rc011 @ rc012)
    #         buildup_module = ret_elem
    #         # logger.debug("mul_module")
    #         # logger.debug(mul_module)
    #         # logger.debug("ret_elem_first")
    #         # logger.debug(ret_elem)
    #         # if len(self) > 1:
    #         #     mul_module = RCGraph(self.rowrange(1, len(self))).coproduct()
    #         #     ret_elem = ret_elem*mul_module

    #     # logger.debug("Final result for")
    #     # logger.debug(self)
    #     # logger.debug("is")
    #     # logger.debug(ret_elem)
    #     # ret_elem *= self.rowrange(1, len(self)).coproduct()
    #     # logger.debug("Result is")
    #     # logger.debug(ret_elem)

    #     return ret_elem

    @classmethod
    def principal_rc(cls, perm, length):
        cd = perm.trimcode
        graph = []
        for i in range(len(cd)):
            row = tuple(range(i + cd[i], i, -1))
            graph.append(row)
        graph = [*graph, *[()] * (length - len(graph))]
        return cls(graph)

    # def pad_with_zeros(self, num_zeros):
    #     return RCGraph.pad_module_with_zeros(1 * self, num_zeros)

    # @staticmethod
    # def pad_module_with_zeros(rc_graph_module, num_zeros):
    #     if num_zeros == 0:
    #         return rc_graph_module
    #     build_module = rc_graph_module
    #     for _ in range(num_zeros):
    #         new_build_module = 0
    #         for rc, coeff in build_module.items():
    #             new_rc = RCGraphModule(dict.fromkeys(rc.right_zero_act(), coeff))
    #             new_build_module += new_rc
    #         build_module = new_build_module
    #     return build_module

    # @staticmethod
    # def pad_tensor_module_with_zeros(rc_graph_module, num_zeros):
    #     build_module = rc_graph_module
    #     for _ in range(num_zeros):
    #         new_build_module = 0
    #         for (rc1, rc2), coeff in build_module.items():
    #             new_rc = RCGraphModule(dict.fromkeys(rc1.right_zero_act(), coeff)) @ RCGraphModule(dict.fromkeys(rc2.right_zero_act(), coeff))
    #             new_build_module += new_rc
    #         build_module = new_build_module
    #     return build_module

    # def prod_with_rc(self, other):
    #     # logger.debug("Mulling")
    #     # logger.debug(self)
    #     # logger.debug("and")
    #     # logger.debug(other)
    #     # if len(other) == 0:
    #     #     return 1 * self
    #     # orig_len = len(other)

    #     # base_rc_set = {self}
    #     # for _ in range(dff):
    #     #     new_base_rc_set = set()
    #     #     for rc in base_rc_set:
    #     #         new_base_rc_set.update(rc.right_zero_act())
    #     #     base_rc_set = new_base_rc_set
    #     # #base_rc = RCGraph([*self, *[()]*max(dff, 0)])
    #     dff2 = len(other) if other.perm.inv == 0 else max(max(other.perm.descents()) + 1 - len(other), 0)
    #     # base_other_rc_set = {other}

    #     # ret_module = RCGraphModule()
    #     # for base_rc in base_rc_set:
    #     #     for base_other_rc in base_other_rc_set:
    #     dff = 0 if self.perm.inv == 0 else max((max(self.perm.descents()) + 1)-len(self),0)
    #     # logger.debug(dff)
    #     # logger.debug(self)
    #     upd = self
    #     if dff > 0:
    #         upd = RCGraph([*self, *[()]*dff])

    #     #         #buildup_module = 1*base_rc
    #     #         buildup_module = 1 * base_rc
    #     #         num_zeros = min(orig_len - dff, 0)
    #     #         for _ in range(num_zeros):
    #     #             new_buildup_module = 0
    #     #             for rc, coeff in buildup_module.items():
    #     #                 new_buildup_module += RCGraphModule(dict.fromkeys(base_rc.right_zero_act(), coeff))
    #     #             buildup_module = new_buildup_module

    #     #         new_buildup_module = 0
    #     interim_ret_module = 0
    #     new_rc = upd.pad_with_zeros(len(other)-dff)
    #     for update_rc, cf in new_rc.items():
    #         new_new_rc = RCGraph([*update_rc[:len(self)], *other.shiftup(len(self))])
    #         if new_new_rc.is_valid:
    #             interim_ret_module += cf * new_new_rc
    #     ret_module = interim_ret_module
    #     # ret_module = 0
    #     # for rc, coeff in interim_ret_module.items():
    #     #     ret_module += coeff * rc.rowrange(0,len(self)+len(other))
    #     #             # logger.debug("new_rc")
    #     #             # logger.debug(new_rc)
    #     #             # logger.debug("Trying to match")
    #     #             # logger.debug(new_rc)
    #     #             num_zeros = dff2
    #     #             new_base_other_rc_set = {base_other_rc_set}
    #     #             for _ in range(dff2):
    #     #                 new_new_base_other_rc_set = set()
    #     #                 for rc in base_other_rc_set:
    #     #                     new_base_other_rc_set.update(rc.right_zero_act())
    #     #                 new_base_other_rc_set = new_base_other_rc_set
    #     #             for rc, coeff in buildup_module.items():
    #     #                 new_buildup_module += RCGraphModule(dict.fromkeys(base_rc.right_zero_act(), coeff))
    #     #             buildup_module = new_buildup_module
    #     #             if new_rc.is_valid:
    #     #                 # logger.debug("Matched")
    #     #                 new_buildup_module += coeff * new_rc
    #     #         # for rc2, coeff2 in new_buildup_module.items():
    #     #         #     ret_module += coeff2 * rc2.normalize_and_remove_last_row()
    #     return ret_module

    def zero_rectify(self):
        from schubmult.utils.perm_utils import has_bruhat_ascent

        if self.is_valid:
            return self
        working_rc = self
        while not working_rc.is_valid:
            mid_rc = RCGraph()
            index = 1
            while index <= len(self) and mid_rc.is_valid:
                mid_rc = working_rc.rowrange(0, index)
                index += 1
            # logger.debug("Got")
            # logger.debug(mid_rc)
            assert not mid_rc.is_valid
            for index2 in range(len(mid_rc.perm_word()) - 1, -1, -1):
                if mid_rc.is_valid:
                    break
                a, b = mid_rc.left_to_right_inversion(index2)
                if a > b and mid_rc.perm[a - 1] > mid_rc.perm[b - 1]:
                    row, col = mid_rc.left_to_right_inversion_coord(index2)
                    assert col > 1
                    mid_rc = mid_rc.toggle_ref_at(row, col)
                    for col2 in range(col, 0, -1):
                        a2, b2 = mid_rc.right_root_at(row, col2)
                        # logger.debug(f"{a2, b2}")
                        if a2 < b2 and has_bruhat_ascent(mid_rc.perm, a2 - 1, b2 - 1) and not mid_rc.has_element(row, col2):
                            # logger.debug("got it")
                            mid_rc = mid_rc.toggle_ref_at(row, col2)
                            break

            # logger.debug("After")
            # logger.debug(mid_rc)
            working_rc = RCGraph([*mid_rc, *working_rc[len(mid_rc) :]])
        return working_rc

    # THE ZERO MAKES SCHUB PROD
    def prod_with_rc(self, other):
        # if len(other) == 0:
        #     return 1 * self
        if self.perm.inv == 0:
            return {RCGraph([*self, *other.shiftup(len(self))]): 1}
        num_zeros = len(other) if other.perm.inv == 0 else max(len(other), len(other.perm))
        base_rc = self
        buildup_module = {base_rc: 1}
        # num_zeros += len(self) - len(base_rc)
        for _ in range(num_zeros):
            new_buildup_module = {}
            for rc, coeff in buildup_module.items():
                new_buildup_module = add_perm_dict(new_buildup_module, dict.fromkeys(rc.right_zero_act(), coeff))
            buildup_module = new_buildup_module

        # logger.debug("Buildup is")
        # logger.debug(buildup_module)
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

    # def as_str_lines(self):
    #     lines = []
    #     for i in range(len(self)):
    #         row = self[i]
    #         splug_row = [*row, i]
    #         row = [str(splug_row[j]) + ("  ") * (splug_row[j] - splug_row[j + 1] - 1) for j in range(len(splug_row) - 1)]
    #         line = ""
    #         line += " ".join([str(r) for r in row])
    #         lines += [line]
    #     if len(lines) == 0:
    #         lines = [""]
    #     lines2 = []
    #     ml = max([len(line) for line in lines])
    #     if len(lines[0]) < ml:
    #         lines[0] = " " * (ml - len(lines[0])) + lines[0]
    #     for line in lines:
    #         lines2 += [" " * (ml - len(line)) + line]
    #     return lines2

    def __ge__(self, other):
        return not (self < other)

    def __gt__(self, other):
        return not (self <= other)

    # def __leq__(self, other):
    #     if not isinstance(other, RCGraph):
    #         return NotImplemented
    #     if len(self) != len(other):
    #         return NotImplemented
    #     for i in range(len(self)):
    #         perm1 = Permutation.ref_product(*self[i])
    #         perm2 = Permutation.ref_product(*other[i])
    #         if not perm1.bruhat_leq(perm2):
    #             return False
    #     return True

    @property
    def is_principal(self):
        return self.perm == uncode(self.length_vector())

    # def __str__(self):
    #     lines = self.as_str_lines()
    #     lines2 = [line for line in lines]
    #     return "\n".join(lines2)

    @property
    def inv(self):
        return self.perm.inv

    def _sympystr(self, printer=None):
        from sympy.printing.str import StrPrinter

        if not printer:
            printer = StrPrinter()
        return printer._print_MatrixBase(self)

    @property
    def rows(self):
        return len(self)

    @cached_property
    def cols(self):
        return max(max((row[i] + i for i in range(len(row))) if len(row) > 0 else [0]) for row in self) if len(self) > 0 else 0


    def __getitem__(self, key):
        # if len(args) == 1:
        #     logger.debug(repr(args),len(args))
        #     return tuple(self)[args[0]]
        if isinstance(key, int):
            return tuple(self)[key]
        if isinstance(key, tuple):
            i, j = key
            if not self.has_element(i + 1, j + 1 ):
                return None
            return i + j + 1
        is_slice = isinstance(key, slice)

        if is_slice:
            values = tuple(tuple(self)[n] for n in range(len(self))[key])
            return values
        raise ValueError(f"Bad indexing {key=}")

    def _format_str(self, printer=None) -> str:
        from sympy.printing import StrPrinter

        if not printer:
            printer = StrPrinter()
        # Handle zero dimensions:
        if self.rows == 1:
            return "RCGraph([%s])" % self.table(printer, rowsep=",\n")
        return "RCGraph([\n%s])" % self.table(printer, rowsep=",\n")

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

    # def __repr__(self):
    #     return f"{self.__class__.__name__}(" + ", ".join([repr(k) for k in self]) + ")"

    def __hash__(self):
        return hash(tuple(self))

    def __lt__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        # return self.perm.bruhat_leq(other.perm) and self.perm != other.perm
        if self.perm.trimcode < other.perm.trimcode:
            return True
        if self.perm.trimcode > other.perm.trimcode:
            return False
        for i in range(self.perm.inv):
            a, b = self.perm.right_root_at(i)
            # if self.inversion_label(a - 1, b - 1) < other.inversion_label(a - 1, b - 1):
            #     return True
            # if self.inversion_label(a - 1, b - 1) > other.inversion_label(a - 1, b - 1):
            #     return False

            if self.inversion_label(a - 1, b - 1) < other.inversion_label(a - 1, b - 1):
                return True
            if self.inversion_label(a - 1, b - 1) > other.inversion_label(a - 1, b - 1):
                return False
        return False

    def __le__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return self < other or self == other
