from functools import cache, cached_property

import schubmult.utils.schub_lib as schub_lib
from schubmult.rings.free_algebra import ASx, FreeAlgebra, FreeAlgebraElement, WordBasis
from schubmult.rings.nil_hecke import NilHeckeRing
from schubmult.schub_lib.perm_lib import Permutation, uncode
from schubmult.symbolic import S, prod
from schubmult.utils._grid_print import GridPrint

# from schubmult.utils.bitfield_row import BitfieldRow
from schubmult.utils.logging import get_logger, init_logging
from schubmult.utils.perm_utils import add_perm_dict

from .crystal_graph import CrystalGraph, CrystalGraphTensor
from .nilplactic import NilPlactic
from .plactic import Plactic

init_logging(debug=False)
logger = get_logger(__name__)


def _is_row_root(row, root):
    return root[0] <= row and root[1] > row


FA = FreeAlgebra(WordBasis)


def debug_print(*args, debug=False):
    if debug:
        print(*args)


def _crystal_isomorphic(c1, c2, cutoff=None):
    hw_1, _ = c1.to_highest_weight(length=cutoff)
    hw_2, _ = c2.to_highest_weight(length=cutoff)

    stack = [(c1, c2)]
    if cutoff is None:
        cutoff = c1.crystal_length()
    if hw_1.crystal_weight != hw_2.crystal_weight:
        return False
    while len(stack) > 0:
        c1_test, c2_test = stack.pop()
        for i in range(1, cutoff):
            c1_test0 = c1_test.lowering_operator(i)
            if c1_test0 is not None:
                c2_test0 = c2_test.lowering_operator(i)
                if c2_test0 is None:
                    return False
                stack.append((c1_test0, c2_test0))
    return True


class RCGraph(GridPrint, tuple, CrystalGraph):
    def div_diff(self, i):
        if i >= self.crystal_length():
            return None
        tst = self.lowering_operator(i)
        if tst is None:
            return self.exchange_property(i)
        return None

    def __eq__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return tuple(self) == tuple(other)

    # def __iter__(self):
    #     for row in super().__iter__():
    #         yield tuple(row)

    def flat_elem_sym_mul(self, k):
        from schubmult.utils.schub_lib import elem_sym_perms

        elem_graph = RCGraph([(i,) for i in range(1, k + 1)])
        mul_graph = self
        if len(elem_graph) != len(mul_graph):
            length = max(len(elem_graph), len(mul_graph))
            elem_graph = elem_graph.resize(length)
            mul_graph = mul_graph.resize(length)
        tensor = CrystalGraphTensor(elem_graph, mul_graph)
        hw_rc, raise_seq = tensor.to_highest_weight()
        perm_list = [perm for (perm, w) in elem_sym_perms(mul_graph.perm, k, k) if w == k]
        results = set()
        for perm in perm_list:
            for rc in RCGraph.all_rc_graphs(perm, length=len(mul_graph), weight=hw_rc.crystal_weight):
                if rc.is_highest_weight and _crystal_isomorphic(hw_rc, rc):
                    results.add(rc.reverse_raise_seq(raise_seq))
                    break
        if results:
            assert len(results) == 1, (
                f"Ambiguous flat elem crystal multiplication results for k={k} on\n{self} \nResults:\n" + "\n".join([str(r) for r in results]) + "\n".join([str(r.p_tableau) for r in results])
            )
            return next(iter(results))
        return None

    def monk_crystal_mul(self, p, k, warn=True):
        if p > k:
            raise ValueError("p must be less than or equal to k")
        if k > len(self):
            return self.extend(k - len(self)).monk_crystal_mul(p, k)

        def _crystal_isomorphic(c1, c2, cutoff=None):
            hw_1, _ = c1.to_highest_weight(length=cutoff)
            hw_2, _ = c2.to_highest_weight(length=cutoff)

            stack = [(c1, c2)]
            if cutoff is None:
                cutoff = c1.crystal_length()
            if hw_1.crystal_weight != hw_2.crystal_weight:
                return False
            while len(stack) > 0:
                c1_test, c2_test = stack.pop()
                for i in range(1, cutoff):
                    c1_test0 = c1_test.lowering_operator(i)
                    if c1_test0 is not None:
                        c2_test0 = c2_test.lowering_operator(i)
                        if c2_test0 is None:
                            return False
                        stack.append((c1_test0, c2_test0))
            return True

        from schubmult.utils.schub_lib import elem_sym_perms

        if k > len(self):
            return self.extend(k - len(self)).monk_crystal_mul(p, k)

        monk_rc = next(iter(RCGraph.all_rc_graphs(Permutation([]).swap(k - 1, k), len(self), weight=(*([0] * (p - 1)), 1, *([0] * (len(self) - p))))))

        results = set()
        lv = [*self.length_vector]
        lv[p - 1] += 1
        up_perms = [pperm for pperm, L in elem_sym_perms(self.perm, 1, k) if L == 1]
        for up_perm in up_perms:
            for rc2 in RCGraph.all_rc_graphs(up_perm, length=len(self), weight=lv):
                try:
                    good = True
                    for start in range(p):
                        pp = p - start
                        for cut in range(pp, k + 1 - start):
                            rc2_hw, _ = rc2.rowrange(start).vertical_cut(cut)[0].to_highest_weight()
                            tensor_cut, raise_seq = CrystalGraphTensor(self.rowrange(start).vertical_cut(cut)[0], monk_rc.rowrange(start).vertical_cut(cut)[0]).to_highest_weight()
                            if rc2_hw.crystal_weight != tensor_cut.crystal_weight or rc2_hw.reverse_raise_seq(raise_seq) != rc2.rowrange(start).vertical_cut(cut)[0]:
                                good = False
                                break
                            rc2_lw, _ = rc2.rowrange(start).vertical_cut(cut)[0].to_lowest_weight()
                            tensor_cut_low, lower_seq = CrystalGraphTensor(self.rowrange(start).vertical_cut(cut)[0], monk_rc.rowrange(start).vertical_cut(cut)[0]).to_lowest_weight()
                            if rc2_lw.crystal_weight != tensor_cut_low.crystal_weight or rc2_lw.reverse_lower_seq(lower_seq) != rc2.rowrange(start).vertical_cut(cut)[0]:
                                good = False
                                break
                    if rc2[k:] != self[k:]:
                        good = False
                    # tensor, raise_seq = CrystalGraphTensor(self, monk_rc).to_highest_weight()
                    # hw_rc, _ = rc2.to_highest_weight()
                    # if hw_rc.crystal_weight != tensor.crystal_weight:# or hw_rc.reverse_raise_seq(raise_seq) != rc2:
                    #     good = False
                    if good:
                        results.add(rc2)
                except Exception as e:  # noqa: F841
                    continue
        try:
            assert len(results) == 1, f"Ambiguous monk crystal multiplication results for p={p}, k={k} on\n{self} \nResults:\n" + "\n".join([str(r) for r in results])
        except AssertionError as e:
            print(e)
            if not warn:
                raise
        return next(iter(results))

    @cached_property
    def crystal_weight(self):
        return self.length_vector

    # UNIQUE
    def tableau_decomp(self):
        descs = self.perm.descents()
        if len(descs) == 0:
            return (self,)
        dscs = sorted([d + 1 for d in descs], reverse=True)

        tup = (self,)
        for d in dscs:
            tup = (*tup[0].vertical_cut(d), *tup[1:])

        return tup[:-1]


    @property
    def is_rc(self):
        for i, row in enumerate(self):
            for a in row:
                if a < i + 1:
                    return False
        return True

    @property
    def is_valid(self):
        if self.perm.inv != len(self.perm_word):
            return False
        
        # if len(self.perm.trimcode) > len(self):
        #     return False
        return True

    def shiftup(self, shift=1):
        rc = self
        # if len(self) < len(self.perm.trimcode) + shift:
        #     rc = rc.extend(len(self.perm.trimcode) + shift - len(self))

        ret = RCGraph([tuple([a + shift for a in rrow]) for rrow in rc])
        assert ret.is_valid
        return ret

    @cache
    def right_root_at(self, i, j):
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

    def reverse_kogan_kumar_insert(self, descent, reflection_path, return_rows=False):
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
            if root[0] not in pair_dict:
                if root[0] in pair_dict_rev and root[1] in pair_dict_rev and pair_dict_rev[root[0]] == pair_dict_rev[root[1]]:
                    return True
                return False
            return root[0] in pair_dict and root[1] in pair_dict[root[0]]

        # may have to add q, s or a_i, q
        def is_relevant_noncrossing(root):
            top, bottom = max(root[1], root[0]), min(root[0], root[1])
            return (bottom <= descent and descent < top and top not in pair_dict_rev) or (bottom in pair_dict_rev and top > descent and top not in pair_dict_rev)

        # Add this intersection. If we are in the first case, insert (s, q) into the sequence (ai, bi) in the rightmost position, such that ai’s remain nondecreasing in the
        # sequence. ((s, q) are the rows where the two strands shown in Figure 3 originate.) If
        # we are in the second case, add (ai, q) just before where (a, bi) is in the sequence.

        working_rc = self

        rows = []

        for row in range(1, len(self) + 1):
            for col in range(1, max(descent, self.cols) + 1):
                if working_rc.has_element(row, col):
                    a, b = working_rc.right_root_at(row, col)
                    if is_relevant_crossing((a, b), working_rc.perm):
                        working_rc = working_rc.toggle_ref_at(row, col)
                        a2 = a
                        if a2 in pair_dict_rev:
                            a2 = pair_dict_rev[a2]

                        pair_dict[a2].remove(b)
                        del pair_dict_rev[b]
                        if len(pair_dict[a2]) == 0:
                            del pair_dict[a2]

                        rows.append(row)

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
        for index in range(len(self.perm_word)):
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

    @property
    def reduced_word(self):
        return self.perm_word

    def is_dom_perm_yamanouchi(self, dom_perm, perm):
        from schubmult.rings.schubert_ring import Sx

        if (Sx(self.perm) * Sx(dom_perm)).get(perm, 0) == 0:
            return False
        length = max(len(perm.trimcode), len(dom_perm.trimcode))
        rc = RCGraph.principal_rc(perm, length)
        dom_rc = RCGraph.principal_rc(dom_perm, length)
        weight = tuple([rc.length_vector[i] - dom_rc.length_vector[i] for i in range(len(rc))])
        if self.length_vector != weight:
            return False
        outer_shape = rc.p_tableau.shape
        inner_shape = dom_rc.p_tableau.shape
        if NilPlactic.exists_ed_tableau_equiv(rc.p_tableau, inner_shape, outer_shape) and Plactic.exists_ss_tableau_equiv(rc.weight_tableau, inner_shape, outer_shape):
            return True  # should match highest weight of tensor
        return False

    @property
    def shape(self):
        P = self.edelman_greene()[0]
        shape = tuple(len(P[i]) for i in range(len(P)))
        return shape

    # product = demazure crystal to demazure crystal
    # transpose is weight preserving

    def __invert__(self):
        new_rc = RCGraph([()] * self.cols)
        for i in range(1, self.rows + 1):
            for j in range(1, self.cols + 1):
                if self.has_element(i, j):
                    new_rc = new_rc.toggle_ref_at(j, i)
        return new_rc

    def normalize(self):
        return self.resize(len(self.perm.trimcode))

    def resize(self, new_length):
        if new_length < len(self):
            return self.rowrange(0, new_length)
        return self.extend(new_length - len(self))

    def edelman_greene(self):
        word1, word2 = (), ()
        index = 0
        for index, (row, col) in enumerate(list(reversed([self.left_to_right_inversion_coords(i) for i in range(self.perm.inv)]))):
            to_insert = row + col - 1
            word1, word2 = NilPlactic.ed_insert_rsk(word1, word2, to_insert, len(self) - row + 1)
            index += 1
        P = word1
        Q = word2

        # reg._rc_graph = self
        return (NilPlactic(P), Plactic(Q))

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
        return i <= len(self) and i + j  - 1 in self[i - 1]

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

    def rowrange(self, start, end=None):
        if not end:
            end = len(self)
        if start == end:
            return type(self)(())
        return type(self)([tuple([a - start for a in row]) for row in self[start:end]])

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
            cls._graph_cache[(perm, length)] = {cls([()] * length if length > 0 else [])}
            return cls._graph_cache[(perm, length)]
        ret = set()
        pm = perm
        L = schub_lib.pull_out_var(1, pm)
        for _, new_perm in L:
            new_row = [new_perm[i] for i in range(max(len(pm), len(new_perm))) if new_perm[i] == pm[i + 1]]
            if weight and len(new_row) != weight[0]:
                continue
            new_row.sort(reverse=True)
            if weight is not None:
                oldset = cls.all_rc_graphs(new_perm, length=length - 1, weight=weight[1:])
            else:
                oldset = cls.all_rc_graphs(new_perm, length=length - 1)
            for old_rc in oldset:
                nrc = cls([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in old_rc]])
                assert nrc.perm == perm
                assert len(nrc) == length
                if weight:
                    assert nrc.length_vector == tuple(weight)
                ret.add(nrc)
        if weight:
            cls._cache_by_weight[(perm, tuple(weight))] = ret
        else:
            cls._graph_cache[(perm, length)] = ret
        return ret

    def extend(self, extra_rows):
        return type(self)([*self, *tuple([()] * extra_rows)])

    def _kogan_kumar_insert_row(self, row, descent, dict_by_a, dict_by_b, num_times, start_index=-1, backwards=True):
        working_rc = self
        if row > descent:
            raise ValueError("All rows must be less than or equal to descent")

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
                i = working_rc.cols + descent + 5
            if not backwards:
                i -= 1
            else:
                i += 1
            flag = False

            if not working_rc.has_element(row, i):
                a, b = working_rc.right_root_at(row, i)
                if a < b:
                    flag = False
                    if _is_row_root(descent, (a, b)) and b not in dict_by_b:
                        working_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[a] = dict_by_a.get(a, set())
                        dict_by_a[a].add(b)
                        dict_by_b[b] = a
                        flag = True

                    elif a in dict_by_b and b > descent and b not in dict_by_b:
                        working_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[dict_by_b[a]].add(b)
                        dict_by_b[b] = dict_by_b[a]
                        flag = True

                if flag:
                    num_done += 1
                if row > 1 and not working_rc.is_valid:
                    working_rc = working_rc._kogan_kumar_rectify(row - 1, descent, dict_by_a, dict_by_b, backwards=backwards)  # minus one?
        return working_rc

    def _kogan_kumar_rectify(self, row_below, descent, dict_by_a, dict_by_b, backwards=True):
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

            if working_rc.has_element(row_below, j):
                a, b = working_rc.right_root_at(row_below, j)

                top, bottom = max(a, b), min(a, b)

                if a < b:
                    continue

                if bottom in dict_by_a and top in dict_by_a[bottom]:
                    new_rc = working_rc.toggle_ref_at(row_below, j)
                    dict_by_a[bottom].remove(top)
                    if len(dict_by_a[bottom]) == 0:
                        del dict_by_a[bottom]
                    del dict_by_b[top]
                    working_rc = new_rc
                    flag = True

                elif bottom in dict_by_b and top in dict_by_b and dict_by_b[top] == dict_by_b[bottom]:
                    new_rc = working_rc.toggle_ref_at(row_below, j)
                    dict_by_a[dict_by_b[bottom]].remove(top)
                    del dict_by_b[top]
                    if len(dict_by_a[dict_by_b[top]]) == 0:
                        del dict_by_a[dict_by_b[top]]
                    flag = True
                    working_rc = new_rc

                else:
                    raise ValueError(f"Could not rectify at {(row_below, j)} with root {(a, b)}")
                if flag:
                    working_rc = working_rc._kogan_kumar_insert_row(row_below, descent, dict_by_a, dict_by_b, num_times=1, backwards=backwards)
        return working_rc._kogan_kumar_rectify(row_below - 1, descent, dict_by_a, dict_by_b, backwards=backwards)

    # VERIFY
    def kogan_kumar_insert(self, descent, rows, return_reflections=False, backwards=True):
        dict_by_a = {}
        dict_by_b = {}
        # row is descent
        # inserting times

        working_rc = type(self)([*self])
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
        for row in sorted(rows_grouping.keys(), reverse=True):
            num_times = rows_grouping[row]
            last_working_rc = working_rc
            working_rc = working_rc._kogan_kumar_insert_row(row, descent, dict_by_a, dict_by_b, num_times, backwards=backwards)

            if row > 1 and not working_rc.is_valid:
                working_rc = working_rc._kogan_kumar_rectify(row - 1, descent, dict_by_a, dict_by_b, backwards=backwards)  # minus one?

            try:
                assert len(working_rc[row - 1]) == len(last_working_rc[row - 1]) + num_times
            except AssertionError:
                raise
        if return_reflections:
            reflections = []
            for a in sorted(dict_by_a.keys()):
                reflections = [*reflections, *[(a, b) for b in dict_by_a[a]]]
            return working_rc, tuple(reflections)
        return working_rc

    @property
    def weight(self):
        wt = []
        for i, row in enumerate(self):
            wt.extend([i + 1] * len(row))
        return tuple(wt)

    @property
    def perm(self):
        perm = Permutation([])
        for row in self:
            for p in row:
                perm = perm.swap(p - 1, p)
        return perm

    def transpose(self, length=None):
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
        new_rc = (type(self)(newrc)).normalize()

        assert new_rc.perm == ~self.perm
        if length is not None:
            new_rc = new_rc.resize(length)
        return new_rc

    @classmethod
    def one_row(cls, p):
        return cls((tuple(range(p, 0, -1)),))

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

    def toggle_ref_at(self, i, j):
        if i <= 0 or j <= 0:
            raise IndexError()
        new_row = [*self[i - 1]]
        if i + j - 1 in new_row:
            index = new_row.index(i + j - 1)
            new_row = [*new_row[:index], *new_row[index + 1 :]]
        else:
            index = 0
            if len(new_row) > 0:
                while index < len(new_row) and new_row[index] > i + j - 1:
                    index += 1
            new_row.insert(index, i + j - 1)
        return type(self)([*self[: i - 1], tuple(new_row), *self[i:]])

    # # THIS IS KEY
    # # EXCHANGE PROPERTY GOES TO UNIQUE PERMUTATION
    # # KOGAN INSERT ENSURES WE GO UP PROPERLY

    _z_cache = {}  # noqa: RUF012

    @cache
    def zero_out_last_row(self):
        # this is important!
        # transition formula
        if len(self[-1]) != 0:
            raise ValueError("Last row not empty")
        if self.perm.inv == 0:
            return self.rowrange(0, len(self) - 1)
        interim = type(self)([*self])

        diff_rows = []
        descs = []
        extend_amount = 1

        while len(interim.perm.trimcode) > len(self) - 1:
            descs += [len(interim.perm.trimcode)]
            interim, row = interim.exchange_property(len(interim.perm.trimcode), return_row=True)
            diff_rows += [row]

        interim2 = type(self)([*interim[:-1], tuple(sorted(descs, reverse=True))])
        interim = interim2.kogan_kumar_insert(len(self.perm.trimcode) - extend_amount, diff_rows)

        return interim.rowrange(0, len(self) - 1)

    def crystal_length(self):
        return len(self)

    def lowering_operator(self, row):
        # RF word is just the RC word backwards
        if row >= len(self):
            return None
        row_i = [*self[row - 1]]
        row_ip1 = [*self[row]]

        # pair the letters
        pairings = []
        unpaired = []
        unpaired_b = [*row_ip1]

        for letter in row_i:
            st = [letter2 for letter2 in unpaired_b if letter2 > letter]
            if len(st) == 0:
                unpaired.append(letter)
            else:
                pairings.append((letter, min(st)))
                unpaired_b.remove(min(st))
        if len(unpaired) == 0:
            return None
        b = min(unpaired)
        t = min([j for j in range(b) if b - j - 1 not in row_i])

        if b - t < row + 1:
            return None
        new_row_i = [s for s in row_i if s != b]
        new_row_ip1 = sorted([b - t, *row_ip1], reverse=True)
        ret_rc = type(self)([*self[: row - 1], tuple(new_row_i), tuple(new_row_ip1), *self[row + 1 :]])
        if ret_rc.perm != self.perm:
            return None
        return ret_rc

    # weak edelman green correspondence
    # better understanding of crystal
    # Following [Assa], for P a semi-standard Young tableau with strictly increasing rows, define the lift of P,
    # denoted by lift(P), to be the tableau of key shape obtained by raising each entry in the first column of P
    # until it equals its row index, and, once columns 1 through c − 1 have been lifted, raising entries in column
    # c from top to bottom, maintaining their relative order, placing each entry in the highest available row such
    # that there is an entry in column c − 1 that is strictly smaller.

    #     Definition 5.6 ([Assa]). For ρ a reduced expression, define the weak insertion tableau Pb(ρ) by Pb(ρ) =
    # lift(P(ρ)), where P(ρ) is the insertion tableau under the Edelman–Greene insertion. In addition, define the
    # weak recording tableau Qb(ρ) to be the unique standard key tableau of the same key shape as Pb(ρ) such that

    #     Theorem 5.11. The operators fi and ei for 1 6 i < n define a Demazure crystal structure on RFC(w).
    # More precisely,
    # RFC(w) ∼=
    # [
    # r∈RFC(w)
    # eir=0 ∀16i<n
    # Bw(r)(wt(r)),
    # where w(r) is the shortest permutation that sorts sh(Pb(r)).

    #     Definition 5.6 ([Assa]). For ρ a reduced expression, define the weak insertion tableau Pb(ρ) by Pb(ρ) =
    # lift(P(ρ)), where P(ρ) is the insertion tableau under the Edelman–Greene insertion. In addition, define the
    # weak recording tableau Qb(ρ) to be the unique standard key tableau of the same key shape as Pb(ρ) such that

    #     Following [Assa], for P a semi-standard Young tableau with strictly increasing rows, define the lift of P,
    # denoted by lift(P), to be the tableau of key shape obtained by raising each entry in the first column of P
    # until it equals its row index, and, once columns 1 through c − 1 have been lifted, raising entries in column
    # c from top to bottom, maintaining their relative order, placing each entry in the highest available row such
    # that there is an entry in column c − 1 that is strictly smaller.

    def raising_operator(self, row):
        # RF word is just the RC word backwards
        if row >= len(self):
            return None
        row_i = [*self[row - 1]]
        row_ip1 = [*self[row]]

        # pair the letters
        pairings = []
        unpaired = []
        unpaired_b = [*row_ip1]

        for letter in row_i:
            st = [letter2 for letter2 in unpaired_b if letter2 > letter]
            if len(st) == 0:
                unpaired.append(letter)
            else:
                pairings.append((letter, min(st)))
                unpaired_b.remove(min(st))
        if len(unpaired_b) == 0:
            return None
        a = max(unpaired_b)
        s = 0
        while a + s + 1 in row_ip1:
            s += 1

        if a + s < row:
            return None
        new_row_ip1 = [let for let in row_ip1 if let != a]
        new_row_i = sorted([a + s, *row_i], reverse=True)
        ret_rc = type(self)([*self[: row - 1], tuple(new_row_i), tuple(new_row_ip1), *self[row + 1 :]])
        if ret_rc.perm != self.perm:
            return None
        return ret_rc

    # preserves plactic class, but not Coxeter-Knuth
    # preserves crystal structure. Decomposes RC graphs into key polynomials
    def vertical_cut(self, row):
        if row < 0:
            raise ValueError("Row out of range")
        if row >= len(self):
            return self, RCGraph()
        front = type(self)([*self[:row]])
        front = front.extend(max(len(self), len(front.perm.trimcode)) - row)
        flen = len(front)
        for _ in range(flen - row):
            front = front.zero_out_last_row()
        if row == len(self):
            back = type(self)()
        else:
            back = self.rowrange(row, len(self))
        return (front, back)

    def right_zero_act(self):
        # NOTE THAT THIS IS STILL USING THE OLD METHOD
        if self.perm.inv == 0:
            return {type(self)([*self, ()])}

        if self in RCGraph._z_cache:
            return RCGraph._z_cache[self]

        up_perms = ASx(self.perm, len(self)) * ASx(uncode([0]), 1)

        rc_set = set()

        for perm, _ in up_perms.keys():
            for rc in type(self).all_rc_graphs(perm, len(self) + 1, weight=(*self.length_vector, 0)):
                if rc.zero_out_last_row() == self:
                    rc_set.add(rc)

        RCGraph._z_cache[self] = rc_set
        return rc_set

    def __hash__(self):
        return hash(tuple(self))

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

    def exchange_property(self, descent, left=False, return_row=False):
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

    @cache
    def left_to_right_inversion(self, index):
        return self.right_root_at(*self.left_to_right_inversion_coords(index))

    @cache
    def left_to_right_inversion_coords(self, index):
        if index < 0 or index >= len(self.perm_word):
            raise ValueError(f"Index {index} out of range {self.perm.inv}")
        index_find = 0
        for i in range(len(self)):
            index_find2 = index_find + len(self[i])
            if index_find2 > index:
                break
            index_find = index_find2

        return (i + 1, self[i][(index - index_find)] - i)

    @classmethod
    def principal_rc(cls, perm, length):
        cd = perm.trimcode
        graph = []
        for i in range(len(cd)):
            row = tuple(range(i + cd[i], i, -1))
            graph.append(row)
        graph = [*graph, *[()] * (length - len(graph))]
        return cls(graph)

    @cached_property
    def p_tableau(self):
        return self.edelman_greene()[0]

    @cached_property
    def q_tableau(self):
        return self.edelman_greene()[1]

    @cached_property
    def weight_tableau(self):
        if self.is_highest_weight:
            tb = Plactic.yamanouchi(self.p_tableau.shape)
            trimmed_lv = list(self.length_vector)
            while len(trimmed_lv) > 0 and trimmed_lv[-1] == 0:
                trimmed_lv.pop()
            trimmed_lv = tuple(trimmed_lv)
            assert tb.shape == trimmed_lv, f"{tb.shape=}, {trimmed_lv=}"
            return tb
        rc_hw, raise_seq = self.to_highest_weight()
        w_tab = rc_hw.weight_tableau
        tb = w_tab.reverse_raise_seq(raise_seq)
        assert tb is not None, f"Could not reverse raise seq {raise_seq} on {w_tab=} {rc_hw=} {self=}"
        return tb

    # THE ZERO MAKES SCHUB PROD
    @cache
    def prod_with_rc(self, other):
        if self.perm.inv == 0:
            return {type(self)([*self, *other.shiftup(len(self))]): 1}
        num_zeros = max(len(other), len(other.perm))
        assert len(self.perm.trimcode) <= len(self), f"{self=}, {self.perm=}"
        base_rc = self
        buildup_module = {base_rc: 1}

        for _ in range(num_zeros):
            new_buildup_module = {}
            for rc, coeff in buildup_module.items():
                new_buildup_module = add_perm_dict(new_buildup_module, dict.fromkeys(rc.right_zero_act(), coeff))
            buildup_module = new_buildup_module
        ret_module = {}

        for rc, coeff in buildup_module.items():
            new_rc = type(rc)([*rc[: len(self)], *other.shiftup(len(self))])
            assert len(new_rc) == len(self) + len(other)
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

    def act(self, p):
        pm = self.perm
        elem = FA(pm, len(self))
        bumpup = FA(uncode([p]), 1) * elem
        ret = set()
        for k, v in bumpup.items():
            perm2 = k[0]
            new_row = [pm[i] for i in range(max(len(pm), len(perm2))) if pm[i] == perm2[i + 1]]
            new_row.sort(reverse=True)
            nrc = type(self)([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in self]])
            assert nrc.perm == perm2
            ret.add(nrc)
        assert ret == self.iterative_act(p), f"{ret=}\n{self.iterative_act(p)=}"
        return ret

    def iterative_act(self, p, insert=True):
        if p == 0:
            if insert:
                return {type(self)([(), *[tuple([row[i] + 1 for i in range(len(row))]) for row in self]])}
            return {type(self)([*self])}
        last = self.iterative_act(p - 1, insert=insert)
        ret = set()
        for rc in last:
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
                    new_rc = type(rc)([tuple(new_top_row), *rc[1:]])
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

    @property
    def rows(self):
        return len(self)

    @property
    def width(self):
        return self.cols

    @property
    def height(self):
        return self.rows

    @property
    def compatible_sequence(self):
        seq = []
        for i in range(len(self)):
            for _ in range(len(self[i])):
                seq.append(i + 1)
        return tuple(seq)

    @property
    def cols(self):
        return max(1, *[self[i][0] - i if len(self[i]) > 0 else 0 for i in range(len(self))]) if len(self) > 0 else 0

    def leibniz_rep(self):
        if len(self) == 0:
            return ()
        w0 = Permutation.w0(len(self) + 1)
        the_perm = (~self.perm) * w0
        cut_rc = self.shiftcut()
        return (*cut_rc.leibniz_rep(), the_perm)

    @classmethod
    @cache
    def all_hw_rcs(cls, perm, length):
        ret = set()
        for rc in cls.all_rc_graphs(perm, length):
            rc_hw, _ = rc.to_highest_weight()
            if rc_hw not in ret:
                ret.add(rc_hw)
        return ret

    @classmethod
    @cache
    def all_lw_rcs(cls, perm, length):
        ret = set()
        for rc in cls.all_rc_graphs(perm, length):
            rc_lw, _ = rc.to_lowest_weight()
            if rc_lw not in ret:
                ret.add(rc_lw)
        return ret

    # @classmethod
    # def from_leibniz_rep(cls, rep):
    #     from schubmult.utils.schub_lib import elem_sym_perms
    #     if len(rep) == 0:
    #         return cls(())
    #     new_rc = [()]*len(rep)
    #     for i in range(len(rep)):
    #         for j in range(i + 1):
    #             if i == 0 and rep[j] == 1:
    #                 new_rc[0] = (*new_rc[0], len(rep))
    #             if i > 0:
    #                 prev_perm
    #     return cls(new_rc_rows)

    def shiftcut(self):
        cut_rc = RCGraph([tuple([a for a in row if a > i]) for i, row in enumerate(self.shiftup(-1)[:-1])])
        return cut_rc

    def divdiff_desc(self, desc):
        ret = set()
        the_rc = self
        rc, row = the_rc.exchange_property(desc, return_row=True)
        rc = rc.normalize()
        if row != desc:
            return ret
        if rc.raising_operator(desc) is not None:
            return ret
        ret.add(rc)
        while rc is not None:
            rc = rc.lowering_operator(desc)
            if rc is not None:
                ret.add(rc)
        return ret

    def divdiff_perm(self, u):
        v = self.perm
        perm2 = v * (~u)
        if perm2.inv != v.inv - u.inv:
            return set()
        # return perm2
        ret = {self}
        working_perm = u
        while working_perm.inv > 0:
            working_set = set()
            desc = max(working_perm.descents()) + 1
            working_perm = working_perm.swap(desc - 1, desc)
            for the_rc in ret:
                assert desc - 1 in the_rc.perm.descents()
                working_set.update(the_rc.divdiff_desc(desc))

            ret = working_set
        return ret

    def dualpieri(self, mu, w):
        from schubmult.rings.rc_graph_ring import RCGraphRing
        from schubmult.utils.schub_lib import pull_out_var

        rc_ring = RCGraphRing()
        if mu.inv == 0:
            return set({((), (), rc) for rc in self.divdiff_perm(w)})

        cycle = Permutation.cycle
        lm = (~mu).trimcode
        cn1w = (~w).trimcode
        if len(cn1w) < len(lm):
            return set()
        for i in range(len(lm)):
            if lm[i] > cn1w[i]:
                return set()
        c = Permutation([])
        for i in range(len(lm), len(cn1w)):
            c = cycle(i - len(lm) + 1, cn1w[i]) * c

        res = {((), (), self)}

        for i in range(len(lm)):
            res2 = set()
            for vlist, perm_list, self_0 in res:
                vp = self_0

                vpl_list = vp.divdiff_perm(cycle(lm[i] + 1, cn1w[i] - lm[i]))

                if len(vpl_list) == 0:
                    continue
                for vpl in vpl_list:
                    vl = pull_out_var(lm[i] + 1, vpl.perm)

                    if lm[i] + 1 > len(vpl.perm.trimcode):
                        rcs = {vpl}
                    else:
                        vpl_bottom, vpl_top = vpl.vertical_cut(lm[i])
                        rcs = rc_ring(vpl_bottom) * rc_ring(vpl_top.rowrange(1))

                    for vpl_new in rcs:
                        if vpl_new.perm not in {vv[-1] for vv in vl}:
                            continue
                        pw = tuple(next(vv[0] for vv in vl if vv[-1] == vpl_new.perm))
                        res2.add(((*vlist, tuple(pw)), (*perm_list, vpl.perm), vpl_new))

            res = res2
        if len(lm) == len(cn1w):
            return res
        res2 = set()
        for vlist, perm_list, self_0 in res:
            vp = self_0
            vpl_list = vp.divdiff_perm(c)
            if len(vpl_list) == 0:
                continue
            for vpl in vpl_list:
                res2.add((tuple(vlist), perm_list, vpl))
        return res2

    @staticmethod
    def divdiff_act_dict(dct, *s_list):
        ret = {**dct}
        for s in reversed(s_list):
            new_ret = {}
            for rc, coeff in ret.items():
                act_set = rc.divdiff_action(s)
                new_ret = add_perm_dict(new_ret, dict.fromkeys(act_set, coeff))
            ret = new_ret
        return ret

    def __getitem__(self, key):
        # FLIPPED FOR PRINTING
        if isinstance(key, int):
            return tuple(self)[key]
        if isinstance(key, tuple):
            i, j = key
            if not self.has_element(i + 1, self.cols - j):
                return None
            return i + self.cols - j
        is_slice = isinstance(key, slice)

        if is_slice:
            return tuple(tuple(self)[n] for n in range(len(self))[key])

        raise ValueError(f"Bad indexing {key=}")

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

    # WE CAN DO STUFF WITH THIS
    def weight_reflection(self, i):
        try:
            rc = self.crystal_reflection(i)
        except Exception:
            rc = None
        if rc is None:
            rc = self.extend(1).shiftup(1).crystal_reflection(i)
        return rc

    @property
    def inverse_crystal(self):
        return InverseRCGraph(self)


class InverseRCGraph(CrystalGraph):
    def __init__(self, base_graph):
        self.base_graph = base_graph

    @property
    def crystal_weight(self):
        return self.base_graph.transpose().crystal_weight

    def raising_operator(self, index):
        lowered = self.base_graph.transpose().lowering_operator(index)
        if lowered is None:
            return None
        return InverseRCGraph(lowered.resize(len(self.base_graph)))

    def lowering_operator(self, index):
        raised = self.base_graph.transpose().raising_operator(index)
        if raised is None:
            return None
        return InverseRCGraph(raised.resize(len(self.base_graph)))

    def crystal_length(self):
        return self.base_graph.transpose().crystal_length()
