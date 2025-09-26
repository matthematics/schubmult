from functools import cache
from itertools import zip_longest

from symengine import SympifyError

import schubmult.schub_lib.schub_lib as schub_lib
from schubmult.perm_lib import Permutation, uncode
from schubmult.symbolic import S, prod, sympify
from schubmult.utils.perm_utils import add_perm_dict

from .free_algebra import FreeAlgebra, FreeAlgebraElement
from .free_algebra_basis import SchubertBasis, WordBasis
from .nil_hecke import NilHeckeRing
from .tensor_ring import TensorRing, TensorRingElement

FAS = FreeAlgebra(basis=SchubertBasis)


def _multiline_join(str1, str2, joiner=" "):
    lines1 = str1.split("\n")
    lines2 = str2.split("\n")
    width1 = max(len(line) for line in lines1)
    width2 = max(len(line) for line in lines2)
    if len(lines1) > len(lines2):
        lines2 += [""] * (len(lines1) - len(lines2))
    elif len(lines2) > len(lines1):
        lines1 += [""] * (len(lines2) - len(lines1))
    initial_line = f"{lines1[0]:<{width1}}{joiner}{lines2[0]:<{width2}}"
    # print(f"{lines1=} {lines2=}")
    result = "\n".join([initial_line, *[f"{l1:<{width1}}{'':<{len(joiner)}}{l2:<{width2}}" for l1, l2 in zip_longest(lines1[1:], lines2[1:], fillvalue=f"{'':<{width2}}")]])
    # print(result)
    return result


def full_multiline_join(*args, joiner=" "):
    if len(args) == 0:
        return ""
    if len(args) == 1:
        return args[0]
    buildup = args[0]
    for arg in args[1:]:
        buildup = _multiline_join(buildup, arg, joiner=joiner)
    return buildup


def _value_dict(dct_of_keys, keytype=None):
    if keytype is None:
        raise ValueError("Must provide keytype")
    dct = {}
    for k, v in dct_of_keys.items():
        if v != 0 and k.value != 0:
            dct[keytype(k)] = v * k.value
    return dct


class KeyType:
    def __init__(self, *args, value=1, **kwargs):
        self.value = value

    def ring_act(self, elem): ...

    def polyvalue(self, x, y=None): ...


class ModuleType:
    def __getitem__(self, key):
        return self.value_dict[key]

    def __radd__(self, other):
        if int(other) == 0:
            return self
        return NotImplemented
    
    def get(self, key, default=None):
        return self.value_dict.get(key, default)

    def keytype(self, *args): ...

    def asdtype(self, cls):
        return sum([v * k.asdtype(cls) for k, v in self._dict.items()])

    def items(self):
        return self._dict.items()

    def keys(self):
        return list(self._dict.keys())

    @property
    def generic_key_type(self):
        return self._generic_key_type

    # ASSUME VALUE DICT ALWAYS NORMALIZED
    @property
    def value_dict(self):
        return self._dict
        # return _value_dict(self._dict, keytype=self.keytype)

    @property
    def ring(self):
        return self._ring

    def clone(self, *args):
        ...
        # if len(args) == 0:
        #     return self.__class__(self._dict, generic_key_type=self._generic_key_type, ring=self._ring, **kwargs)
        # return self.__class__(*args, generic_key_type=self._generic_key_type, ring=self._ring, **kwargs)

    def _empty_module(self):
        return self.clone({})

    def __init__(self, *args, generic_key_type=None, ring=None, _clone=True, **kwargs):
        start_dct = dict(*args)

        if generic_key_type is None:
            raise ValueError("Must provide generic_key_type")
        if ring is None:
            raise ValueError("Must provide ring")

        # if _clone:
        #     self._dict = start_dct
        #     self._ring = ring
        #     self._generic_key_type = generic_key_type
        #     return

        self._ring = ring
        self._generic_key_type = generic_key_type
        self._dict = _value_dict(start_dct, keytype=self.keytype)

        # for k, v in start_dct.items():
        #     if not isinstance(k, self._generic_key_type):
        #         raise ValueError(f"Bad key type: {type(k)}, expected {generic_key_type}")
        #     value = v
        #     value *= k.value
        #     k2 = self.keytype(k)
        #     self._dict[k2] = value

    def _addkeys(self, other):
        return self.clone(add_perm_dict(_value_dict(self.value_dict, keytype=self.keytype), _value_dict(other.value_dict, keytype=other.keytype)))

    def __add__(self, other):
        if isinstance(other, self.__class__):
            return self._addkeys(other)
        if isinstance(other, self._generic_key_type):
            return self._addkeys(self.clone({other: 1}))
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, self.__class__):
            return self._addkeys(self.clone({k: -v for k, v in other._dict.items()}))
        if isinstance(other, self._generic_key_type):
            return self._addkeys(self.clone({other: -1}))
        return NotImplemented

    def __neg__(self):
        return self.clone({k: -v for k, v in self._dict.items()})

    def __rmul__(self, other):
        if not hasattr(other, "ring"):
            other = sympify(other)
            return self.clone({k: v * other for k, v in self.value_dict.items()})
        mod = self._empty_module()
        for k, v in self._dict.items():
            mod += v * self.clone(self.keytype(k).ring_act(other))
        return mod

    def __mul__(self, other):
        if not hasattr(other, "ring"):
            other = sympify(other)
            return self.clone({k: v * other for k, v in self.value_dict.items()})
        mod = self._empty_module()
        for k, v in self._dict.items():
            mod += v * self.clone(self.keytype(k).right_ring_act(other))
        return mod

    def __matmul__(self, other):
        return TensorModule(self, other)

    def __str__(self):
        buildups = ""
        if len(self.keys()) == 0:
            return "<zero module>"
        for k, v in self._dict.items():
            st1= self.coeffify(v)
            if len(st1) > 0:
                to_add = _multiline_join(st1, str(k), joiner="*")
            else:
                to_add = str(k)
            if len(buildups) == 0:
                buildups = to_add
            else:
                buildups = _multiline_join(buildups, to_add, joiner=" + ")
        return buildups

    def polyvalue(self, x, y=None):
        return sum([v * k.polyvalue(x, y) for k, v in self._dict.items()])

    def coeffify(self, v):
        if v == 1:
            return ""
        return str(v) + "*"


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


class RCGraph(KeyType, UnderlyingGraph):
    def __eq__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return tuple(self) == tuple(other)

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

    def reverse_kogan_insert(self, descent, pair_sequence):
        # print("Here")
        pair_sequence = sorted(pair_sequence, key=lambda x: x[0])
        pair_dict = {}
        pair_dict_rev = {}
        for a, b in pair_sequence:
            assert a <= descent and descent < b
            pair_dict[a] = pair_dict.get(a, set())
            pair_dict[a].add(b)

        for a, b in pair_sequence:
            pair_dict_rev[b] = a


        def is_relevant_crossing(root, prm):
            #min_root = max(pair_dict.keys())
            if root[0] not in pair_dict:
                if root[0] not in pair_dict or  pair_dict_rev.get(root[0], 0) != pair_dict_rev.get(root[1], 0):
                    return False
                return True
            return root[1] in pair_dict[root[0]] and prm[root[0] - 1] > prm[root[1] - 1]

        # may have to add q, s or a_i, q
        def is_relevant_noncrossing(root):
            top, bottom = max(root), min(root)
            return (root[0] <= descent and descent < root[1] and root[1] not in pair_dict_rev) or ((bottom in pair_dict_rev or (bottom not in pair_dict_rev and top in pair_dict_rev)) and bottom > descent)

        # Add this intersection. If we are in the first case, insert (s, q) into the sequence (ai, bi) in the rightmost position, such that ai’s remain nondecreasing in the sequence. ((s, q) are the rows where the two strands shown in Figure 3 originate.) If
        # we are in the second case, add (ai, q) just before where (a, bi) is in the sequence.

        new_rc = self
        i = len(new_rc) - 1
        index = 0
        # print(pair_dict)
        while i >= 0:
            build_graph = [*new_rc]
            did_any = False
            # print("starting with")
            # print(new_rc)
            if len(new_rc[i]) > 0:
                # print("Row", i)
                while index < max(new_rc[i]) - i:
                    #col_index = max(new_rc[i]) - i - index - 1
                    col_index = index
                    refl = col_index + i + 1
                    # index = col_index + 1 + 1
                    #assert index != 0 or new_rc.has_element(i + 1, col_index + 1)
                    print("Starting")
                    print(new_rc)
                    if new_rc.has_element(i + 1, col_index + 1):
                        print(f"Found at {col_index + 1}")
                        root = new_rc.right_root_at(i + 1, col_index + 1)
                        print(f"{root=}")
                        if is_relevant_crossing(root, new_rc.perm):
                            print("Relevant")
                            did_any = True
                            if root[0] in pair_dict:
                                pair_dict[root[0]].remove(root[1])
                                del pair_dict_rev[root[1]]
                                if len(pair_dict[root[0]]) == 0:
                                    del pair_dict[root[0]]
                            # may have to remember this root
                            build_graph[i] = tuple(a for a in build_graph[i] if a != refl)
                            start_spot = len(new_rc[i])
                            new_rc = RCGraph(build_graph)
                            print(f"removed")
                            #print(build_graph[i])
                            print(new_rc)
                            for index2 in range(index, start_spot):
                                #col_index2 = max(new_rc[i]) - i - index2 - 1
                                col_index2 = index2
                                refl = col_index2 + i + 1
                                if not new_rc.has_element(i + 1, col_index2 + 1):
                                    a, b = new_rc.right_root_at(i + 1, col_index2 + 1)
                                    root = (a, b)
                                    if is_relevant_noncrossing(root):
                                        # print(f"Putting it back")
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
                                        print("added back")
                                        print(new_rc)
                                        did_any = False
                    index += 1
                    if did_any:
                        break
            # print(f"{did_any=} {index=}")
            if not did_any:
                i -= 1
            else:
                i = len(new_rc) - 1
            index = 0
                # print(f"{i=}")
                # print(f"{index=}")
        assert len(pair_dict) == 0, f"{pair_dict=}, {pair_dict_rev=}, {new_rc=}, {pair_sequence=}"
        # print("Got")
        # print(new_rc)
        return new_rc

    def right_p_act(self, p):
        if len(self) == 0:
            return set((FA(p) * self).value_dict.keys())

        up_perms = ASx(self.perm, len(self)) * ASx(uncode([p]),1)

        old_rc_set = self.rowrange(1, len(self)).right_p_act(p)
        rc_set = set()

        for rc in old_rc_set:
            new_rc_set = FA(len(self[0])) * rc
            for rc2 in new_rc_set.value_dict.keys():
                if (rc2.perm, len(self) + 1) in up_perms.keys():
                    rc_set.add(rc2)
        return rc_set


    @staticmethod
    def full_rc_coproduct(perm, length):
        if length < len(perm.trimcode):
            raise ValueError("Length too short")
        #print(f"Trying {perm.trimcode}, {length}")
        if perm.inv == 0:
            if length == 0:
                return RCGraph([()*length])@RCGraph([()*length])
        if length == 1:
            return RCGraph.principal_rc(perm, length).coproduct()
        lower_perm = uncode(perm.trimcode[1:])
        elem = ASx(lower_perm, length - 1)
        lower_module1 = RCGraph.full_rc_coproduct(lower_perm, length - 1)
        #print("lower_module1")
        #print(lower_module1)
        ret_elem = RCGraph.full_rc_coproduct(uncode([perm.trimcode[0]]),1) * lower_module1

        ret_elem = ret_elem.clone({k: v for k, v in ret_elem.items() if k[0].perm <= perm and k[1].perm <= perm})

        if length == 1:
            return ret_elem
        keys = set(ret_elem.keys())
        # print(f"{repr(keys)=} {perm=}")
        up_elem = ASx(uncode([perm.trimcode[0]]), 1) * elem
        # print(f"{up_elem=}")
        for key, coeff in up_elem.items():
            if key[0] != perm:
                assert coeff == 1
                ret_elem -= RCGraph.full_rc_coproduct(key[0], length)
                # for (rc1_bad, rc2_bad), cff2 in try_lr_module(key[0], length).items():
                #     keys2 = set(ret_elem.keys())
                #     for rc1, rc2 in keys2:
                #         if (rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm) and (rc1.length_vector() >= rc1_bad.length_vector() or rc2.length_vector() >= rc2_bad.length_vector()):
                #             try:
                #                 keys = set(keys2)
                #                 keys.remove((rc1, rc2))
                #             except KeyError:
                #                 # print(repr(keys))
                #                 raise
                #             break
        # print(f"Done {perm}")
        return ret_elem


    @cache
    def inversion_label(self, i, j):
        if i >= j:
            raise ValueError("i must be less than j")
        if self.perm[i] < self.perm[j]:
            raise ValueError("Not an inversion")
        index = 0
        for i0, row in enumerate(self):
            for j0, a in enumerate(row):
                if(self.perm.right_root_at(index, word=self.perm_word()) == (i + 1, j + 1)):
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

    def edelman_greene(self):
        from schubmult.perm_lib import NilPlactic

        word1 = []
        word2 = []
        index = 0
        evil_self = list(reversed([list(reversed(row)) for row in self]))
        for i in range(len(evil_self)):
            for a in evil_self[i]:
                to_insert = len(self) - i
                word1, word2 = NilPlactic.ed_insert_rsk(word1, word2, a, to_insert)
                index += 1
        P = Tableau(word1)
        Q = Tableau(word2)
        # reg._rc_graph = self
        return (P, Q)

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

    def __matmul__(self, other):
        if isinstance(other, RCGraph):
            return RCGraphModule({self: 1}, generic_key_type=self.__class__) @ RCGraphModule({other: 1}, generic_key_type=self.__class__)
        if isinstance(other, ModuleType):
            return RCGraphModule({self: 1}, generic_key_type=self.__class__) @ other
        return NotImplemented

    def __rmul__(self, other):
        return other * RCGraphModule({self: 1}, generic_key_type=self.__class__)

    def __mul__(self, other):
        return RCGraphModule({self: 1}, generic_key_type=self.__class__).__mul__(other)

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

    @classmethod
    @cache
    def all_rc_graphs(cls, perm, length=-1):
        if perm.inv == 0:
            return {RCGraph([()] * length if length >= 0 else [()])}
        if len(perm.trimcode) == 1:
            nrc = RCGraph((tuple(range(perm.code[0], 0, -1)),))
            if len(nrc) < length:
                nrc = RCGraph((*nrc, *tuple([()] * (length - len(nrc)))))
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
                ret.add(nrc)
        return ret

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

    def dualact(self, p):
        pm = ~self.perm
        elem = FAS(pm, len(self))
        bumpup = FAS(uncode([p]), 1) * elem
        ret = set()
        for k, v in bumpup.items():
            perm2 = k[0]
            new_row = [pm[i] for i in range(max(len(pm), len(perm2))) if pm[i] == perm2[i + 1]]
            new_row.sort()
            lst = [tuple([a + 1 for a in row]) for row in self]

            for index in range(max(len(self) + 1, *new_row)):
                if index < len(lst):
                    if index + 1 in new_row:
                        lst[index] = (*lst[index], index + 1)
                else:
                    if index + 1 in new_row:
                        lst += [(index + 1,)]
                    else:
                        lst += [()]
            nrc = RCGraph(lst)
            assert nrc.perm == ~perm2
            ret.add(nrc)
        return ret

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

    @cache
    def coproduct(self):
        if len(self) == 0:
            return self @ self
        if len(self) == 1:
            m = self.perm.inv
            mod = FA(m).coproduct() * (RCGraph()@RCGraph())
            # print("mod")
            # print(mod)
            return mod
        # if len(self) > len(self.perm.trimcode):
        #     return (self.rowrange(0, len(self.perm.trimcode))).coproduct() * RCGraph([()]*(len(self)-len(self.perm.trimcode))).coproduct()
        # if len(self.perm) - 1< len(self):
        #     forba = self.rowrange(0, len(self.perm) - 1).coproduct()
        #     ret_elem = 0
        #     for (rc1, rc2), coeff in forba.items():
        #         ret_elem += (RCGraph([*rc1, ()*(len(self) - len(rc1))]) @ RCGraph([*rc2, ()*(len(self) - len(rc2))]))
        #     return ret_elem
        if self.is_principal:
            lower_graph = self.rowrange(1, len(self))
            lower_module1 = lower_graph.coproduct()
            #print("lower_module1")
            #print(lower_module1)
            ret_elem = FA(len(self[0])).coproduct() * lower_module1

            

            # print(f"{repr(keys)=} {perm=}")
            
            #print("ret_elem")
            #print(ret_elem)

            up_elem = RCGraph.one_row(len(self[0])) * lower_graph
            for key, coeff in up_elem.items():
                if key.perm != self.perm:
                    assert coeff == 1
                    assert key.perm != self.perm
                    #print("key")
                    #print(key)
                    #print("self")
                    #print(self)
                    ret_elem -= key.coproduct()
                    # for (rc1_bad, rc2_bad), cff2 in try_lr_module(key[0], length).items():
                    #     keys2 = set(ret_elem.keys())
                    #     for rc1, rc2 in keys2:
                    #         if (rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm) and (rc1.length_vector() >= rc1_bad.length_vector() or rc2.length_vector() >= rc2_bad.length_vector()):
                    #             try:
                    #                 keys = set(keys2)
                    #                 keys.remove((rc1, rc2))
                    #             except KeyError:
                    #                 # print(repr(keys))
                    #                 raise
                    #             break
            # print(f"Done {perm}")
            #print(ret_elem)
        else:
            prin_elem = RCGraph.principal_rc(self.perm, len(self)).coproduct()
            ret_elem = 0
            rc_selection = set((FA(*self.length_vector()).coproduct() * (RCGraph()@RCGraph())).value_dict.keys())
            for (rc1, rc2), coeff in prin_elem.items():
                if coeff == 0:
                    continue
                new_set = set(rc_selection)
                for rc01, rc02 in sorted(new_set):
                    if rc01.lehmer_partial_leq(rc1) and rc02.lehmer_partial_leq(rc2):
                        ret_elem += rc01 @ rc02
                        rc_selection.remove((rc01, rc02))
                        break
        ret_elem = ret_elem.clone({k: v for k, v in ret_elem.items() if k[0].perm <= self.perm and k[1].perm <= self.perm})
        return ret_elem
    

    @classmethod
    def principal_rc(cls, perm, length):
        cd = perm.trimcode
        graph = []
        for i in range(len(cd)):
            row = tuple(range(i + cd[i], i, -1))
            graph.append(row)
        graph = [*graph, *[()] * (length - len(graph))]
        return cls(graph)

    def prod_with_rc(self, other):
        buildup_module = 1*self
        num_zeros = len(other)
        for row in other:
            new_buildup_module = RCGraphModule()
            for rc, coeff in buildup_module.items():
                new_buildup_module += RCGraphModule(dict.fromkeys(rc.right_p_act(len(row)), coeff))
            buildup_module = new_buildup_module
        ret_module = RCGraphModule()
        up_perms = ASx(self.perm, len(self)) * ASx(other.perm, 1)
        vset = {RCGraph([*rc[:len(self)], () * len(other)]) for rc in buildup_module.value_dict.keys() if rc.perm.inv == len(rc.perm_word())}
        for rc in vset:
            new_rc = RCGraph([*rc[:len(self)], *[tuple([a + len(self) for a in row]) for row in other]])
            if new_rc.perm.inv == len(new_rc.perm_word()):
                if (new_rc.perm, len(new_rc)) in up_perms.keys():
                    coeff = up_perms[(new_rc.perm, len(new_rc))]
                    ret_module += coeff * new_rc
        return ret_module

    def ring_act(self, elem):
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
            # top row has some stuff
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

    def as_str_lines(self):
        lines = []
        for i in range(len(self)):
            row = self[i]
            splug_row = [*row, i]
            row = [str(splug_row[j]) + ("  ") * (splug_row[j] - splug_row[j + 1] - 1) for j in range(len(splug_row) - 1)]
            line = ""
            line += " ".join([str(r) for r in row])
            lines += [line]
        if len(lines) == 0:
            lines = [""]
        lines2 = []
        ml = max([len(line) for line in lines])
        if len(lines[0]) < ml:
            lines[0] = " " * (ml - len(lines[0])) + lines[0]
        for line in lines:
            lines2 += [" " * (ml - len(line)) + line]
        return lines2

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

    def __str__(self):
        lines = self.as_str_lines()
        lines2 = [line for line in lines]
        return "\n".join(lines2)

    def __repr__(self):
        return f"{self.__class__.__name__}(" + ", ".join([repr(k) for k in self]) + ")"

    def __hash__(self):
        return hash(tuple(self))

    def __lt__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
           #return self.perm.bruhat_leq(other.perm) and self.perm != other.perm
        try:
            for i in range(self.perm.inv):
                a, b = self.perm.right_root_at(i)
                # if self.inversion_label(a - 1, b - 1) < other.inversion_label(a - 1, b - 1):
                #     return True
                # if self.inversion_label(a - 1, b - 1) > other.inversion_label(a - 1, b - 1):
                #     return False
                
                if self.lehmer_label(a - 1, b - 1) < other.lehmer_label(a - 1, b - 1):
                    return True
                if self.lehmer_label(a - 1, b - 1) > other.lehmer_label(a - 1, b - 1):
                    return False
        except ValueError:
            return False
        return False

    def __le__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return self < other or self == other


class Tableau(RCGraph):

    def __mul__(self, other):
        from schubmult.perm_lib import Plactic
        if not isinstance(other, Tableau):
            return NotImplemented
        # Plactic
        working_tab = Plactic(tuple([tuple([-a for a in row]) for row in self]))
        for row in reversed(other):
            for a in row:
                working_tab = working_tab.rs_insert(-a)
        working_tab2 = (tuple([-a for a in row]) for row in working_tab._word)
        return Tableau(working_tab2)

    def __le__(self, other):
        if not isinstance(other, Tableau):
            return NotImplemented
        return self.length_vector() < other.length_vector() or (self.length_vector() == other.length_vector() and tuple(self) <= tuple(other))

    def __lt__(self, other):
        if not isinstance(other, Tableau):
            return NotImplemented
        return self <= other and self != other

    def __new__(cls, *args, value=1):
        new_args = tuple(tuple(arg) for arg in args)
        obj = Tableau.__xnew_cached__(cls, *new_args)
        obj.value = value
        return obj

    def __init__(self, *args, value=1):
        pass

    @staticmethod
    @cache
    def __xnew_cached__(_class, *args):
        return Tableau.__xnew__(_class, *args)

    @staticmethod
    def __xnew__(_class, *args):
        return UnderlyingGraph.__new__(_class, *args)

    def as_str_lines(self):
        lines = []
        # print(tuple(self))
        if len(self) == 0:
            return [""]
        for i in range(len(self)):
            row = self[i]
            line = ""
            line += " ".join([str(r) for r in row])
            lines += [line]
        if len(lines) == 0:
            lines = [""]
        lines2 = []
        ml = max([len(line) for line in lines])
        if len(lines[0]) < ml:
            lines[0] = lines[0] + " " * (ml - len(lines[0]))
        for line in lines:
            lines2 += [line + " " * (ml - len(line))]
        return lines2

    @property
    def perm(self):
        perm = Permutation([])
        for row in self:
            for p in reversed(row):
                perm = perm.swap(p - 1, p)
        return perm
    
    def maxvalue(self):
        mv = 0
        for row in self:
            for a in row:
                mv = max(mv, a)
        return mv


class RCGraphEG(RCGraph):

    def __new__(cls, *args, value=1):
        if len(args) == 0:
            return RCGraph.__new__(cls, value=value)
        if isinstance(args[0], RCGraphEG):
            return RCGraph.__new__(cls, *args, value=value)
        if isinstance(args[0], (tuple, list)):
            return cls.from_rc_graph(RCGraph(args[0], value=value))
        raise ValueError(f"Don't know what to do with {args}")

    def __init__(self, *args, value=1):
        if len(args) == 0:
            self.value = value
            self._rc_graph = RCGraph()
            self.P = Tableau()
            self.Q = Tableau()
            return
        # if len(args) == 1:
        if isinstance(args[0], RCGraphEG):
            self.P, self.Q = args[0].P, args[0].Q
            self._rc_graph = args[0]._rc_graph
            self.value = value
            return
        return

        # elif len(args) == 2:
        #     P, Q = args
        #     self._length = len(P)
        # else:
        #     raise ValueError("Must provide 1 or 2 arguments")
        # if isinstance(P, tuple):
        #     P = Tableau(P)
        # if isinstance(Q, tuple):
        #     Q = Tableau(Q)
        # if not isinstance(P, Tableau) or not isinstance(Q, Tableau):
        #     raise ValueError("Must provide tableaux")
        # if P.length_vector() != Q.length_vector():
        #     raise ValueError("Tableaux must have same shape")
        # self.P = P
        # self.Q = Q
        # self.value = value
        #self._length = 

    @classmethod
    def from_rc_graph(cls, rc_graph):
        sporg = cls(rc_graph.edelman_greene())
        sporg._rc_graph = rc_graph
        sporg.value = rc_graph.value
        return sporg

    def __len__(self):
        return len(self.P)

    def __str__(self):
        return full_multiline_join("\n".join(self.P.as_str_lines()), "\n".join(self.Q.as_str_lines()), joiner=" | ")

    def __repr__(self):
        return f"RCGraphEG({repr(self.P)}, {repr(self.Q)})"

    def __eq__(self, other):
        if not isinstance(other, RCGraphEG):
            return NotImplemented
        return self.P == other.P and self.Q == other.Q

    def __hash__(self):
        return hash((self.P, self.Q))

    @property
    def perm(self):
        return self.P.perm

    def iterative_act(self, p):
        raise NotImplementedError("Use act instead")

    # def ring_act(self, elem):
    #     if isinstance(elem, FreeAlgebraElement):
    #         wd_dict = elem.change_basis(WordBasis)
    #         ret = {}
    #         for k, v in wd_dict.items():
    #             # print(f"{k=} {v=}")
    #             acted_element = {self: v}
    #             # not reversed
    #             for a in k:
    #                 acted_element2 = {}
    #                 for k2, v2 in acted_element.items():
    #                     # print(f"{dict.fromkeys(k2.act(a), v2)=}")
    #                     acted_element2 = add_perm_dict(acted_element2, dict.fromkeys(k2.act(a), v2))
    #                 acted_element = acted_element2
    #             # print(f"{acted_element=}")
    #             ret = add_perm_dict(ret, acted_element)
    #         # print(f"{ret=}")
    #         return ret
    #     raise ValueError(f"Cannot act by {type(elem)} {elem=}")

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
    #     return RCGraphEG(Tableau(word1), Tableau(word2))

    def act(self, p):
        from schubmult.perm_lib import NilPlactic

        # pm = self.perm
        # elem = FAS(pm)
        # bumpup =  elem * FAS(uncode([p]), 1)
        # ret = set()

        # for k, v in bumpup.items():
        #     perm2 = k[0]
        #     new_row = [pm[i] for i in range(max(len(pm), len(perm2))) if pm[i] == (perm2)[i + 1]]
        #     new_row.sort(reverse=True)
        #     newP = Tableau([tuple([a + 1 for a in row]) for row in self.P])
        #     newQ = Tableau([tuple([a + 1 for a in row]) for row in self.Q])
        #     newP_tup, newQ_tup = tuple(newP), tuple(newQ)
        #     for char in new_row:
        #         newP_tup, newQ_tup = NilPlactic.ed_insert_rsk(newP_tup, newQ_tup, char, 1)
        #     newP = Tableau(newP_tup)
        #     newQ = Tableau(newQ_tup)
        #     nrc = RCGraphEG(newP, newQ)
        #     # print(nrc)
        #     #assert nrc.perm == perm2, f"{nrc.perm=} {perm2=}"
        #     ret.add(nrc)
        # assert ret == self.iterative_act(p), f"{ret=}\n{self.iterative_act(p)=}"
        # pm = ~self.perm
        # elem = FAS(pm)
        # bumpup =   FAS(uncode([p]), 1)*elem
        # ret = set()
        # print(f"INSERTING {self.perm.trimcode=}")
        # for k, v in bumpup.items():
        #     perm2 = k[0]
        #     new_row = [(pm)[i] for i in range(max(len(pm), len(perm2))) if (pm)[i] == (perm2)[i + 1]]
        #     new_row.sort()
        #     assert len(new_row) == p
        #     newP = Tableau([tuple([a for a in row]) for row in self.P])
        #     newQ = Tableau([tuple([a for a in row]) for row in self.Q])
        #     newP_tup, newQ_tup = tuple(newP), tuple(newQ)
        #     for char in new_row:
        #         newP_tup, newQ_tup = NilPlactic.ed_insert_rsk(newP_tup, newQ_tup, char, self._length + 1)
        #     newP = Tableau(newP_tup)
        #     newQ = Tableau(newQ_tup)
        #     nrc = RCGraphEG(newP, newQ)
        #     nrc._length = self._length + 1
        #     # print(nrc)
        #     assert nrc.perm.inv == self.perm.inv + p, f"{nrc.perm.trimcode=} {self.perm.trimcode=}"
        #     #assert nrc.perm == perm2, f"{nrc.perm=} {perm2=}"
        #     ret.add(nrc)
        # return ret

        st = self._rc_graph.act(p)
        return {RCGraphEG.from_rc_graph(rc_graph) for rc_graph in st}

    def __lt__(self, other):
        if not isinstance(other, RCGraphEG):
            return NotImplemented
        if self.P != other.P:
            return self.P < other.P
        return self.Q < other.Q

    def __le__(self, other):
        return self == other or self < other


class RCGraphModule(ModuleType):
    def clone(self, *args):
        if len(args) == 0:
            return RCGraphModule(self._dict, generic_key_type=self._generic_key_type, ring=self._rings)
        return RCGraphModule(*args, generic_key_type=self._generic_key_type, ring=self._ring)

    def keytype(self, k, value=1):
        return self.generic_key_type(k, value=value)

    def __init__(self, *args, generic_key_type=RCGraph, ring=FreeAlgebra(WordBasis), **kwargs):
        super().__init__(*args, generic_key_type=generic_key_type, ring=ring, **kwargs)

    def __mul__(self, other):
        if isinstance(other, RCGraph):
            result = RCGraphModule()
            for rc_graph, coeff in self._dict.items():
                result += coeff * rc_graph.prod_with_rc(other)
            return result
        if isinstance(other, RCGraphModule):
            result = RCGraphModule()
            for rc_graph2, coeff2 in other._dict.items():
                result += coeff2 * (self * rc_graph2)
            return result
        return super().__mul__(other)

    def schubvalue(self, sring):
        ret = S.Zero
        for k, v in self._dict.items():
            ret += v * sring(k.perm)
        return ret

    def _produce_poly_dict(self, poly, genset, length):
        from .variables import genset_dict_from_expr

        dct = genset_dict_from_expr(poly, genset)
        dct2 = {}
        for vec, coeff in dct.items():
            if len(vec) > length:
                return {}
            dct2[tuple([*vec] + [0] * (length - len(vec)))] = coeff
        return dct2

    def coproduct(self):
        res = 0
        for rc, coeff in self._dict.items():
            res += coeff * rc.coproduct()
        return res

    def apply(self, poly, genset, length):
        from . import ASx

        dct = self._produce_poly_dict(poly, genset, length)

        if len(dct) == 0:
            return RCGraphModule()
        res = RCGraphModule()

        for rc, coeff_rc in self._dict.items():
            fa_elem = ASx(rc.perm, length).change_basis(WordBasis)
            for vec, coeff_vec in fa_elem.items():
                res += coeff_vec * coeff_rc * dct.get(vec, 0) * RCGraphModule({rc: 1})

        return res

    def apply_product(self, poly1, poly2, genset, length):
        from . import ASx

        dct1 = self._produce_poly_dict(poly1, genset, length)
        dct2 = self._produce_poly_dict(poly2, genset, length)

        res = RCGraphModule()

        for rc, coeff_rc in self._dict.items():
            fa_elem = ASx(rc.perm, length).change_basis(WordBasis).coproduct()
            for (vec1, vec2), coeff_vec in fa_elem.items():
                res += coeff_vec * coeff_rc * dct1.get(vec1, 0) * dct2.get(vec2, 0) * RCGraphModule({rc: 1})

        return res

    @classmethod
    def schub_module(cls, perm):
        return cls(dict.fromkeys(RCGraph.all_rc_graphs(perm), 1))

    @classmethod
    def schub_injection(cls, poly, genset):
        from .schubert_ring import SingleSchubertRing

        ring = SingleSchubertRing(genset)
        dct = ring(poly)
        ret = cls()
        for perm, coeff in dct.items():
            ret += coeff * RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(perm), 1))
        return ret

    @classmethod
    def from_poly(cls, poly, genset, length):
        from .free_algebra import FreeAlgebra
        from .variables import genset_dict_from_expr

        FA = FreeAlgebra()
        md = cls()
        md0 = cls({RCGraph(): 1})
        gdct = genset_dict_from_expr(poly, genset)
        for key, val in gdct.items():
            if len(key) > length:
                raise ValueError("Key too long")
            if len(key) < length:
                key = (*key, *(0,) * (length - len(key)))
            md += val * FA(*key) * md0
        return md

    # def __new__(cls, *args):
    #     obj = ModuleType.__new__(cls, *args)
    #     obj.keytype = type("RCGraph", (RCGraph,), {"module": obj})
    #     return obj

    def polyvalue(self, x, y=None):
        ret = S.Zero
        for k, v in self._dict.items():
            ret += v * k.polyvalue(x, y)
        return ret

    def as_nil_hecke(self, x, y=None):
        ret = S.Zero
        for k, v in self._dict.items():
            ret += v * k.as_nil_hecke(x, y)
        return ret

    def transpose(self):
        return RCGraphModule({graph.transpose(): v for graph, v in self._dict.items()})

    # def as_str_lines(self):

    # return full_multiline_join(full_multiline_join(k.as_str_lines()) + f"  (coeff {v})" for k, v in self._dict.items()])
    # # if len(self.keys()) == 0:
    # #     return ["0"]
    # # lines = [""]
    # # first = True
    # # for k, v in self._dict.items():
    # #     lines2 = k.as_str_lines()
    # #     if len(lines) < len(lines2):
    # #         upstr = ""
    # #         if len(lines[0]) > 0:
    # #             upstr = " " * len(lines[0])
    # #         lines += [upstr] * (len(lines2) - len(lines))
    # #     padlen = 0
    # #     for i in range(len(lines2)):
    # #         coeffstr = ""
    # #         if not first:
    # #             if i == 0:
    # #                 coeffstr += " + "
    # #             else:
    # #                 coeffstr += "   "
    # #         if i == 0:
    # #             if v == -1:
    # #                 coeffstr += "-"
    # #             elif v != 1:
    # #                 coeffstr += str(v) + " * "
    # #             padlen = len(coeffstr)
    # #         else:
    # #             coeffstr = " " * padlen

    # #         lines[i] += coeffstr + lines2[i]
    # #     first = False
    # # return lines
    # buildup_str = ""
    # for k,v in self._dict.items():
    #     lines2 = k.as_str_lines()
    #     width1 = max(len(line) for line in lines2)
    #     width2 = len(str(v)) + 3
    #     lines2 = _multiline_join(str(v) + " *", "\n".join(lines2), width2, width1)
    #     if len(lines2) == 0:
    #         lines2 = "0"
    #     yield lines2
    #     yield ""

    # def __str__(self):
    #     return "\n".join(self.as_str_lines())


# class DualRCGraphModule(RCGraphModule):
#     def __new__(cls, *args):
#         return RCGraphModule.__new__(cls, *args)

#     def __init__(self, *args):
#         pass

#     def __add__(self, other):
#         if isinstance(other, DualRCGraphModule):
#             return DualRCGraphModule(add_perm_dict(self, other))
#         return NotImplemented

#     def __neg__(self):
#         return DualRCGraphModule({k: -v for k, v in self._dict.items()})

#     def __sub__(self, other):
#         if isinstance(other, DualRCGraphModule):
#             return DualRCGraphModule(add_perm_dict(self, -other))
#         return NotImplemented

#     def __mul__(self, other):
#         if isinstance(other, FreeAlgebraElement):
#             wd_dict = other.change_basis(WordBasis)
#             ret = {}
#             for k0, v0 in self._dict.items():
#                 for k, v in wd_dict.items():
#                     addup = DualRCGraphModule({k0: v0 * v})
#                     for a in k:
#                         new_addup = DualRCGraphModule()
#                         for kr, vv in addup.items():
#                             new_addup += DualRCGraphModule(dict.fromkeys(kr.dualact(a), vv))
#                         addup = new_addup

#                     ret = add_perm_dict(ret, addup)
#             return DualRCGraphModule(ret)
#         try:
#             other = sympify(other)
#         except SympifyError:
#             return NotImplemented
#         return DualRCGraphModule({k: v * other for k, v in self._dict.items()})

#     def __rmul__(self, other):
#         try:
#             other = sympify(other)
#         except SympifyError:
#             return NotImplemented
#         return DualRCGraphModule({k: v * other for k, v in self._dict.items()})

#     def as_str_lines(self):
#         if len(self.keys()) == 0:
#             return ["0"]
#         lines = [""]
#         first = True
#         for k, v in self._dict.items():
#             lines2 = k.as_str_lines()
#             if len(lines) < len(lines2):
#                 upstr = ""
#                 if len(lines[0]) > 0:
#                     upstr = " " * len(lines[0])
#                 lines += [upstr] * (len(lines2) - len(lines))
#             padlen = 0
#             for i in range(len(lines2)):
#                 coeffstr = ""
#                 if not first:
#                     if i == 0:
#                         coeffstr += " + "
#                     else:
#                         coeffstr += "   "
#                 if i == 0:
#                     if v == -1:
#                         coeffstr += "-"
#                     elif v != 1:
#                         coeffstr += str(v) + " * "
#                     padlen = len(coeffstr)
#                 else:
#                     coeffstr = " " * padlen

#                 lines[i] += coeffstr + lines2[i] + "^^"
#             first = False
#         return lines

#     def pairing(self, other):
#         if not isinstance(other, RCGraphModule):
#             return NotImplemented
#         ret = S.Zero
#         for k1, v1 in self._dict.items():
#             for k2, v2 in other.items():
#                 if k1 == k2:
#                     ret += v1 * v2
#         return ret

#     def __call__(self, other):
#         return self.pairing(other)

#     def __str__(self):
#         return "\n".join(self.as_str_lines())


class RCGraphTensor(KeyType):
    def __lt__(self, other):
        return tuple(self) < tuple(other)

    def __le__(self, other):
        return tuple(self) <= tuple(other)

    def polyvalue(self, x, y=None):
        return self[0].polyvalue(x, y) * self[1].polyvalue(x, y)

    def asdtype(self, cls):
        return cls.dtype().ring.from_rc_graph_tensor(self)

    def __init__(self, *keys, modules=None, value=1):
        if modules is None:
            new_keys = []
            for key in keys:
                if isinstance(key, RCGraphTensor):
                    for key in keys:
                        assert isinstance(key, RCGraph)
                        new_keys.append(key.__class__(key, value=1))
                else:
                    new_keys.append(key)
            self._keys = tuple(new_keys)
        else:
            new_keys = []
            for key in keys:
                if isinstance(key, RCGraphTensor):
                    new_keys += [key._modules[i].keytype(key[i]) for i, k in enumerate(key)]
                else:
                    new_keys.append(key.__class__(key, value=1))
            new_value = value * prod([k.value for k in keys])
            super().__init__(value=new_value)
            self._modules = modules
            self._keys = tuple(new_keys)
            assert all(k.value == 1 for k in self._keys)
        #  #  # print(f"BOOGERS {self._keys=}")

    def ring_act(self, elem):
        assert self.value == 1
        if isinstance(elem, TensorRingElement) and len(elem.ring.rings) == len(self):
            ret = {}
            for elem_key, value in elem.items():
                modules = []
                for i, key in enumerate(elem_key):
                    graph = self[i]
                    modules.append(self._modules[i].clone(dict.fromkeys(graph.ring_act(elem.ring.rings[i](*key)), 1)))
                ret = add_perm_dict(ret, TensorModule(value * modules[0], *modules[1:], ring=elem.ring).value_dict)
            return ret
        raise TypeError("TOASNTOSN")

    def right_ring_act(self, elem):
        assert self.value == 1
        if isinstance(elem, TensorRingElement) and len(elem.ring.rings) == len(self):
            ret = {}
            for elem_key, value in elem.items():
                modules = []
                for i, key in enumerate(elem_key):
                    graph = self[i]
                    modules.append(self._modules[i].clone(dict.fromkeys(graph.right_ring_act(elem.ring.rings[i](*key)), 1)))
                ret = add_perm_dict(ret, TensorModule(value * modules[0], *modules[1:], ring=elem.ring).value_dict)
            return ret
        raise TypeError("TOASNTOSN")

    def as_str_lines(self):
        lines1 = self[0].as_str_lines()
        lines2 = self[1].as_str_lines()

        ml = max(len(lines1), len(lines2))
        if len(lines1) < ml:
            lines1 = [" " * len(lines1[0])] * (ml - len(lines1)) + lines1
        if len(lines2) < ml:
            lines2 = [" " * len(lines2[0])] * (ml - len(lines2)) + lines2

        lines = [lines1[0] + "  #  " + lines2[0]]
        for i in range(1, ml):
            lines += [lines1[i] + "     " + lines2[i]]

        return lines

    def __repr__(self):
        return "RCGraphTensor(" + ", ".join([repr(k) for k in self._keys]) + ")"

    def __str__(self):
        return full_multiline_join(*[str(k) for k in self._keys], joiner="  #  ")

    def __hash__(self):
        return hash(self._keys)

    def __getitem__(self, index):
        return self._keys[index]

    def __len__(self):
        return len(self._keys)

    def __iter__(self):
        return iter(self._keys)

    def __eq__(self, other):
        if not isinstance(other, (RCGraphTensor, tuple)):
            return False
        return (isinstance(other, tuple) and tuple(self) == other) or self._keys == other._keys


class TensorModule(ModuleType):
    def __mul__(self, other):
        if isinstance(other, RCGraphTensor) and len(other) == len(self._modules):
            result = 0
            for rc_graph_tensor, coeff in self._dict.items():
                result0 = None
                for i in range(len(rc_graph_tensor)):
                    firstten = rc_graph_tensor[i] * other[i]
                    if result0 is None:
                        result0 = firstten
                    else:
                        result0 = result0 @ firstten
                result += coeff * result0
            return result

        if isinstance(other, TensorModule) and len(other._modules) == len(self._modules):
            result = 0
            for rc_graph2, coeff2 in other._dict.items():
                result += coeff2 * (self * rc_graph2)
            return result
        return super().__mul__(other)
    # def apply_product(self, poly1, poly2, genset, length):
    #     res = TensorModule()
    #     for (rc1, rc2), coeff in self._dict.items():
    #         res += TensorModule.ext_multiply(RCGraphModule({rc1: coeff}).apply(poly1, genset, length), RCGraphModule({rc2: 1}).apply(poly2, genset, length))

    #     return res

    # def __new__(cls, *args, **kwargs):
    #     obj = ModuleType.__new__(cls, *args)
    #     obj.ringtype = type("TensorRing", (TensorRing,))
    #     obj.keytype = type("RCGraphTensor", (RCGraphTensor,), {"module": obj}})
    #     return obj

    # def __reduce__(self):
    #     return (self.__class__, (self._dict, self._generic_key_type, self._ring, self._modules))

    def filter_keys(self, lambd):
        return self.clone({k: v for k, v in self._dict.items() if lambd(k)})

    def keytype(self, *args):
        key = RCGraphTensor(*args, modules=self._modules)
        assert key.value == 1
        return key

    def clone(self, *args):
        return TensorModule(*args, generic_key_type=self._generic_key_type, ring=self._ring, _modules=self._modules, _clone=True)

    def _empty_module(self):
        return TensorModule(generic_key_type=self._generic_key_type, ring=self._ring, _modules=self._modules, _clone=True)

    def __init__(self, *args, generic_key_type=RCGraphTensor, _modules=None, ring=None, _clone=False):
        if _clone:
            self._modules = _modules
            super().__init__(*args, generic_key_type=generic_key_type, ring=ring, _clone=False)
            return
        _dict = {}
        # put together the keys and values
        # if len(args) == 0 or (len(args) == 1 and isinstance(args[0], dict)):
        #     super().__init__(*args, generic_key_type=generic_key_type, **kwargs)
        #     return
        modules = _modules
        rings = []
        if _modules is None:
            modules = []
            for arg in args:
                if isinstance(arg, ModuleType):
                    modules.append(arg)
                else:
                    raise ValueError(f"Bad argument type: {type(arg)}, expecting Modules")

            if len(modules) > 0:
                if len(modules) != len(args):
                    raise ValueError("Cannot mix modules and non-modules in TensorModule")
                new_modules = []
                for i, mod in enumerate(modules):
                    if isinstance(mod, TensorModule):
                        new_modules += [*mod._modules]
                        rings += [*mod.ring.rings]
                    else:
                        new_modules.append(mod)
                        rings.append(mod.ring)

                for i, mod in enumerate(new_modules):
                    if i == 0:
                        _dict = mod._dict
                        continue
                    new_dict = {}
                    for key1, val1 in _dict.items():
                        for key2, val2 in mod._dict.items():
                            new_key = RCGraphTensor(key1, key2, modules=new_modules[: i + 1])
                            new_dict[new_key] = new_dict.get(new_key, 0) + val1 * val2
                    _dict = new_dict
            self._dict = _dict
            self._modules = tuple(modules)
        else:
            self._modules = tuple(modules)
            rings = [mod.ring for mod in modules]
            _dict = {self.keytype(k): v for k, v in dict(*args).items()}
            self._dict = _dict

        super().__init__(_dict, generic_key_type=generic_key_type, ring=TensorRing(*rings), _clone=False)

    @classmethod
    def ext_multiply(cls, elem1, elem2, strict=False):
        ret = cls()
        for key, val in elem1.items():
            if strict:
                if not isinstance(key, RCGraph):
                    raise ValueError("Not RCGraph in ext_multiply")
            for key2, val2 in elem2.items():
                if not isinstance(key2, RCGraph):
                    raise ValueError("Not RCGraph in ext_multiply")
                ret += cls({RCGraphTensor(key, key2): val * val2})
        return ret

    # def __rmul__(self, other):
    #     return super().__rmul__(self, other)
    # if isinstance(other, TensorRingElement):
    #     ret = {}
    #     for k0, v0 in self._dict.items():
    #         for k, v in other.items():
    #             addup = TensorModule({k0: v * v0})
    #             new_addup = TensorModule()
    #             for kr, vv in addup.items():
    #                 elem1 = other.ring.rings[0](*k[0]) * RCGraphModule({kr[0]: 1})
    #                 if isinstance(other.ring.rings[1], TensorRing):
    #                     raise ValueError("Not implemented")
    #                 elem2 = other.ring.rings[1](*k[1]) * RCGraphModule({kr[1]: 1})

    #                 new_addup += vv * TensorModule.ext_multiply(elem1, elem2, strict=True)
    #             addup = new_addup

    #             ret = add_perm_dict(ret, addup)
    #     return TensorModule(ret)
    # try:
    #     other = sympify(other)
    #     return TensorModule({k: v * other for k, v in self._dict.items()})
    # except Exception:
    #     return NotImplemented


ASx = FreeAlgebra(SchubertBasis)
FA = FreeAlgebra(WordBasis)


@cache
def try_lr_module(perm, length=None):
    # print(f"Starting {perm}")
    if length is None:
        length = len(perm.trimcode)
    elif length < len(perm.trimcode):
        raise ValueError("Length too short")
    if perm.inv == 0:
        if length == 0:
            mod = RCGraph() @ RCGraph()
            #  #  # print(f"MOASA!! {mod=} {type(mod)=}")
            return mod
        return FA(*([0] * length)).coproduct() * (RCGraph() @ RCGraph())
    lower_perm = uncode(perm.trimcode[1:])
    elem = ASx(lower_perm, length - 1)
    lower_module1 = try_lr_module(lower_perm, length - 1)
    assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
    #  #  # print(f"Coproducting {ASx(uncode([perm.trimcode[0]]), 1).coproduct()=}")
    #  #  # print(ASx(uncode([perm.trimcode[0]]), 1).coproduct())
    #  #  # print("Going for it")
    #  #  # print(f"{type(lower_module1)=} {lower_module1=}")
    #  #  # print(f"{type(ASx(uncode([perm.trimcode[0]]), 1).coproduct())=}")
    ret_elem = ASx(uncode([perm.trimcode[0]]), 1).coproduct() * lower_module1
    #  #  # print(f"{ret_elem=}")
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"

    ret_elem = ret_elem.clone({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)})

    if length == 1:
        return ret_elem
    keys = set(ret_elem.keys())
    # print(f"{repr(keys)=} {perm=}")
    up_elem = ASx(uncode([perm.trimcode[0]]), 1) * elem
    # print(f"{up_elem=}")
    for key, coeff in up_elem.items():
        if key[0] != perm:
            assert coeff == 1
            for (rc1_bad, rc2_bad), cff2 in try_lr_module(key[0], length).items():
                keys2 = set(ret_elem.keys())
                for rc1, rc2 in keys2:
                    if (rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm) and (rc1.length_vector() >= rc1_bad.length_vector() or rc2.length_vector() >= rc2_bad.length_vector()):
                        try:
                            keys = set(keys2)
                            keys.remove((rc1, rc2))
                        except KeyError:
                            # print(repr(keys))
                            raise
                        break
    # print(f"Done {perm}")
    ret_elem = ret_elem.clone({k: v for k, v in ret_elem.items() if k in keys})
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {perm.trimcode=}"
    return ret_elem


@cache
def try_lr_module_inject(perm, length=None):
    # print(f"Starting {perm}")
    if length is None:
        length = len(perm.trimcode)
    elif length < len(perm.trimcode):
        raise ValueError("Length too short")
    if perm.inv == 0:
        if length == 0:
            mod = RCGraph() @ RCGraph()
            #  #  # print(f"MOASA!! {mod=} {type(mod)=}")
            return mod
        return FA(*([0] * length)).coproduct() * (RCGraph() @ RCGraph())
    lower_perm = uncode(perm.trimcode[1:])
    elem = ASx(lower_perm, length - 1)
    lower_module1 = try_lr_module_inject(lower_perm, length - 1)
    assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
    #  #  # print(f"Coproducting {ASx(uncode([perm.trimcode[0]]), 1).coproduct()=}")
    #  #  # print(ASx(uncode([perm.trimcode[0]]), 1).coproduct())
    #  #  # print("Going for it")
    #  #  # print(f"{type(lower_module1)=} {lower_module1=}")
    #  #  # print(f"{type(ASx(uncode([perm.trimcode[0]]), 1).coproduct())=}")
    ret_elem = ASx(uncode([perm.trimcode[0]]), 1).coproduct() * lower_module1
    #  #  # print(f"{ret_elem=}")
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"

    ret_elem = ret_elem.clone({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)})

    if length == 1:
        return ret_elem
    keys = set(ret_elem.keys())
    # print(f"{repr(keys)=} {perm=}")
    up_elem = ASx(perm, length).change_basis(WordBasis)
    elem = ASx(lower_perm, length - 1)
    up_elem = ASx(uncode([perm.trimcode[0]]), 1) * elem

    # print(f"{up_elem=}")
    def pad_vector(vec, length):
        if len(vec) < length:
            return (*vec, *(0,) * (length - len(vec)))
        return tuple(vec)

    loopkeys = [pad_vector(key[0].trimcode, length) for key, _ in up_elem.items() if key[0] != perm]
    loopkeys.sort()
    for key in loopkeys:
        code_key = key
        for (rc1_bad, rc2_bad), cff2 in (FA(*code_key).coproduct() * (RCGraph() @ RCGraph())).items():
            keys2 = set(ret_elem.keys())
            for rc1, rc2 in keys2:
                if (rc1.perm == rc1_bad.perm and rc2_bad.perm == rc2.perm) and (rc1.length_vector() >= rc1_bad.length_vector() or rc2.length_vector() >= rc2_bad.length_vector()):
                    try:
                        keys = set(keys2)
                        keys.remove((rc1, rc2))
                    except KeyError:
                        # print(repr(keys))
                        raise
                    break
    # print(f"Done {perm}")
    ret_elem = ret_elem.clone({k: v for k, v in ret_elem.items() if k in keys})
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {perm.trimcode=}"
    return ret_elem



# jump through hoops to make this polynomial

def fa_elem(perm, length):
    return FA(*((perm).trimcode), *((0,) * (length - len((perm).trimcode))))

@cache
def try_lr_module_biject(perm, length):
    from schubmult import schubmult_py

   
    rc_set = set((fa_elem(perm,length)*RCGraph()).value_dict.keys())
    #prin_rc
    consideration_set = [(k[0], k[1]) for k in ((fa_elem(perm,length).coproduct()*(RCGraph() @ RCGraph()))).value_dict.keys()]

    

    if len(rc_set) == 1:
        # dominant
        return consideration_set
    consider_dict = {}
    
    def lehmer_partial_pair(pair1, pair2):
        return pair1[0].lehmer_partial_leq(pair2[0]) or pair1[1].lehmer_partial_leq(pair2[1])
        #(pair1[0] == pair2[0] and pair1[1] <= pair2[1])

    for (rc1, rc2) in consideration_set:
        consider_dict[(rc1.perm, rc2.perm)] = consider_dict.get((rc1.perm, rc2.perm), set())
        consider_dict[(rc1.perm, rc2.perm)].add((rc1, rc2))
            # if any(lehmer_partial_pair((rc1, rc2), (rc1_dom, rc2_dom)) for (rc1_dom, rc2_dom) in dom_rcs if rc1_dom.perm == rc1.perm and rc2_dom.perm == rc2.perm) and rc1.perm <= perm and rc2.perm <= perm:
            #     consider_dict[(rc1.perm, rc2.perm)].add((rc1, rc2))
    # relation = {}
    # for (perm1, perm2), st in consider_dict.items():
    #     for (rc1, rc2) in st:
    #         for (rc1b, rc2b) in st:
    #             if rc1.lehmer_partial_leq(rc1b) or rc2.lehmer_partial_leq(rc2b):
    #                 relation[(rc1, rc2)] = relation.get((rc1, rc2), set())
    #                 relation[(rc1, rc2)].add((rc1b, rc2b))
    
    # for key in relation:
    #     # print(key, relation[key])
    # ret_elem = []

    ret_elem = []

    for v in consider_dict.values():
        ret_elem = [*ret_elem, *list(v)]

    ret_elem = set(ret_elem)

    bijection = {}
    for rc_graph in sorted(rc_set):
        if rc_graph.perm == perm:
            continue
        old_set = try_lr_module_biject(rc_graph.perm, length)
        for rc1_bad, rc2_bad in sorted(old_set):
            for rc1, rc2 in sorted(consideration_set):
                if (rc1, rc2) not in consider_dict[(rc1_bad.perm, rc2_bad.perm)]:
                    continue
                if (rc1, rc2) not in bijection and (lehmer_partial_pair((rc1_bad, rc2_bad), (rc1, rc2)) or (lehmer_partial_pair((rc1, rc2), (rc1_bad, rc2_bad)))):
                #and (rc1.lehmer_partial_leq(rc1_bad) or rc2.lehmer_partial_leq(rc2_bad)):
                    bijection[(rc1, rc2)] = (rc1_bad, rc2_bad)
                    break
                # meet = next(iter(sorted(spitzu)))
                # bijection[meet] = (rc1_bad, rc2_bad)
                # ret_elem.remove(meet)

    for good_key, bad_key in bijection.items():
        ret_elem.remove(good_key)
    # coll_biject = {}

    # for (rc1_bad, rc2_bad), st in bijection.items():
    #     coll_biject[rc1_bad.perm, rc2_bad.perm] = coll_biject.get((rc1_bad.perm, rc2_bad.perm), set())
    #     coll_biject[rc1_bad.perm, rc2_bad.perm].update(st)

    # iter_biject = {k: iter(v) for k, v in coll_biject.items()}

    # for (rc1, rc2), st in bijection.items():
    #     for rcc in iter_biject[(rc1.perm, rc2.perm)]:
    #         if rcc in ret_elem:
    #             ret_elem.remove(rcc)
    #             break

    return tuple(ret_elem)


def try_lr_module_biject_freeze(perm, length):
    from schubmult import schubmult_py

   
    if perm.inv == 0:
        if length == 0:
            mod = [(RCGraph(), RCGraph())]
        else:
            mod = [(RCGraph([() * length]), RCGraph([() * length]))]
        return mod
    
    rc_set_left = set((fa_elem(perm,length)*RCGraph()).value_dict.keys())
    prin_rc = next(iter({rc for rc in rc_set_left if rc.is_principal}))

    rc_set = set((RCGraph()*fa_elem(perm,length)).value_dict.keys())
    consideration_set = {(k[0], k[1]): v for k, v in ((RCGraph() @ RCGraph())*fa_elem(perm,length).coproduct()).value_dict.items()}


    consideration_set_left = {(k[0], k[1]): v for k, v in ((fa_elem(perm,length).coproduct()*(RCGraph() @ RCGraph()))).value_dict.items()}

    # print(f"{rc_set_left=}")
    # print(f"{rc_set=}")
    first_mod = None
    for rc in rc_set_left:
        if first_mod is None:
            first_mod = 1*rc
        else:
            first_mod += rc
    second_mod = None
    for rc in rc_set:
        if second_mod is None:
            second_mod = 1*rc
        else:
            second_mod += rc
    # print("FIRST")
    # print(first_mod)
    # print("SECOND")
    # print(second_mod)

    perm_set = {rc.perm for rc in rc_set}
    perm_set_left = {rc.perm for rc in rc_set_left}
    #assert perm_set == perm_set_left, f"Perm sets not equal {perm=} {perm_set} {perm_set_left}"

    actual_dict_right = {}
    actual_dict_left = {}

    consider_dict = {}
    consider_dict_left = {}
    for (rc1, rc2), v in consideration_set.items():
        consider_dict[(rc1.perm, rc2.perm)] = consider_dict.get((rc1.perm, rc2.perm), set())
        #consider_dict[(rc1.perm, rc2.perm)].add((rc1.edelman_greene(), rc2.edelman_greene()))
        consider_dict[rc1.perm, rc2.perm].add((rc1, rc2))
        # actual_dict_right[(rc1.perm, rc2.perm)] = actual_dict_right.get((rc1.perm, rc2.perm), set())
        # actual_dict_right[(rc1.perm, rc2.perm)].add((rc1, rc2))

    for (rc1, rc2), v in consideration_set_left.items():
        consider_dict_left[(rc1.perm, rc2.perm)] = consider_dict_left.get((rc1.perm, rc2.perm), set())
        consider_dict_left[(rc1.perm, rc2.perm)].add((rc1.edelman_greene(), rc2.edelman_greene()))
        actual_dict_left[(rc1.perm, rc2.perm)] = actual_dict_left.get((rc1.perm, rc2.perm), set())
        actual_dict_left[(rc1.perm, rc2.perm)].add((rc1, rc2))

    ret_elem = []

    # for key in actual_dict_left:
    #     perm1, perm2 = key
    #     if key in actual_dict_right:
    #         #print("Left set only:")
    #         for rc1, rc2 in sorted(actual_dict_left[key]):
    #             # print(rc1)
    #             # print("-------")
    #             # print(rc2)
    #             # print("=======")
    #             if (rc1, rc2) in actual_dict_right[key]:
    #                 if rc1.perm <= perm and rc2.perm <= perm:
    #                     # print(f"Shared {rc1,rc2}")
    #                     ret_elem.append((rc1.edelman_greene(), rc2.edelman_greene()))
    #         # print("Right set:")
    #         # for rc1, rc2 in sorted(actual_dict_right[key]):
    #         #     # print(rc1)
    #         #     # print("-------")
    #         #     # print(rc2)
    #         #     # print("=======")
    # return ret_elem
    # # ret_elem = []
    # # for key in actual_dict_left:
    # #     if key not in actual_dict_right:
    # #         #print("Left set only:")
    # #         for rc1, rc2 in sorted(actual_dict_left[key]):
    # #             # print(rc1)
    # #             # print("-------")
    # #             # print(rc2)
    # #             # print("=======")
    # #             ret_elem.append((rc1.edelman_greene(), rc2.edelman_greene()))
    # #         # print("Right set:")
    # #         # for rc1, rc2 in sorted(actual_dict_right[key]):
    # #         #     # print(rc1)
    # #         #     # print("-------")
    # #         #     # print(rc2)
    # #         #     # print("=======")

    # return ret_elem

    ret_elem = None

    exclusions = {}
    secret_element = {}
    # rc_set = {rc_graph.perm for rc_graph in rc_set if rc_graph.perm != perm}
    for perm1, perm2 in consider_dict:
        for rc_graph in sorted(rc_set):
            # print(f"{perm.trimcode=}")
            # print(rc_graph)
            if rc_graph.perm != perm:
                        # exclusions[(perm1, perm2)].add((rc1, rc2))
                val = int(schubmult_py({perm1: S.One}, perm2).get(rc_graph.perm, 0))
                # old_set = set(try_lr_module_biject_cache(perm0, lock, shared_cache_dict=shared_cache_dict, length=length))
                exclusions[(perm1, perm2)] = exclusions.get((perm1, perm2), 0) + val
            # else:
            #     # print(f"Principal {perm} yo yo yo for {perm1, perm2}")
            #     for (rc1, rc2) in consider_dict[(perm1, perm2)]:
            #         # print(f"{rc1.lehmer_partial_leq(rc_graph)=} {rc2.lehmer_partial_leq(rc_graph)=}")
            #     val = int(schubmult_py({perm1: S.One}, perm2).get(rc_graph.perm, 0))
            #     # print(f"btw {val=}")
            #     exclusions[(perm1, perm2)].update(old_set)

    # def comp_them(pair1, pair2):
    #     return (pair1[0], pair2[0], pair1[1], pair2[1])

    def lehmer_partial_pair(pair1, pair2):
        return pair1[0].lehmer_partial_leq(pair2[0]) and pair1[1].lehmer_partial_leq(pair2[1])
    minplack_elem = {}
    minindex = {}
    chains = {}
    for (perm1, perm2), st in consider_dict.items():
        lst = sorted(st)
        v = exclusions.get((perm1, perm2), 0)
        try:
            secret_element[(perm1, perm2)] = lst[v]
        except IndexError:
            pass
        chains[(perm1, perm2)] = []
        # print(f"Checking {perm1.trimcode, perm2.trimcode}")
        for j, elem in enumerate(reversed(lst)):
            if len(chains[(perm1, perm2)]) == 0:
                chains[(perm1, perm2)].append([elem])
            else:
                found = False
                for c in chains[(perm1, perm2)]:
                    if lehmer_partial_pair(elem, c[0]):
                        c.insert(0,elem)
                        found = True
                        break
                if not found:
                    chains[(perm1, perm2)].append([elem])
        for c in chains[(perm1, perm2)]:
            if (perm1, perm2) in secret_element and secret_element[(perm1, perm2)] in c:
                assert secret_element[(perm1, perm2)] == c[0]
                # print(f"Good chain: {c=}")
        # print(f"All chains: {chains[(perm1, perm2)]=}")
            # if all(lehmer_partial_pair(elem, other) for other in save_set):
            #     save_set.add(elem)
            #     #last_element = elem
            # else:
            #     try:
            #         minplack_elem[(perm1, perm2)] = lst[j + 1]
            #         minindex[(perm1, perm2)] = j + 1
            #         # print(f"Compare: {v=} {minindex[(perm1, perm2)]=} {len(lst)=} {perm1, perm2=}")
            #         break
            #     except IndexError:
            #         break
        
        # if last_element is not None:
        #     secret_element[(perm1, perm2)] = last_element
    

    ret_elem = []
    for (perm1, perm2), st in consider_dict.items():
        for pair in st:
            if (perm1, perm2) in secret_element:
                if lehmer_partial_pair(secret_element[(perm1, perm2)], pair):
                    # print(f"TAB FOR {pair[0][0].perm, pair[1][0].perm}")
                    # for tab in comp_them(*pair):
                    #     # print(tab)
                    #     # print("-------")
                    ret_elem.append(pair)
                # else:
                    # print(f"Didn't match but is this true? {lehmer_partial_pair(pair, secret_element[(perm1, perm2)])}")
                    # print(f"THIS TAB NO GOOD FOR {pair[0][0].perm, pair[1][0].perm}")
                    # for tab in comp_them(*pair):
                    #     # print(tab)
                    #     # print("-------")
    # prin_rc = next(iter(k for k in (FA(*perm.trimcode,*((0,)*(length-len(perm.trimcode))))*RCGraph()).value_dict.keys() if k.is_principal))
    # st1 = {(rc1, rc2) for (rc1, rc2) in (FA(*prin_rc.length_vector()).coproduct()*(RCGraph()@RCGraph())).value_dict.keys()}
    # st2 = {(rc1.transpose(), rc2.transpose()) for (rc1, rc2) in (FA(*list(reversed(prin_rc.transpose().length_vector()))).coproduct()*(RCGraph()@RCGraph())).value_dict.keys()}
    # ret_elem = []
    # for k in st1:
    #     rc1, rc2 = k
    #     if rc1.perm<=perm and rc2.perm<= perm and k in st2:
    #         ret_elem.append(k)

    return ret_elem


# def nonrecursive_lr_module(perm, length=None):
#     # print(f"Starting {perm}")
#     if length is None:
#         length = len(perm.trimcode)
#     elif length < len(perm.trimcode):
#         raise ValueError("Length too short")
#     if perm.inv == 0:
#         if length == 0:
#             mod = RCGraph() @ RCGraph()
#             #  #  # print(f"MOASA!! {mod=} {type(mod)=}")
#             return mod
#         return FA(*([0] * length)).coproduct() * (RCGraph() @ RCGraph())
#     lower_perm = uncode(perm.trimcode[1:])
#     lower_module1 = try_lr_module(lower_perm, length - 1)
#     assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
#     #  #  # print(f"Coproducting {ASx(uncode([perm.trimcode[0]]), 1).coproduct()=}")
#     #  #  # print(ASx(uncode([perm.trimcode[0]]), 1).coproduct())
#     #  #  # print("Going for it")
#     #  #  # print(f"{type(lower_module1)=} {lower_module1=}")
#     #  #  # print(f"{type(ASx(uncode([perm.trimcode[0]]), 1).coproduct())=}")
#     ret_elem = ASx(uncode([perm.trimcode[0]]), 1).coproduct() * lower_module1
#     assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"

#     ret_elem = TensorModule({RCGraphTensor(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm)})

#     if length == 1:
#         return ret_elem

#     up_elem = ASx(uncode([perm.trimcode[0]]), 1) * elem
#     leftover = {}
#     for key, coeff in up_elem.items():
#         if key[0] != perm:
#             assert coeff == 1
#             for (rc1_bad, rc2_bad), cff2 in lr_module(key[0], length).items():
#                 to_subtract = cff2
#                 if leftover.get((rc1_bad.perm, rc2_bad.perm), 0) > 0:
#                     to_subtract += leftover[(rc1_bad.perm, rc2_bad.perm)]
#                     del leftover[(rc1_bad.perm, rc2_bad.perm)]
#                 keys2 = set(ret_elem.keys())
#                 for rc1, rc2 in keys2:
#                     if rc1.perm == rc1_bad.perm and rc2_bad.perm == rc2.perm:
#                         cff = ret_elem.get((rc1, rc2), 0)

#                         if cff < 0:
#                             raise ValueError(f"Negative {rc1_bad.perm=} {rc2_bad.perm=} {to_subtract=} {cff=}")

#                         assert rc1.length_vector() >= rc1_bad.length_vector() or rc2.length_vector() >= rc2_bad.length_vector()
#                         ret_elem -= TensorModule({RCGraphTensor(rc1, rc2): to_subtract})
#                         break
#     if len(leftover) > 0:
#         raise ValueError(f"Leftover {leftover}")

#     assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {perm.trimcode=}"
#     return ret_elem


@cache
def lr_module(perm, length=None):
    if length is None:
        length = len(perm.trimcode)
    elif length < len(perm.trimcode):
        raise ValueError("Length too short")
    if perm.inv == 0:
        if length == 0:
            return TensorModule({RCGraphTensor(RCGraph(), RCGraph()): 1})
        return FA(*([0] * length)).coproduct() * TensorModule({RCGraphTensor(RCGraph(), RCGraph()): 1})

    lower_perm = uncode(perm.trimcode[1:])
    elem = ASx(lower_perm, length - 1)
    lower_module1 = lr_module(lower_perm, length - 1)
    assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
    ret_elem = ASx(uncode([perm.trimcode[0]]), 1).coproduct() * lower_module1
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"

    ret_elem = TensorModule({RCGraphTensor(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm)})

    if length == 1:
        return ret_elem

    up_elem = ASx(uncode([perm.trimcode[0]]), 1) * elem
    leftover = {}
    for key, coeff in up_elem.items():
        if key[0] != perm:
            assert coeff == 1
            for (rc1_bad, rc2_bad), cff2 in lr_module(key[0], length).items():
                to_subtract = cff2
                if leftover.get((rc1_bad.perm, rc2_bad.perm), 0) > 0:
                    to_subtract += leftover[(rc1_bad.perm, rc2_bad.perm)]
                    del leftover[(rc1_bad.perm, rc2_bad.perm)]
                keys2 = set(ret_elem.keys())
                for rc1, rc2 in keys2:
                    if rc1.perm == rc1_bad.perm and rc2_bad.perm == rc2.perm:
                        cff = ret_elem.get((rc1, rc2), 0)

                        if cff < 0:
                            raise ValueError(f"Negative {rc1_bad.perm=} {rc2_bad.perm=} {to_subtract=} {cff=}")

                        assert rc1.length_vector() >= rc1_bad.length_vector() or rc2.length_vector() >= rc2_bad.length_vector()
                        ret_elem -= TensorModule({RCGraphTensor(rc1, rc2): to_subtract})
                        break
    if len(leftover) > 0:
        raise ValueError(f"Leftover {leftover}")

    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {perm.trimcode=}"
    return ret_elem

@cache
def all_fa_degree(degree, length):
    FA = FreeAlgebra(WordBasis)
    if degree < 0 or length < 0:
        return FA.zero
    if degree == 0 and length == 0:
        return FA.from_dict({(): 1})
    res = FA.zero
    for i in range(degree + 1):
        prev = all_fa_degree(degree - i, length - 1)
        res += FA(i) * prev
    return res


if __name__ == "__main__":
    # test module functionality
    perm1 = Permutation([1, 4, 2, 3])
    graph1 = next(iter(RCGraph.all_rc_graphs(perm1,2)))
    perm2 = Permutation([3, 1, 4, 2])
    graph2 = next(iter(RCGraph.all_rc_graphs(perm2, 3)))

    print(graph1 * graph2)
    #p = 3
    # print(f"Acting with {p} on the right")
    # rc_set = graph1.right_p_act(p)
    

    # for rc in rc_set:
    #     print(rc)
    #     print("-----")
    #     q = 2
    #     print(f"Acting with {q} on the right")
    #     rc_set2 = rc.right_p_act(q)
    #     for rc2 in rc_set2:
    #         print(rc2)
    #         print("=====")