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

        self._dict = _value_dict(start_dct, keytype=self.keytype)
        self._ring = ring
        self._generic_key_type = generic_key_type

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
        #  #  # print(f"In rmul {self=} {other=}")
        if not hasattr(other, "ring"):
            other = sympify(other)
            #  #  # print(f"Sympified {other=}")
            return self.clone({k: v * other for k, v in self.value_dict.items()})
        #  #  # print("I try do this")
        mod = self._empty_module()
        #  #  # print(f"{self=}")
        #  #  # print(f"Is this a tensormod? {type(self)=}")
        #  #  # print(f"{type(mod)=}")
        #  #  # print(f"{other=}")
        #  #  # print(f"{self._dict=}")
        for k, v in self._dict.items():
            #  #  # print("ASTIHASOSATMS")
            #  #  # print(type(self.keytype(k)))
            #  #  # print(f"ASOTN{self.keytype(k).ring_act(other)}")
            mod += v * self.clone(self.keytype(k).ring_act(other))
            #  #  # print(f"{mod=}")
            #  #  # print(f"{type(mod)=}")
        return mod

    def __matmul__(self, other):
        return TensorModule(self, other)

    def __str__(self):
        buildups = ""
        if len(self.keys()) == 0:
            return "<zero module>"
        for k, v in self._dict.items():
            st1, mul_joiner, add_joiner = self.coeffify(v)
            to_add = _multiline_join(st1, str(k), joiner=mul_joiner)
            if len(buildups) == 0:
                buildups = to_add
            else:
                buildups = _multiline_join(buildups, to_add, joiner=add_joiner)
        return buildups
    
    def polyvalue(self, x, y=None):
        return sum([v * k.polyvalue(x, y) for k, v in self._dict.items()])

    def coeffify(self, v):
        if v < 0:
            return str(-v), " * ", " + "
        if v == 1:
            return "", "", " + "
        return str(v), " * ", " + "

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

    def weight_word(self):
        perm = self.perm
        nz = len([a for a in perm.trimcode if a != 0])
        root_dict = {perm.right_root_at(index): index for index in range(perm.inv)}
        result_word = [9]*perm.inv
        index = 0
        perm_word = self.perm_word()
        for i, row in enumerate(self):
            for _ in range(len(row)):
                root = perm.right_root_at(index, word=perm_word)
                result_word[root_dict[root]] = (nz - perm.code_index_of_index(index))
                index += 1
        result_word.reverse()
        return tuple(result_word)

    def __matmul__(self, other):
        if isinstance(other, RCGraph):
            return RCGraphModule({self: 1}) @ RCGraphModule({other: 1})
        if isinstance(other, ModuleType):
            return RCGraphModule({self: 1}) @ other
        return NotImplemented

    def __rmul__(self, other):
        return RCGraphModule({self: 1}).__rmul__(other)

    def asdtype(self, cls):
        return cls.dtype().ring.from_rc_graph(self)

    def as_nil_hecke(self, x, y=None):
        R = NilHeckeRing(x)
        return self.polyvalue(x, y) * R(self.perm)

    def has_element(self, i, j):
        return i <= len(self) and j + i in self[i - 1]

    def right_root_at(self, i, j):
        from bisect import bisect_left

        start_root = (i, j + 1)
        if i > len(self):
            return start_root
        row = self[i - 1]
        revved = [*row]
        revved.reverse()

        index = bisect_left(revved, i + j - 1)
        perm = Permutation.ref_product(*revved[:index])
        start_root = (perm[start_root[0] - 1], perm[start_root[1] - 1])
        lower_perm = Permutation([])

        for rrow in self[i:]:
            lower_perm *= ~Permutation.ref_product(*rrow)

        return (lower_perm[start_root[0] - 1], lower_perm[start_root[1] - 1])

    def length_vector(self):
        return tuple([len(row) for row in self])

    # def __new__(cls, *args, **kwargs):
    #     return tuple.__new__(cls, *args)

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

    def __le__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return self.length_vector() <= other.length_vector()

    def __lt__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return self.length_vector() < other.length_vector()

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

    def coproduct(self):
        from . import FA, ASx

        new_set_of_perms = ASx(self.perm, len(self)).coproduct()
        rc_set = FA(*self.length_vector()).coproduct() * TensorModule({RCGraphTensor(RCGraph(), RCGraph()): 1})
        result = TensorModule()
        for (rc1, rc2), coeff in rc_set.items():
            if ((rc1.perm, len(self)), (rc2.perm, len(self))) in new_set_of_perms:
                result += TensorModule({RCGraphTensor(rc1, rc2): new_set_of_perms[((rc1.perm, len(self)), (rc2.perm, len(self)))]})
        return result

    def prod_with_rc(self, other):
        from . import FA, ASx

        new_set_of_perms = ASx(self.perm, len(self)) * ASx(other.perm, len(other))
        rc_set = set((FA(*self.length_vector()) * RCGraphModule({other: 1})).keys())
        result = RCGraphModule()
        for (perm, length), coeff in new_set_of_perms.items():
            result += RCGraphModule({rc: coeff for rc in rc_set if rc.perm == perm})
        return result

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
            for i in range(last_desc+1, mx + 1):
                if old_perm[i-1] < old_perm[i]:
                    new_perm = old_perm.swap(i-1, i)
                    if max((~new_perm).descents()) + 1 > len(rc):
                        continue
                    new_top_row = [i, *rc[0]]
                    new_rc = RCGraph([tuple(new_top_row), *rc[1:]])
                    ret.add(new_rc)
        return ret

        # pm = self.perm
        # elem = FAS(pm, len(self))
        # bumpup = FAS(uncode([p]), 1) * elem
        # ret = set()
        # for k, v in bumpup.items():
        #     perm2 = k[0]
        #     new_row = [pm[i] for i in range(max(len(pm), len(perm2))) if pm[i] == perm2[i + 1]]
        #     new_row.sort(reverse=True)
        #     nrc = RCGraph([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in self]])
        #     assert nrc.perm == perm2
        #     ret.add(nrc)
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

    def __leq__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        if len(self) != len(other):
            return NotImplemented
        for i in range(len(self)):
            perm1 = Permutation.ref_product(*self[i])
            perm2 = Permutation.ref_product(*other[i])
            if not perm1.bruhat_leq(perm2):
                return False
        return True

    @property
    def is_principal(self):
        return self.perm == uncode(self.length_vector())

    def __str__(self):
        lines = self.as_str_lines()
        lines2 = [line for line in lines]
        return "\n".join(lines2)

    def __repr__(self):
        return "RCGraph(" + ", ".join([repr(k) for k in self]) + ")"

    def __hash__(self):
        return hash(tuple(self))


class RCGraphModule(ModuleType):
    def clone(self, *args):
        if len(args) == 0:
            return RCGraphModule(self._dict, generic_key_type=self._generic_key_type, ring=self._rings)
        return RCGraphModule(*args, generic_key_type=self._generic_key_type, ring=self._ring)

    def keytype(self, k, value=1):
        return RCGraph(k, value=value)

    def __init__(self, *args, generic_key_type=RCGraph, ring=FreeAlgebra(WordBasis), **kwargs):
        super().__init__(*args, generic_key_type=generic_key_type, ring=ring, **kwargs)

    # def __mul__(self, other):
    #     if isinstance(other, RCGraphModule):
    #         result = RCGraphModule()
    #         for rc_graph, coeff in self._dict.items():
    #             for other_rc_graph, coeff2 in other.items():
    #                 result += coeff * coeff2 * rc_graph.prod_with_rc(other_rc_graph)
    #     return result

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
        res = RCGraphModule()
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
                        new_keys.append(key)
                else:
                    new_keys.append(RCGraph(key))
            self._keys = tuple(new_keys)
        else:
            new_keys = []
            for key in keys:
                if isinstance(key, RCGraphTensor):
                    new_keys += [key._modules[i].keytype(key[i]) for i, k in enumerate(key)]
                else:
                    new_keys.append(RCGraph(key))
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
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {lower_perm=} {length=}"

    ret_elem = ret_elem.clone({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)})

    if length == 1:
        return ret_elem
    keys = set(ret_elem.keys())
    # print(f"{repr(keys)=} {perm=}")
    up_elem = ASx(uncode([perm.trimcode[0]]), 1) * elem
    # print(f"{up_elem=}")
    for key, coeff in up_elem.items():
        if key[0] != perm:
            assert coeff == 1, f"failed coeff 1 {coeff=}"
            # print(f"Iteration {key[0]}")
            for (rc1_bad, rc2_bad), cff2 in try_lr_module(key[0], length).items():
                keys2 = set(keys)
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


def co_rc(rc11, rc222):
    from schubmult.abc import y
    from schubmult.rings.variables import ZeroGeneratingSet
    from schubmult.symbolic import S, expand
    rc1, rc2 = rc11
    rc12, rc22 = rc222
    if len(rc1) != len(rc12):
        return False
    if len(rc2) != len(rc22):
        return False
    if rc1.perm != rc12.perm:
        return False
    if rc2.perm != rc22.perm:
        return False
    return True
    # if len(rc1) <= 1:
    #     return True
    # left_vector = [a + b for a, b in zip(rc1.transpose().length_vector(), rc2.transpose().length_vector())]
    # left_vector.reverse()

    # stanky = FA(*left_vector).coproduct() * (RCGraph() @ RCGraph())
    # for rc_bob_1, rc_bob_2 in stanky.keys():
    #     if rc_bob_1.transpose() == rc12 and rc_bob_2.transpose() == rc22:
    #         return True
    # # da_vec1 = list(reversed(rc1.transpose().length_vector()))
    # # da_vec11 = list(reversed(rc12.transpose().length_vector()))
    # # da_vec2 = list(reversed(rc2.transpose().length_vector()))
    # # da_vec22 = list(reversed(rc22.transpose().length_vector()))

    # # da_da_vec = [a+b for a,b in zip(da_vec1, da_vec2)]
    # # da_da_vec2 = [a+b for a,b in zip(da_vec11,da_vec22)]
    # # mod1 = FA(*da_da_vec).coproduct() * (RCGraph()@RCGraph())
    # # for rc0, rc00 in mod1.keys():
    # #     if rc0.transpose() == rc12 and rc00.transpose() == rc22:
    # #         return True
    # # left is the prin
    # # bottom_left = rc1.rowrange(1,len(rc1))
    # # bottom_right = rc2.rowrange(1, len(rc2))
    # # for i in range(len(rc1)-1):
    # #     row1 = bottom_left[i]
    # #     row2 = bottom_right[i]
    # #     if len(row1) > len(row2):
    # #         return False
    # #     rev_row1 = tuple(reversed(row1))
    # #     rev_row2 = tuple(reversed(row2))
    # #     for i in range(len(rev_row1)):
    # #         if rev_row2[i] != rev_row1[i]:
    # #             return False
    # # print(bottom_left)
    # # print("Matches")
    # # print(bottom_right)
    # return False


def rank(rc1, rc2, n):
    """
    Rank of the RC graphs by lex weight/order
    """
    #from schubmult.utils.perm_utils import artin_sequences
    length = len(rc1)
    assert len(rc1) == len(rc2)
    deg = rc1.perm.inv + rc2.perm.inv
    seqs = all_fa_degree(deg, length)
    weight_vec = tuple([a+b for a,b in zip(rc1.length_vector(),rc2.length_vector())])
    seqs = tuple(sorted([seq for seq in seqs]))
    rank = 0
    for seq in seqs:
        if len(uncode(seq)) > n:
            continue
        if seq != weight_vec:
            rank += len([(rc[0], rc[1]) for rc in (FA(*seq).coproduct()* (RCGraph()@RCGraph())).value_dict.keys() if rc[0].perm == rc1.perm and rc[1].perm == rc2.perm and len(rc1.perm)<=n and len(rc2.perm)<=n])
        else:
            rank  += len([(rc[0], rc[1]) for rc in (FA(*seq).coproduct()* (RCGraph()@RCGraph())).value_dict.keys() if rc[0].perm == rc1.perm and rc[1].perm == rc2.perm and (rc[0],rc[1]) <= (rc1, rc2) ])
            break
    return rank


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
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {lower_perm=} {length=}"

    ret_elem = ret_elem.clone({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)})

    if length == 1:
        return ret_elem
    keys = set(ret_elem.keys())
    # print(f"{repr(keys)=} {perm=}")
    up_elem = ASx(uncode([perm.trimcode[0]]), 1) * elem
    # print(f"{up_elem=}")
    for key in sorted(up_elem.keys(), key=lambda k: (k[0].trimcode)):
        coeff = up_elem[key]
        if key[0] != perm:
            assert coeff == 1, f"failed coeff 1 {coeff=}"
            key_module = try_lr_module_inject(key[0], length=length)
            lst = list(sorted(key_module.value_dict.keys()))
            perm_count_sorted = {}
            perm_count_sorted2 = {}
            rank = {}
            rank2 = {}
            for (rc1, rc2) in lst:
                rank[(rc1, rc2)] = perm_count_sorted.get((rc1.perm, rc2.perm), 0)
                perm_count_sorted[(rc1.perm, rc2.perm)] = perm_count_sorted.get((rc1.perm, rc2.perm), 0) + 1
            for (rc1, rc2) in sorted(keys):
                rank2[(rc1, rc2)] = perm_count_sorted2.get((rc1.perm, rc2.perm), 0)
                perm_count_sorted2[(rc1.perm, rc2.perm)] = perm_count_sorted2.get((rc1.perm, rc2.perm), 0) + 1
            keys2 = list(sorted(keys))
            for (rc1, rc2) in lst:
                for (rc1_check, rc2_check) in sorted(keys2):
                    if co_rc((rc1,rc2),(rc1_check,rc2_check)) and rank[(rc1,rc2)] == rank2[(rc1_check,rc2_check)]:
                        keys.remove((rc1_check,rc2_check))
    ret_elem = ret_elem.clone({k: v for k, v in ret_elem.items() if k in keys})
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {perm.trimcode=}"
    
    return ret_elem

@cache
def try_lr_module_biject(perm):
    # print(f"Starting {perm}")
    if perm.inv == 0:
        mod = RCGraph() @ RCGraph()
        return mod

    rc_set = FA(*perm.trimcode)*RCGraph()
    consideration_set = {(k[0],k[1]) for k in (FA(*perm.trimcode).coproduct() * (RCGraph() @RCGraph())).value_dict.keys() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)}

    consideration_list = list(sorted(consideration_set))

    ret_elem = None

    for rc_graph in rc_set.value_dict.keys():
        if rc_graph.perm != perm:
            rcs = try_lr_module_biject(rc_graph.perm)
            for rc in rcs.keys():
                for (rc1, rc2) in consideration_list:
                    if rc1.perm == rc[0].perm and rc2.perm == rc[1].perm and (rc1, rc2) in consideration_set:
                        consideration_set.remove((rc1, rc2))
                        break

        # if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm) and rank(rc1, rc2) > lr_rank[(rc1.perm, rc2.perm)]:
        #     if ret_elem is None:
        #         ret_elem = rc1 @ rc2
        #     else:
        #         ret_elem += rc1 @ rc2
    for (rc1, rc2) in consideration_set:
        if ret_elem is None:
            ret_elem = rc1 @ rc2
        else:
            ret_elem += rc1 @ rc2
    return ret_elem



def try_lr_module_biject_cache(perm, lock, shared_cache_dict, length):
    # print(f"Starting {perm}")
    from schubmult import schubmult_py
    ret_elem = None
    with lock:
        if perm in shared_cache_dict:
            ret_elem = shared_cache_dict[perm]

    if ret_elem is not None:
        return ret_elem

    if perm.inv == 0:
        if length == 0:
            mod = [(RCGraph(),RCGraph())]
        else:
            mod = [(RCGraph(() * length),RCGraph(() * length))]
        with lock:
            shared_cache_dict[perm] = mod
        return mod
    rc_set = {rc for rc in (FA(*perm.trimcode, *((0,)*(length-len(perm.trimcode))))*RCGraph()).value_dict.keys()}
    consideration_set = {(k[0],k[1]): v for k, v in (FA(*perm.trimcode, *((0,)*(length-len(perm.trimcode)))).coproduct() * (RCGraph() @RCGraph())).value_dict.items()}

    consider_dict = {}
    for (rc1, rc2), v in consideration_set.items():
        consider_dict[(rc1.perm, rc2.perm)] = consider_dict.get((rc1.perm, rc2.perm), set())
        consider_dict[(rc1.perm, rc2.perm)].add((rc1, rc2))
    #consideration_list = {rc1.permlist(sorted(consideration_set, key=lambda rc: (rc[0].perm.trimcode, rc[1].perm.trimcode)))

    ret_elem = None

    for (perm1, perm2) in consider_dict:
        for rc_graph in sorted(rc_set):
            if rc_graph.perm != perm:
                val = int(schubmult_py({perm1: S.One}, perm2).get(rc_graph.perm, 0))
                lst = list(sorted(consider_dict[(perm1, perm2)]))
                for i in range(val):
                    consider_dict[(perm1, perm2)].remove(lst[i])


        # if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm) and rank(rc1, rc2) > lr_rank[(rc1.perm, rc2.perm)]:
        #     if ret_elem is None:
        #         ret_elem = rc1 @ rc2
        #     else:
        #         ret_elem += rc1 @ rc2
    for k, v in consider_dict.items():
        #if len(v)!= schubmult_py({k[0]: S.One}, k[1]).get(perm, 0):
            # print("OH NO")

        if ret_elem is None:
            ret_elem = list(v)
        else:
            ret_elem += list(v)
    # print(consider_dict)
    with lock:
        shared_cache_dict[perm] = ret_elem
    return ret_elem


def nonrecursive_lr_module(perm, length=None):
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
    lower_module1 = try_lr_module(lower_perm, length - 1)
    assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
    #  #  # print(f"Coproducting {ASx(uncode([perm.trimcode[0]]), 1).coproduct()=}")
    #  #  # print(ASx(uncode([perm.trimcode[0]]), 1).coproduct())
    #  #  # print("Going for it")
    #  #  # print(f"{type(lower_module1)=} {lower_module1=}")
    #  #  # print(f"{type(ASx(uncode([perm.trimcode[0]]), 1).coproduct())=}")
    ret_elem = ASx(uncode([perm.trimcode[0]]), 1).coproduct() * lower_module1
    #  #  # print(f"{ret_elem=}")
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {lower_perm=} {length=}"

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
        for (rc1_bad, rc2_bad), cff2 in (FA(*code_key).coproduct()*(RCGraph() @ RCGraph())).items():
            keys2 = set(keys)
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


# class DualTensorModule(DualRCGraphModule):
#     def __new__(cls, *args):
#         return DualRCGraphModule.__new__(cls, *args)

#     def __add__(self, other):
#         if isinstance(other, DualTensorModule):
#             return DualTensorModule(add_perm_dict(self, other))
#         return NotImplemented

#     def __rmul__(self, other):
#         try:
#             other = sympify(other)
#         except SympifyError:
#             return NotImplemented
#         return DualTensorModule({k: v * other for k, v in self._dict.items()})

#     def __mul__(self, other):
#         if isinstance(other, TensorRingElement):
#             ret = {}
#             for k0, v0 in self._dict.items():
#                 for k, v in other.items():
#                     addup = DualTensorModule({k0: v * v0})
#                     new_addup = DualTensorModule()
#                     for kr, vv in addup.items():
#                         elem1 = DualRCGraphModule({kr[0]: 1}) * other.ring.rings[0](*k[0])
#                         if isinstance(other.ring.rings[1], TensorRing):
#                             elem2 = DualTensorModule({kr[1]: 1}) * other.ring.rings[1](k[1])
#                         else:
#                             elem2 = DualRCGraphModule({kr[1]: 1}) * other.ring.rings[1](*k[1])

#                         new_addup += vv * DualTensorModule(TensorModule.ext_multiply(elem1, elem2))
#                     addup = new_addup
#                     ret = add_perm_dict(ret, addup)
#             return DualTensorModule(ret)
#         try:
#             other = sympify(other)
#             return DualTensorModule({k: v * other for k, v in self._dict.items()})
#         except Exception:
#             return NotImplemented


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
    perm1 = Permutation([4, 1, 2, 5, 3])
    graph1 = next(iter(RCGraph.all_rc_graphs(perm1)))
    perm2 = Permutation([3, 1, 2])
    graph2 = next(iter(RCGraph.all_rc_graphs(perm2)))

    mod1 = FA(2) * graph1
    mod2 = FA(3) * graph2

    #  #  # print(mod1)
    #  #  # print(mod2)

    tmod = mod1 @ mod2
