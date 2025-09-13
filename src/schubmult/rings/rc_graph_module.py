import sys
from functools import cache

from symengine import SympifyError

import schubmult.schub_lib.schub_lib as schub_lib
from schubmult.perm_lib import Permutation, uncode
from schubmult.rings.nil_hecke import NilHeckeRing
from schubmult.rings.schubert_ring import DoubleSchubertElement, SingleSchubertRing
from schubmult.symbolic import S, expand, sympify
from schubmult.utils.perm_utils import add_perm_dict, artin_sequences

from .free_algebra import FreeAlgebra, FreeAlgebraElement
from .free_algebra_basis import FreeAlgebraBasis, SchubertBasis, WordBasis
from .nil_hecke import NilHeckeRing
from .tensor_ring import TensorRing, TensorRingElement

FAS = FreeAlgebra(basis=SchubertBasis)


class RCGraph(tuple):
    # def monk_insert(self, row, reps=1):
    #     if reps > 1:
    #         return self.monk_insert(row, reps-1).monk_insert(row)

    #     if len(self) == 0:
    #         tup = [()]* row
    #         tup[row - 1] = (row,)
    #         return RCGraph(tup)
    #     if len(self) < row:
    #         tup = [*self, *([()]*(row - len(self)))]
    #         tup[row - 1] = (row,)
    #         return RCGraph(tup)

    #     if row > 1:
    #         return RCGraph([*self[:row-1],*[tuple([a + row - 1 for a in roww]) for roww in self.rowrange(row - 1, len(self)).monk_insert(1)]])

    #     for j in range(1,100):
    #         if not self.has_element(row, j):
    #             a, b = self.right_root_at(row, i)
    #             if a == 1:
    #                 new_row = tuple(sorted([*self[0], j]))
    #                 #NOTDONE

    # def covers(self, other):
    #     if not isinstance(other, RCGraph):
    #         return NotImplemented
    #     return all(any(a <= b for a in row_a for b in row_b) for row_a, row_b in zip(self, other))

    def __matmul__(self, other):
        return TensorModule({RCGraphTensor(self, other): 1})

    def __rmul__(self, other):
        if isinstance(other, FreeAlgebraElement):
            return other * RCGraphModule({self: 1})
        if isinstance(other, (int, float, S.__class__)):
            return RCGraphModule({self: other})
        raise ValueError("Can't multiply")

    def asdtype(self, cls):
        if cls == FreeAlgebraElement:
            return cls.dtype().ring.from_rc_graph(self)

    def as_nil_hecke(self, x, y=None):
        R = NilHeckeRing(x)
        return self.polyvalue(x, y) * R(self.perm)

    def has_element(self, i, j):
        return i <= len(self) and j + i in self[i - 1]

    def right_root_at(self, i, j):
        # row i column j
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

    def __new__(cls, *args):
        obj = tuple.__new__(cls, *args)
        return obj

    def rowrange(self, start, end):
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
            rcg_row = RCGraph(tuple(new_row))
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
        return tuple(self) <= tuple(other)

    def __lt__(self, other):
        if not isinstance(other, RCGraph):
            return NotImplemented
        return tuple(self) < tuple(other)

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

        # print(tuple(new_rc))
        # print(tuple(self))
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
            # print(f"{pm=} {k=} {new_row=}")
            # nrc = RCGraph([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in self]])
            # assert nrc.perm == perm2
            lst = [tuple([a + 1 for a in row]) for row in self]

            for index in range(max(len(self) + 1, max(new_row))):
                if index < len(lst):
                    if index + 1 in new_row:
                        lst[index] = tuple([*lst[index], index + 1])
                else:
                    if index + 1 in new_row:
                        lst += [(index + 1,)]
                    else:
                        lst += [()]
            nrc = RCGraph(lst)
            # print(nrc)
            assert nrc.perm == ~perm2
            # except AssertionError:
            #     # print(self)
            #     # print(perm2)
            #     # print(nrc.perm)
            #     # print(nrc)
            #     # print(f"{new_row=}")
            #     raise
            ret.add(nrc)
        return ret
        # elem rc

        # if len(self) == 0:
        #     return {}
        # ret = set()
        # pm = self.perm
        # L = schub_lib.pull_out_var(1, pm)
        # perms_to_try = set()
        # for index_list, new_perm in L:
        #     if len(index_list) == p:
        #         new_row = [new_perm[i] for i in range(max(len(pm), len(new_perm))) if new_perm[i] == pm[i + 1]]
        #         new_row.sort(reverse=True)
        #         if tuple(new_row) != self[0]:
        #             continue
        #         perms_to_try.add(new_perm)

        # return ret

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
            # except AssertionError:
            #     # print(self)
            #     # print(perm2)
            #     # print(nrc.perm)
            #     # print(nrc)
            #     # print(f"{new_row=}")
            #     raise
            ret.add(nrc)
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
        return "\n".join(self.as_str_lines())

    def __hash__(self):
        return hash((tuple(self), "RCGRAPH"))


# rc graph module is a module of tuples
# backwards inserted, left action by free algebra


class RCGraphModule(dict):
    def asdtype(self, cls):
        return sum([v * k.asdtype(cls) for k, v in self.items()])

    def __matmul__(self, other):
        return TensorModule.ext_multiply(self, other)

    def __mul__(self, other):
        if isinstance(other, RCGraphModule):
            result = RCGraphModule()
            for rc_graph, coeff in self.items():
                for other_rc_graph, coeff2 in other.items():
                    result += coeff * coeff2 * rc_graph.prod_with_rc(other_rc_graph)
        return result

    def schubvalue(self, sring):
        ret = S.Zero
        for k, v in self.items():
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
        for rc, coeff in self.items():
            res += coeff * rc.coproduct()
        return res

    def apply(self, poly, genset, length):
        from . import ASx

        dct = self._produce_poly_dict(poly, genset, length)

        if len(dct) == 0:
            return RCGraphModule()
        res = RCGraphModule()

        for rc, coeff_rc in self.items():
            fa_elem = ASx(rc.perm, length).change_basis(WordBasis)
            for vec, coeff_vec in fa_elem.items():
                res += coeff_vec * coeff_rc * dct.get(vec, 0) * RCGraphModule({rc: 1})

        # for vec, coeff in dct2.items():
        #     res += RCGraphModule({rc: v * coeff for rc, v in self.items() if rc.length_vector() == vec})
        return res

    def apply_product(self, poly1, poly2, genset, length):
        from . import ASx

        dct1 = self._produce_poly_dict(poly1, genset, length)
        dct2 = self._produce_poly_dict(poly2, genset, length)

        res = RCGraphModule()

        for rc, coeff_rc in self.items():
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

    def __new__(cls, *args, length_limit=-1):
        obj = dict.__new__(cls)
        dct = dict(*args)
        obj.update({k: v for k, v in dct.items() if v != 0})
        obj._length_limit = length_limit
        return obj

    def polyvalue(self, x, y=None):
        ret = S.Zero
        for k, v in self.items():
            ret += v * k.polyvalue(x, y)
        return ret

    def as_nil_hecke(self, x, y=None):
        ret = S.Zero
        for k, v in self.items():
            ret += v * k.as_nil_hecke(x, y)
        return ret

    def __init__(self, *args):
        pass

    def __add__(self, other):
        if isinstance(other, RCGraphModule):
            return RCGraphModule(add_perm_dict(self, other))
        return NotImplemented

    def __neg__(self):
        return RCGraphModule({k: -v for k, v in self.items()})

    def __sub__(self, other):
        if isinstance(other, RCGraphModule):
            return RCGraphModule(add_perm_dict(self, -other))
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, FreeAlgebraElement):
            wd_dict = other.change_basis(WordBasis)
            ret = {}
            for k0, v0 in self.items():
                for k, v in wd_dict.items():
                    addup = RCGraphModule({k0: v0 * v})
                    for a in reversed(k):
                        new_addup = RCGraphModule()
                        for kr, vv in addup.items():
                            new_addup += RCGraphModule(dict.fromkeys(kr.act(a), vv))
                        addup = new_addup

                    ret = add_perm_dict(ret, addup)
            return RCGraphModule(ret)
        if isinstance(other, (int, float, S.__class__)):
            other = sympify(other)
            return RCGraphModule({k: v * other for k, v in self.items()})
        return NotImplemented

    def transpose(self):
        return RCGraphModule({graph.transpose(): v for graph, v in self.items()})

    def as_str_lines(self):
        if len(self.keys()) == 0:
            return ["0"]
        lines = [""]
        first = True
        for k, v in self.items():
            lines2 = k.as_str_lines()
            if len(lines) < len(lines2):
                upstr = ""
                if len(lines[0]) > 0:
                    upstr = " " * len(lines[0])
                lines += [upstr] * (len(lines2) - len(lines))
            padlen = 0
            for i in range(len(lines2)):
                coeffstr = ""
                if not first:
                    if i == 0:
                        coeffstr += " + "
                    else:
                        coeffstr += "   "
                if i == 0:
                    if v == -1:
                        coeffstr += "-"
                    elif v != 1:
                        coeffstr += str(v) + " * "
                    padlen = len(coeffstr)
                else:
                    coeffstr = " " * padlen

                lines[i] += coeffstr + lines2[i]
            first = False
        return lines

    def __str__(self):
        return "\n".join(self.as_str_lines())


class DualRCGraphModule(RCGraphModule):
    def __new__(cls, *args):
        return RCGraphModule.__new__(cls, *args)

    def __init__(self, *args):
        pass

    def __add__(self, other):
        if isinstance(other, DualRCGraphModule):
            return DualRCGraphModule(add_perm_dict(self, other))
        return NotImplemented

    def __neg__(self):
        return DualRCGraphModule({k: -v for k, v in self.items()})

    def __sub__(self, other):
        if isinstance(other, DualRCGraphModule):
            return DualRCGraphModule(add_perm_dict(self, -other))
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, FreeAlgebraElement):
            wd_dict = other.change_basis(WordBasis)
            ret = {}
            for k0, v0 in self.items():
                for k, v in wd_dict.items():
                    addup = DualRCGraphModule({k0: v0 * v})
                    for a in k:
                        new_addup = DualRCGraphModule()
                        for kr, vv in addup.items():
                            new_addup += DualRCGraphModule(dict.fromkeys(kr.dualact(a), vv))
                        addup = new_addup

                    ret = add_perm_dict(ret, addup)
            return DualRCGraphModule(ret)
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        return DualRCGraphModule({k: v * other for k, v in self.items()})

    def __rmul__(self, other):
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        return DualRCGraphModule({k: v * other for k, v in self.items()})

    def as_str_lines(self):
        if len(self.keys()) == 0:
            return ["0"]
        lines = [""]
        first = True
        for k, v in self.items():
            lines2 = k.as_str_lines()
            if len(lines) < len(lines2):
                upstr = ""
                if len(lines[0]) > 0:
                    upstr = " " * len(lines[0])
                lines += [upstr] * (len(lines2) - len(lines))
            padlen = 0
            for i in range(len(lines2)):
                coeffstr = ""
                if not first:
                    if i == 0:
                        coeffstr += " + "
                    else:
                        coeffstr += "   "
                if i == 0:
                    if v == -1:
                        coeffstr += "-"
                    elif v != 1:
                        coeffstr += str(v) + " * "
                    padlen = len(coeffstr)
                else:
                    coeffstr = " " * padlen

                lines[i] += coeffstr + lines2[i] + "^^"
            first = False
        return lines

    def pairing(self, other):
        if not isinstance(other, RCGraphModule):
            return NotImplemented
        ret = S.Zero
        for k1, v1 in self.items():
            for k2, v2 in other.items():
                if k1 == k2:
                    ret += v1 * v2
        return ret

    def __call__(self, other):
        return self.pairing(other)

    def __str__(self):
        return "\n".join(self.as_str_lines())


class RCGraphTensor(tuple):
    def polyvalue(self, x, y=None):
        return self[0].polyvalue(x, y) * self[1].polyvalue(x, y)

    def asdtype(self, cls):
        return cls.dtype().ring.from_rc_graph_tensor(self)

    def __new__(cls, graph1, graph2):
        obj = tuple.__new__(cls, (graph1, graph2))
        return obj

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

    def __str__(self):
        return "\n".join(self.as_str_lines())

    def __hash__(self):
        return hash((tuple(self), "RCGRAPHTENSOR"))


class TensorModule(RCGraphModule):
    # def apply(self, poly1, poly2genset, length):
    #     from . import ASx
    #     from .variables import genset_dict_from_expr
    #     dct = genset_dict_from_expr(poly, genset)
    #     dct2 = {}
    #     for vec, coeff in dct.items():
    #         if len(vec) > length:
    #             return RCGraphModule()
    #         dct2[tuple([0]*(length-len(vec)) + [*vec])] = coeff

    #     res = RCGraphModule()

    #     for (rc1, rc2), coeff_rc in self.items():
    #         fa_elem = ASx(rc.perm, length).change_basis(WordBasis)
    #         for vec, coeff_vec in fa_elem.items():
    #             res += coeff_vec * coeff_rc * dct2.get(vec, 0) * RCGraphModule({rc: 1})

    #     # for vec, coeff in dct2.items():
    #     #     res += RCGraphModule({rc: v * coeff for rc, v in self.items() if rc.length_vector() == vec})
    #     return res

    def apply_product(self, poly1, poly2, genset, length):
        from . import ASx

        res = TensorModule()
        for (rc1, rc2), coeff in self.items():
            res += TensorModule.ext_multiply(RCGraphModule({rc1: coeff}).apply(poly1, genset, length), RCGraphModule({rc2: 1}).apply(poly2, genset, length))

        return res

    def __new__(cls, *args):
        return RCGraphModule.__new__(cls, *args)

    def __add__(self, other):
        if isinstance(other, TensorModule):
            return TensorModule(add_perm_dict(self, other))
        return NotImplemented

    def __neg__(self):
        return TensorModule({k: -v for k, v in self.items()})

    def __sub__(self, other):
        if isinstance(other, TensorModule):
            return TensorModule(add_perm_dict(self, -other))
        return NotImplemented

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

    def __rmul__(self, other):
        if isinstance(other, TensorRingElement):
            ret = {}
            for k0, v0 in self.items():
                for k, v in other.items():
                    addup = TensorModule({k0: v * v0})
                    new_addup = TensorModule()
                    for kr, vv in addup.items():
                        elem1 = other.ring.rings[0](*k[0]) * RCGraphModule({kr[0]: 1})
                        if isinstance(other.ring.rings[1], TensorRing):
                            raise ValueError("Not implemented")
                        elem2 = other.ring.rings[1](*k[1]) * RCGraphModule({kr[1]: 1})

                        new_addup += vv * TensorModule.ext_multiply(elem1, elem2, strict=True)
                    addup = new_addup

                    ret = add_perm_dict(ret, addup)
            return TensorModule(ret)
        try:
            other = sympify(other)
            return TensorModule({k: v * other for k, v in self.items()})
        except Exception:
            return NotImplemented


ASx = FreeAlgebra(SchubertBasis)
FA = FreeAlgebra(WordBasis)


@cache
def try_lr_module(perm, length=None):
    if length is None:
        length = len(perm.trimcode)
    elif length < len(perm.trimcode):
        raise ValueError("Length too short")
    if perm.inv == 0:
        if length == 0:
            return TensorModule({RCGraphTensor(RCGraph(), RCGraph()): 1})
        unit = FA(*([0] * length)).coproduct() * TensorModule({RCGraphTensor(RCGraph(), RCGraph()): 1})
        return unit
    lower_perm = uncode(perm.trimcode[1:])
    elem = ASx(lower_perm, length - 1)
    lower_module1 = try_lr_module(lower_perm, length - 1)
    assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
    ret_elem = ASx(uncode([perm.trimcode[0]]), 1).coproduct() * lower_module1
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"

    ret_elem = TensorModule({RCGraphTensor(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm)})

    if length == 1:
        return ret_elem
    # lower_module2 = FA(0).coproduct() * ret_elem

    # trim_module1 = TensorModule({RCGraphTensor(rc1.rowrange(1, len(rc1)),rc2.rowrange(1,len(rc2))): v for (rc1, rc2), v in ret_elem.items()})
    # trim_module2 = TensorModule({RCGraphTensor(rc1.rowrange(1, len(rc1)),rc2.rowrange(1,len(rc2))): v for (rc1, rc2), v in lower_module2.items() if rc1.perm.bruhat_leq(uncode([0,*perm.trimcode])) and rc2.perm.bruhat_leq(uncode([0,*perm.trimcode]))})

    # print("TRIM1 compare")
    # print(trim_module1)
    # print(lower_module1)

    # ret_elem2 = TensorModule({RCGraphTensor(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm) and rc2.is_principal and (rc1, rc2) not in ret_elem1})

    # ret_elem += TensorModule({RCGraphTensor(rc2, rc1): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm) and rc1.is_principal and not rc2.is_principal})

    up_elem = ASx(uncode([perm.trimcode[0]]), 1) * elem
    leftover = {}
    for key, coeff in up_elem.items():
        if key[0] != perm:
            assert coeff == 1
            for (rc1_bad, rc2_bad), cff2 in try_lr_module(key[0], length).items():
                # to_subtract = cff2
                # if leftover.get((rc1_bad.perm, rc2_bad.perm), 0) > 0:
                #     to_subtract += leftover[(rc1_bad.perm, rc2_bad.perm)]
                #     del leftover[(rc1_bad.perm, rc2_bad.perm)]
                keys2 = set(ret_elem.keys())
                for rc1, rc2 in keys2:
                    if (rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm) and (rc1.length_vector() >= rc1_bad.length_vector() or rc2.length_vector() >= rc2_bad.length_vector()):
                        del ret_elem[RCGraphTensor(rc1, rc2)]
                        break
    # if len(leftover) > 0:
    #     raise ValueError(f"Leftover {leftover}")

    # print("RESULT COMPARE")
    # print(trim_module2)
    # print(ret_module)
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
        unit = FA(*([0] * length)).coproduct() * TensorModule({RCGraphTensor(RCGraph(), RCGraph()): 1})
        return unit
    lower_perm = uncode(perm.trimcode[1:])
    elem = ASx(lower_perm, length - 1)
    lower_module1 = lr_module(lower_perm, length - 1)
    assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
    ret_elem = ASx(uncode([perm.trimcode[0]]), 1).coproduct() * lower_module1
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"

    ret_elem = TensorModule({RCGraphTensor(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm)})

    if length == 1:
        return ret_elem
    # lower_module2 = FA(0).coproduct() * ret_elem

    # trim_module1 = TensorModule({RCGraphTensor(rc1.rowrange(1, len(rc1)),rc2.rowrange(1,len(rc2))): v for (rc1, rc2), v in ret_elem.items()})
    # trim_module2 = TensorModule({RCGraphTensor(rc1.rowrange(1, len(rc1)),rc2.rowrange(1,len(rc2))): v for (rc1, rc2), v in lower_module2.items() if rc1.perm.bruhat_leq(uncode([0,*perm.trimcode])) and rc2.perm.bruhat_leq(uncode([0,*perm.trimcode]))})

    # print("TRIM1 compare")
    # print(trim_module1)
    # print(lower_module1)

    # ret_elem2 = TensorModule({RCGraphTensor(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm) and rc2.is_principal and (rc1, rc2) not in ret_elem1})

    # ret_elem += TensorModule({RCGraphTensor(rc2, rc1): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm) and rc1.is_principal and not rc2.is_principal})

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
                        result = cff - to_subtract
                        if cff < 0:
                            # print("Found leftover")
                            raise ValueError(f"Negative {rc1_bad.perm=} {rc2_bad.perm=} {to_subtract=} {cff=}")
                            # ret_elem -= TensorModule({RCGraphTensor(rc1, rc2): cff})
                            # leftover[(rc1_bad.perm, rc2_bad.perm)] = -result
                        else:
                            # print("Excluding")
                            # print(RCGraphTensor(rc1, rc2))
                            assert rc1.length_vector() >= rc1_bad.length_vector() or rc2.length_vector() >= rc2_bad.length_vector()
                            ret_elem -= TensorModule({RCGraphTensor(rc1, rc2): to_subtract})
                        break
    if len(leftover) > 0:
        raise ValueError(f"Leftover {leftover}")

    # print("RESULT COMPARE")
    # print(trim_module2)
    # print(ret_module)
    assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {perm.trimcode=}"
    return ret_elem


class DualTensorModule(DualRCGraphModule):
    def __new__(cls, *args):
        return DualRCGraphModule.__new__(cls, *args)

    def __add__(self, other):
        if isinstance(other, DualTensorModule):
            return DualTensorModule(add_perm_dict(self, other))
        return NotImplemented

    # @classmethod
    # def ext_multiply(cls, elem1, elem2):
    #     ret = cls()
    #     for key, val in elem1.items():
    #         for key2, val2 in elem2.items():
    #             ret += cls({RCGraphTensor(key, key2): val * val2})
    #     return ret

    def __rmul__(self, other):
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        return DualTensorModule({k: v * other for k, v in self.items()})

    def __mul__(self, other):
        if isinstance(other, TensorRingElement):
            ret = {}
            for k0, v0 in self.items():
                for k, v in other.items():
                    addup = DualTensorModule({k0: v * v0})
                    new_addup = DualTensorModule()
                    for kr, vv in addup.items():
                        elem1 = DualRCGraphModule({kr[0]: 1}) * other.ring.rings[0](*k[0])
                        if isinstance(other.ring.rings[1], TensorRing):
                            elem2 = DualTensorModule({kr[1]: 1}) * other.ring.rings[1](k[1])
                        else:
                            elem2 = DualRCGraphModule({kr[1]: 1}) * other.ring.rings[1](*k[1])

                        new_addup += vv * DualTensorModule(TensorModule.ext_multiply(elem1, elem2))
                    addup = new_addup
                    ret = add_perm_dict(ret, addup)
            return DualTensorModule(ret)
        try:
            other = sympify(other)
            return DualTensorModule({k: v * other for k, v in self.items()})
        except Exception:
            return NotImplemented


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


def schubert_positive_product(poly1, poly2, genset, degree, length, check=False):
    FA = FreeAlgebra(WordBasis)

    result = RCGraphModule()

    one_module = RCGraphModule({RCGraph(): 1})
    for a, cof in all_fa_degree(degree, length).items():
        module = cof * FA(*a) * one_module
        result += module.apply_product(poly1, poly2, genset, length)
        if check:
            assert all(c > 0 for c in result.values())
    return result


def schubert_act(poly, rc_module, genset, degree, length, check=False):
    from .schubert_ring import SingleSchubertRing

    FA = FreeAlgebra(WordBasis)
    ring = SingleSchubertRing(genset)

    result = RCGraphModule()

    one_module = RCGraphModule({RCGraph(): 1})
    for a, cof in all_fa_degree(degree, length).items():
        for graph, coeff in rc_module.items():
            module = cof * FA(*[a1 + b1 for a1, b1 in zip(a, graph.length_vector())]) * one_module
            result += coeff * module.apply_product(poly, ring(graph.perm).expand(), genset, length)
        if check:
            assert all(c > 0 for c in result.values())
    return result


def nilhecke_power(start, end, length, variable):
    from schubmult.abc import x

    ring = NilHeckeRing(x)
    if length > end - start + 1:
        return 0
    if length < 0:
        return 0
    if length == 0:
        return ring.one
    result = ring.zero
    result += variable * ring(Permutation([]).swap(end - 1, end)) * nilhecke_power(start, end - 1, length - 1, variable)
    result += nilhecke_power(start, end - 1, length, variable)
    return result


def change_free_tensor_basis(tensor, old_basis, new_basis):
    new_ring = TensorRing(FreeAlgebra(new_basis), tensor.ring.rings[1])
    new_tensor = new_ring.zero
    original_ring = FreeAlgebra(old_basis)
    for (key1, key2), coeff in tensor.items():
        new_tensor += coeff * new_ring.ext_multiply(original_ring(*key1).change_basis(new_basis), tensor.ring.rings[1](key2))
    return new_tensor


if __name__ == "__main__":
    import sys

    from schubmult import Permutation
    from schubmult.rings import FA, ASx, SchubertBasis, WordBasis

    perms = Permutation.all_permutations(int(sys.argv[1]))

    for perm in perms:
        elem = ASx(perm, len(perm.trimcode)).change_basis(WordBasis)

        for word, coeff in elem.items():
            if uncode(word) != perm and coeff != 0:
                assert len([rc for rc in RCGraph.all_rc_graphs(perm, len(word)) if rc.length_vector() == word]) == 0, f"Failed for {perm} {word} {coeff}"
        for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
            assert rc.is_principal or elem.get(rc.length_vector(), 0) == 0, f"Failed for {perm} {rc.length_vector()} {rc}"
            if not rc.is_principal:
                mod = FA(*rc.length_vector()) * RCGraphModule({rc: 1})
        

        print("Success for", perm)
