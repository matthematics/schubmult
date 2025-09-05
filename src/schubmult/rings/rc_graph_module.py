from functools import cache

from scipy.fftpack import dct
from symengine import SympifyError

import schubmult.schub_lib.schub_lib as schub_lib
from schubmult.perm_lib import Permutation, uncode
from schubmult.rings.nil_hecke import NilHeckeRing
from schubmult.rings.schubert_ring import DoubleSchubertElement, SingleSchubertRing
from schubmult.symbolic import S, expand, sympify
from schubmult.utils.perm_utils import add_perm_dict

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

    def as_nil_hecke(self, x, y=None):
        R = NilHeckeRing(x)
        return self.polyvalue(x, y) * R(self.perm)


    def has_element(self, i, j):
        return i <= len(self) and j + i in self[i-1]

    def right_root_at(self, i, j):
        # row i column j
        from bisect import bisect_left
        start_root = (i,j+1)
        if i > len(self):
            return start_root
        row = self[i-1]
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
                ret *= x[i+1]**len(row)
            else:
                ret *= prod([x[i+1] - y[row[j] - i] for j in range(len(row))])
        return ret

    @classmethod
    @cache
    def all_rc_graphs(cls, perm, length=-1):
        if perm.inv == 0:
            return {RCGraph(())}
        if len(perm.trimcode) == 1:
            return {RCGraph((tuple(range(perm.code[0], 0, -1)),))}
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
                    nrc = RCGraph((*nrc,*tuple([()]*(length - len(nrc)))))
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
        while len(newrc) < len(self) or any(len(row)>0 for row in trimself):
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
        bumpup = FAS(uncode([p]), 1)* elem
        ret = set()
        for k, v in bumpup.items():
            perm2 = k[0]
            new_row = [pm[i] for i in range(max(len(pm), len(perm2))) if pm[i] == perm2[i+1]]
            new_row.sort()
            # print(f"{pm=} {k=} {new_row=}")
            # nrc = RCGraph([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in self]])
            # assert nrc.perm == perm2
            lst = [tuple([a + 1 for a in row]) for row in self]

            for index in range(max(len(self)+1, max(new_row))):
                if index < len(lst):
                    if index + 1 in new_row:
                        lst[index] = tuple([*lst[index],index+1])
                else:
                    if index + 1 in new_row:
                        lst += [(index + 1,)]
                    else:
                        lst += [()]
            nrc = RCGraph(lst)
            # print(nrc)
            assert nrc.perm == ~perm2
            # except AssertionError:
            #     print(self)
            #     print(perm2)
            #     print(nrc.perm)
            #     print(nrc)
            #     print(f"{new_row=}")
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
        new_set_of_perms = ASx(self.perm,len(self)).coproduct()
        rc_set = (FA(*self.length_vector()).coproduct()*TensorModule({RCGraphTensor(RCGraph(),RCGraph()): 1}))
        result = TensorModule()
        for (rc1, rc2), coeff in rc_set.items():
            if ((rc1.perm, len(self)), (rc2.perm, len(self))) in new_set_of_perms:
                result += TensorModule({RCGraphTensor(rc1, rc2): new_set_of_perms[((rc1.perm, len(self)), (rc2.perm, len(self)))]})
        return result

    def prod_with_rc(self, other):
        from . import FA, ASx
        new_set_of_perms = ASx(self.perm,len(self)) * ASx(other.perm, len(other))
        rc_set = set((FA(*self.length_vector())*RCGraphModule({other: 1})).keys())
        result = RCGraphModule()
        for (perm, length), coeff in new_set_of_perms.items():
            result += RCGraphModule({rc: coeff for rc in rc_set if rc.perm == perm})
        return result

    def act(self, p):
        pm = self.perm
        elem = FAS(pm, len(self))
        bumpup = FAS(uncode([p]), 1)* elem
        ret = set()
        for k, v in bumpup.items():
            perm2 = k[0]
            new_row = [pm[i] for i in range(max(len(pm), len(perm2))) if pm[i] == perm2[i+1]]
            new_row.sort(reverse=True)
            nrc = RCGraph([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in self]])
            assert nrc.perm == perm2
            # except AssertionError:
            #     print(self)
            #     print(perm2)
            #     print(nrc.perm)
            #     print(nrc)
            #     print(f"{new_row=}")
            #     raise
            ret.add(nrc)
        return ret

    def as_str_lines(self):
        lines = []
        for i in range(len(self)):
            row = self[i]
            splug_row = [*row, i]
            row = [str(splug_row[j]) + ("  ")*(splug_row[j] - splug_row[j+1] - 1) for j in range(len(splug_row)-1)]
            line = ""
            line += " ".join([str(r) for r in row])
            lines += [line]
        if len(lines) == 0:
            lines = ["" ]
        lines2 = []
        ml = max([len(line) for line in lines])
        if len(lines[0]) < ml:
            lines[0] = " " * (ml - len(lines[0])) + lines[0]
        for line in lines:
            lines2 += [" "*(ml - len(line)) + line]
        return lines2

    def __str__(self):
        return "\n".join(self.as_str_lines())

    def __hash__(self):
        return hash((tuple(self), "RCGRAPH"))




# rc graph module is a module of tuples
# backwards inserted, left action by free algebra


class RCGraphModule(dict):

    

    def __mul__(self, other):
        if isinstance(other, RCGraphModule):
            result = RCGraphModule()
            for rc_graph, coeff in self.items():
                for other_rc_graph, coeff2 in other.items():
                    result += coeff * coeff2 * rc_graph.prod_with_rc(other_rc_graph)
        return result

    def _produce_poly_dict(self, poly, genset, length):
        from .variables import genset_dict_from_expr
        dct = genset_dict_from_expr(poly, genset)
        dct2 = {}
        for vec, coeff in dct.items():
            if len(vec) > length:
                return {}
            dct2[tuple([*vec]+ [0]*(length-len(vec)))] = coeff
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

    def __new__(cls, *args):
        obj = dict.__new__(cls)
        dct = dict(*args)
        obj.update({k: v for k, v in dct.items() if v != 0})
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
                    addup = RCGraphModule({k0: v0*v})
                    for a in reversed(k):
                        new_addup = RCGraphModule()
                        for kr, vv in addup.items():
                            new_addup += RCGraphModule(dict.fromkeys(kr.act(a), vv))
                        addup = new_addup

                    ret = add_perm_dict(ret, addup)
            return RCGraphModule(ret)
        try:
            other = sympify(other)
        except SympifyError:
            return NotImplemented
        return RCGraphModule({k: v * other for k, v in self.items()})

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
                    addup = DualRCGraphModule({k0: v0*v})
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

    def polyvalue(self, x):
        return self[0].polyvalue(x) * self[1].polyvalue(x)


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

    @classmethod
    def ext_multiply(cls, elem1, elem2):
        ret = cls()
        for key, val in elem1.items():
            for key2, val2 in elem2.items():
                ret += cls({RCGraphTensor(key, key2): val * val2})
        return ret

    def __rmul__(self, other):
        if isinstance(other, TensorRingElement):
            ret = {}
            for k0, v0 in self.items():
                for k, v in other.items():
                    addup = TensorModule({k0: v*v0})
                    new_addup = TensorModule()
                    for kr, vv in addup.items():
                        elem1 = other.ring.rings[0](*k[0]) * RCGraphModule({kr[0]: 1})
                        if isinstance(other.ring.rings[1], TensorRing):
                            elem2 = other.ring.rings[1](k[1]) * TensorModule({kr[1]: 1})
                        else:
                            try:
                                elem2 = other.ring.rings[1](*k[1]) * RCGraphModule({kr[1]: 1})
                            except Exception:
                                elem2 = other.ring.rings[1](k[1]) * RCGraphModule({kr[1]: 1})

                        new_addup += vv * TensorModule.ext_multiply(elem1, elem2)
                    addup = new_addup

                    ret = add_perm_dict(ret, addup)
            return TensorModule(ret)
        try:
            other = sympify(other)
            return TensorModule({k: v * other for k, v in self.items()})
        except Exception:
            return NotImplemented


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
                    addup = DualTensorModule({k0: v*v0})
                    new_addup = DualTensorModule()
                    for kr, vv in addup.items():
                        elem1 = DualRCGraphModule({kr[0]: 1}) *other.ring.rings[0](*k[0])
                        if isinstance(other.ring.rings[1], TensorRing):
                            elem2 = DualTensorModule({kr[1]: 1}) * other.ring.rings[1](k[1])
                        else:
                            elem2 =  DualRCGraphModule({kr[1]: 1}) * other.ring.rings[1](*k[1])

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
    for i in range(degree+1):
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
    result += variable * ring(Permutation([]).swap(end-1,end)) * nilhecke_power(start, end - 1, length - 1, variable)
    result += nilhecke_power(start, end-1, length, variable)
    return result
         



def change_free_tensor_basis(tensor, old_basis, new_basis):
    new_ring = TensorRing(FreeAlgebra(new_basis),tensor.ring.rings[1])
    new_tensor = new_ring.zero
    original_ring = FreeAlgebra(old_basis)
    for (key1, key2), coeff in tensor.items():
        new_tensor += coeff * new_ring.ext_multiply(original_ring(*key1).change_basis(new_basis), tensor.ring.rings[1](key2))
    return new_tensor

def artin_sequences(n):
    if n == 0:
        return set([()])
    old_seqs = artin_sequences(n-1)

    ret = set()
    for seq in old_seqs:
        for i in range(n+1):
            ret.add((i,*seq))
    return ret

def sparkly_sequences(n, L):
    if L == 0:
        return set([()])
    old_seqs = sparkly_sequences(n, L-1)

    ret = set()
    for seq in old_seqs:
        for i in range(n+1):
            ret.add((i,*seq))
    return ret


if __name__ == "__main__":
    from schubmult.abc import E, x, y, z
    from schubmult.rings import FA, ASx, DSx, ElementaryBasis, MonomialBasis, PolynomialAlgebra, Sx
    from schubmult.symbolic import prod

    from .variables import genset_dict_from_expr

    n = 4
    perms = Permutation.all_permutations(n)


    seqs = artin_sequences(n-1)
    seqs2 = sparkly_sequences(n-1,n-1)
    #yring = SingleSchubertRing(y)
    yring = Sx([]).ring
    zring = SingleSchubertRing(z)
    #ring = TensorRing(FreeAlgebra(SchubertBasis)@FreeAlgebra(SchubertBasis),NilHeckeRing(x))
    #ring = TensorRing(FreeAlgebra(SchubertBasis),NilHeckeRing(x))
    #ring = TensorRing(yring, TensorRing(FreeAlgebra(SchubertBasis)@FreeAlgebra(SchubertBasis), NilHeckeRing(x)))
    #ring = TensorRing(PolynomialAlgebra(MonomialBasis(x, n-1)), TensorRing(FreeAlgebra(SchubertBasis), NilHeckeRing(x)))
    ring = TensorRing(FreeAlgebra(SchubertBasis)@FreeAlgebra(SchubertBasis), NilHeckeRing(x))
    ring2 = TensorRing(Sx([]).ring, ring)
    #ring = TensorRing(TensorRing(zring,FreeAlgebra(SchubertBasis)),TensorRing(yring, NilHeckeRing(x)))
    result = ring2.zero
    
    dp = -1
    def descent_pickle(perm, desc):
        if desc < 0:
            return True
        return (perm.inv == 0 or max(perm.descents()) == desc - 1)

    # for perm in perms:
    #     #nil_elem = ring.rings[1].rings[1](perm)
    #     nil_elem=ring.rings[1].rings[1].one
    #     free_elem = ASx(perm,n-1).change_basis(WordBasis)
    #     #result += ring((perm, ((perm,perm))))
    #     result += ring.ext_multiply(ring.rings[0](yring(perm).expand()),ring.rings[1].ext_multiply(free_elem, nil_elem))
    # plus
    # print(result)
    # exit()
    mod1 = RCGraphModule({RCGraph(): 1})
    result0 = ring2.zero
    for seq in seqs:
        modmod = (FA(*seq)*mod1).as_nil_hecke(x)
        print(f"{seq=}")
        print(f"{modmod=}")
        ding = ring.zero
        for kingo, valval in modmod.items():
            ding += ring.ext_multiply(FreeAlgebraBasis.change_tensor_basis(FA(*((0,)*(n-1))).coproduct(),SchubertBasis,SchubertBasis), ring.rings[1](kingo))
                              
        result += ring2.ext_multiply(Sx(prod([x[i+1]**seq[i] for i in range(len(seq))])),ding + ring.ext_multiply(FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(),SchubertBasis,SchubertBasis), ring.rings[1].one))

        #ring.ext_multiply(FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(),SchubertBasis,SchubertBasis), ring.rings[1](Permutation([]))))
    # result2 = ring2.zero
    # for key, coeff in result0.items():
    #     result += ring2.ext_multiply(Sx(coeff), ring(key))

    #result2 = result2.ring.from_dict({k: v for k, v in result2.items() if len(k[0])<=n and len(k[1][0][0])<=n and len(k[1][1]) <= n})

    print(result)
    Permutation.print_as_code=True
    separate = {}
    for key, value in result.items():
        # if any(len(permperm) > n for permperm in (Sx(key[1][0][0][0])*Sx(key[1][0][1][0])).keys()):
        #     continue
        #assert key[0] == key[1][0][0] or key[0] == key[1][1] or key[1][0][0] == key[1][1] or value == 0, f"{key=} {value=}"
        if key[1][1].inv == 0:
            #separate[key[0]] = separate.get(key[0],ring2.rings[1].rings[0].zero) + value*ring2.rings[1].rings[0](key[1][0])
            separate[key[1][0]] = separate.get(key[1][0],ring2.rings[0].zero) + value*ring2.rings[0](key[0])
        # if key[0].inv == (n*(n-1))//2:
        #     separate[key[1][0]] = separate.get(key[1][0],ring2.rings[1].rings[1].zero) + value*ring2.rings[1].rings[1](key[1][1])

    for perm in separate:
        if any(len(permperm) > n for permperm, val in (Sx(perm[0][0])*Sx(perm[1][0])).items() if val != S.Zero):
            continue
        #print(f"{(perm[0][0].trimcode,perm[1][0].trimcode)}: {separate[perm]}")
        print(f"Test {(perm[0][0].trimcode,perm[1][0].trimcode)}: {separate[perm]}")
        testval = (separate[perm] - Sx(perm[0][0])*Sx(perm[1][0]))
        if testval.expand() != S.Zero:
            print(f"Failures for {(perm[0][0].trimcode,perm[1][0].trimcode)}: {[(perm2, key) for perm2, key in testval.items() if key != 0]}")
    
    exit()


    if False:
        nil_elem = ring.rings[1].rings[1].one
        
            #nil_elem *= nilhecke_power(index + 1, n - 1, value, 1)
        interim_result = ring.rings[1].zero

        #for perm, coeff in nil_elem.items():
        #interim_result += ring.rings[1].ext_multiply(yring(prod([y[i+1]**seq[i] for i in range(len(seq))])),nil_elem)
        #interim_result += ring.rings[1].ext_multiply(FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(),SchubertBasis,SchubertBasis),nil_elem)
        interim_result += ring.rings[1].ext_multiply(FA(*seq).change_basis(SchubertBasis),nil_elem)
        #interim_result += ring.rings[1].ext_multiply(yring.one,nil_elem)

        #result += ring.ext_multiply(FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(),SchubertBasis,SchubertBasis), nil_elem)
        result += ring.ext_multiply(ring.rings[0](prod([x[i+1]**seq[i] for i in range(len(seq))])), interim_result)
        #result += ring.ext_multiply(ring.rings[0].ext_multiply(zring(prod([z[i+1]**seq[i] for i in range(len(seq))])),FA(*seq).change_basis(SchubertBasis)), interim_result)
        #result += ring.ext_multiply(FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(),SchubertBasis,SchubertBasis), interim_result)

    #result = ring.from_dict({k: v for k, v in result.items() if len(k[1][0][0])<=n and (all([descent_pickle(k[1][0],dp), descent_pickle(k[1][1],dp)])) })
    print(result)

    # for ((perm1, (perm2, a)), (perm3, perm4)), coeff in result.items():
    #     assert perm1 == perm2 
    
    exit()
    #for (((perm1, a), (perm2, b)), (perm3, perm4)), coeff in result.items():
    #for (perm1, (((perm2, a),(perm3,b)), perm4)), coeff in result.items():
    for (perm1, (((perm2, a)), perm4)), coeff in result.items():
        perm3 = Permutation([])
        bong = Sx(perm2) * Sx(perm3)
        if any(len(permm)>n for permm in bong):
            continue
        try:
            assert perm1 != perm4 or bong.get(perm1, 0) == coeff, f"{bong=} {perm1=} {perm2=} {perm3=} {perm4=} {coeff=}"
        except AssertionError:
            print(f"Assertion failed: {bong.get(perm1,0)=} {perm1=} {perm2=} {perm3=} {perm4=} {coeff=}")
            continue
        
    exit()

    # result = TensorModule()
    # one_mod = RCGraphModule({RCGraph(): 1})
    # free_ring = FreeAlgebra(SchubertBasis)
    # ring2 = TensorRing(Sx([]).ring,SingleSchubertRing(y))
    # ring0 = TensorRing(ASx([]).ring, Sx([]).ring)
    # #DSx((FA(*seq)*one_mod).polyvalue(x,y))
    # tomo = ring0.one
    # w0 = uncode(list(range(n-1,0,-1)))
    # bismark = ring0.zero
    # # elem
    # for perm in perms:
    #     bismark += ring0.ext_multiply(ASx(perm,n-1), ASx(perm,n-1)*Sx(w0))
    
    # print(bismark)
    # exit()

    # for seq in seqs:
    #     #result += TensorModule.ext_multiply(FA(*seq).change_basis(SchubertBasis),Sx((FA(*seq)*one_mod).polyvalue(x,y)))
    #     stink = FA(*seq)*one_mod
    #     result += TensorModule.ext_multiply(FA(*seq).coproduct()*TensorModule({RCGraphTensor(RCGraph(),RCGraph()): 1}),(Sx(prod([x[i+1]**seq[i] for i in range(len(seq))]))*one_mod))
    #     rseq = tuple(reversed(seq))
    #     #result += TensorModule.ext_multiply(FA(*seq).change_basis(SchubertBasis),Sx(prod(E(rseq[i],rseq[i],x[1:],[y[n-1-i] for j in range(10)]) for i in range(len(seq)))))

    # result2 = TensorModule()

    # for ((perm1, a), perm), coeff in result.items():
    #     if len(perm1) > n or len(perm.perm) > n:
    #         continue
    #     if expand(coeff) == S.Zero:
    #         continue
    #     print(f"{perm1.trimcode=} {perm.perm.trimcode=}")
    # #exit()
    
    # # #exit()
    #     #result2 += TensorModule.ext_multiply(ASx(perm1, a),ring2.ext_multiply(Sx(perm),ring2.rings[1](coeff)))
    #     #print(f"{perm1.trimcode=} {perm.perm.trimcode=} {coeff=}")
    # exit()
    # for ((perm1, a), (perm2, perm3)), coeff in result2.items():
    #     print(f"{perm1.trimcode=} {perm2.trimcode=} {perm3.trimcode=} {coeff=}")
    # exit()
    ring = TensorRing(FreeAlgebra(ElementaryBasis), Sx([]).ring)
    ring0 = TensorRing(FreeAlgebra(SchubertBasis), Sx([]).ring)
    ring2 = TensorRing(Sx([]).ring, ring0)
    ring3 = TensorRing(FreeAlgebra(WordBasis), ring2)
    EE=FreeAlgebra(ElementaryBasis)
    zeppoli = ring.zero

    dingbat = ring3.zero
    u = uncode([2,0,1,3])
    doink_poly = Sx([])

    #DSx x and y
    # mul inserts column middle
    result = TensorModule()
    for seq in seqs:
        dingmod = FA(*seq).change_basis(SchubertBasis)
        # stinky = list(dingmod.keys())
        # for rc_graph in stinky:
        #     rt = rc_graph.transpose()
        #     vec = rt.length_vector()
            
        other_elem = Sx([]) * prod([x[i+1]**seq[i] for i in range(len(seq))])
        result += TensorModule.ext_multiply(dingmod, other_elem)
        # elem1 = ring3.rings[0](*seq)
        # elem2 = doink_poly * prod([x[i+1] ** (seq[i]) for i in range(len(seq))])
        # #elem2 = doink_poly * prod([E(seq[i],i+1,x[1:],z[1:]) for i in range(len(seq))])
        # #zeppoli += ring.ext_multiply(elem1, elem2)
        # seq = tuple(reversed(seq))
        # dingbat2 = Sx([])*prod([E(seq[i],i+1,x[1:],[0,0,0,0,0,0,0,0]) for i in range(len(seq))])
        # dingbat1 = EE(seq,n-1)
        # stinker = ring.ext_multiply(dingbat1, dingbat2)
        # farter = change_free_tensor_basis(stinker, ElementaryBasis,SchubertBasis)
        # dingbat += ring3.ext_multiply(elem1, ring2.ext_multiply(elem2,farter))

    # print(zeppoli)
    for (graph,perm), coeff in result.items():
        if len(perm) > n or len(graph[0]) > n:
            continue
        print(perm.trimcode,graph[0].trimcode)
        print(coeff)
        print(graph)
    # exit()
    exit()
    rng = SingleSchubertRing(y)
    da_baby = {}#TensorModule()
    for (rc_graph, schub), coeff in result.items():
        if len(schub) > n or len(rc_graph.perm) > n:
            continue
        da_baby[schub] = da_baby.get(schub, RCGraphModule()) + RCGraphModule({rc_graph: coeff})
    
    da_baby2 = {}
    for schub, module in da_baby.items():
        print(schub.trimcode)
        print(module)

    for (perm1, perm2), mod in da_baby2.items():
        print(perm1.trimcode)
        print(perm2.trimcode)
        print(mod)
    
    exit()
    zop = change_free_tensor_basis(dingbat, WordBasis, SchubertBasis)
    for zingbat, doofy in zop.items():
        if len(zingbat[0][0])>n or doofy == 0:
            continue
        print(f"{zingbat=}")
        print(doofy)
    print("PAINTED")
    # for pants, food in dingbat.items():
    #     if food == 0:
    #         continue
    #     print(f"{pants[0].trimcode} {pants[1].trimcode}")
    #     print(food)

    exit()
    modmod = RCGraphModule({RCGraph(): 1})

    addup = TensorModule()

    for perm in perms:
        addup += TensorModule.ext_multiply(ASx(perm,n-1).coproduct(), RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(perm,n-1), 1)))

    save_bongbong = {}

    for (elem, rc_graph), coeff in addup.items():
        if len(elem[0][0]) <= n and len(elem[1][0]) <= n and len(rc_graph.perm) <= n:
            print(f"{coeff=}")
            save_bongbong[elem] = save_bongbong.get(elem, RCGraphModule()) + RCGraphModule({rc_graph: 1})
            assert (Sx(elem[0][0])*Sx(elem[1][0])).get(rc_graph.perm,0) == coeff

    for dinglestick, mod in save_bongbong.items():
        print(dinglestick)
        print(mod)
    exit()

    ring1 = FreeAlgebra(SchubertBasis)
    ring2 = Sx([]).ring @ Sx([]).ring

    ring = TensorRing(ring1, ring2)

    
    addup = ring.zero
    print(seqs)
    for seq in seqs:
        fa_elem = FA(*seq).change_basis(SchubertBasis)
        fa_elem = ring1.from_dict({k: v for k,v in fa_elem.items() if len(k[0])<=n})
        ring2_elem = Sx(prod([x[i+1]**seq[i] for i in range(len(seq))]))
        ring3_elem = Sx([])*(prod([E(n-1-i-seq[i],n-1-i,x[1:],[0,0,0,0,0]) for i in range(len(seq))]))
        addup += ring.ext_multiply(fa_elem, ring2.ext_multiply(ring2_elem, ring3_elem))


    print(addup)


    save_dict = {}

    for (key1, (key2, key3)), coeff in addup.items():
        if len(key1[0]) > n:
            continue
        save_dict[key1] = save_dict.get(key1,0) + coeff * ring2((key2,key3))
        # assert coeff >=0 and coeff<=1
        # assert coeff == 0 or key1[0] == key2, f"{key1=}, {key2=} {key3=}"
        #assert coeff == 0 or ((~key2)*key3).trimcode == list(range(n-1,0,-1)), f"{key2.trimcode=} {key3.trimcode=}"
    for key1 in save_dict:
        print(f"{key1=}")
        print(f"{save_dict[key1]=}")
    
    # for i, perm1 in enumerate(perms):
    #     poly1 = Sx(perm1).expand()
    #     for perm2 in perms[i:]:
    #         poly2 = Sx(perm2).expand()
    #         print(f"{perm1.trimcode=}, {perm2.trimcode=}")

    #         rc_bungbat = schubert_positive_product(Sx(perm1).expand(), Sx(perm2).expand(), x, perm1.inv + perm2.inv, max(len(perm1.trimcode), len(perm2.trimcode)), check=False)

    #         real_result = Sx(perm1) * Sx(perm2)
    #         print(rc_bungbat)
            # for rc, val in rc_bungbat.items():
            #     #assert real_result.get(rc.perm, 0) == val
            #     print(rc)

    # from . import ASx
    # from .variables import genset_dict_from_expr
    # n = 4
    # perms = Permutation.all_permutations(n)
    # for perm1 in perms:
    #     for perm2 in perms:
    #         length = max(len(perm1.trimcode), len(perm2.trimcode))
    #         noimg = (Sx(perm1)*Sx(perm2))
    #         if any(len(k) > n for k in noimg.keys()):
    #             continue
    #         print(f"{perm1.trimcode}, {perm2.trimcode}")

    #         pl = noimg.expand()
    #         mod = RCGraphModule.from_poly(pl, Sx([]).ring.genset, length)
    #         # pare out
    #         mod00 = RCGraphModule(mod)
    #         mod2 = RCGraphModule()
    #         # found = True
    #         # while found:
    #         #     found = False
    #         #     st = set(mod2.keys())
    #         #     gr = min(st)
    #         #     pm = gr.perm
    #         pl2 = pl
    #         dm = RCGraphModule({RCGraph(): 1})
    #         for perm in perms:
    #             if len(perm.trimcode) > length:
    #                 continue
    #             try:
    #                 if perm not in {k.perm for k in mod.keys()}:
    #                     continue
    #                 assert min([v for k, v in mod.items() if k.perm == perm]) == noimg.get(perm, 0)
    #             except AssertionError:
    #                 # print(f"Failure {perm.trimcode=}")
    #                 # print(RCGraphModule({k: v for k, v in mod.items() if k.perm == perm}))
    #                 assert noimg.get(perm, 0) == 0
    #                 print('OK')
    #             if noimg.get(perm, 0) != 0:
    #                 print(f"Success {perm.trimcode=}")


    #         # while expand(pl2) != S.Zero:
    #         #     polyd = genset_dict_from_expr(pl2, Sx([]).ring.genset)
    #         #     mx = min(polyd.keys())
    #         #     pmm = uncode(mx)
    #         #     val = polyd[mx]
    #         #     mm = RCGraphModule({rc: val for rc, vll in mod.items() if rc.perm == pmm})
    #         #     mod00 -= mm
    #         #     mod2 += mm
    #         #     assert all(v > 0 for v in mod00.values())
    #         #     pl2 -= val * Sx(pmm).expand()
    #         # print(mod2)
    #         # assert expand(mod2.polyvalue(x) - pl) == S.Zero
            


# dual
# if __name__ == "__main_!!_":
#     from itertools import product

#     from schubmult.abc import x, y
#     from schubmult.symbolic import prod

#     mx = 3
#     length = 2
#     rd = TensorModule({RCGraphTensor(RCGraph(()),RCGraph(())): 1})

#     drd = DualTensorModule()
#     ding = TensorModule()
#     FA = FreeAlgebra(basis=WordBasis)
#     ASx = FreeAlgebra(basis=SchubertBasis)

#     lc = product(range(mx + 1), repeat=length)
#     for c in lc:
#         monom = prod([y[i + 1]**c[i] for i in range(length)])
#         #oil =  FA(uncode([i]),1).coproduct() * FA(uncode([j]),1).coproduct()
#         #oil = ring.ext_multiply(FAW(i,j),FAW(i,j).coproduct())
#         soing = tuple(c)
#         oil1 = FA(*soing)
#         oil2 = FA(*soing).coproduct()
#         ding += oil2* rd



#     n = 5
#     perms = list(Permutation.all_permutations(n))

#     for perm1 in perms:
#         if len(perm1.trimcode) > length:
#             continue
#         for perm2 in perms:
#             if len(perm2.trimcode) > length:
#                     continue
#             drd = DualTensorModule({k: v for k, v in ding.items() if {k[0].perm,k[1].perm} == {perm1, perm2}})
#             #drd = DualTensorModule(TensorModule.ext_multiply(DualRCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(perm1),1)), DualRCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(perm2),1))))
#             for perm in perms:
#                 if len(perm.trimcode) > length:
#                     continue
#                 if perm1.inv + perm2.inv != perm.inv:
#                     continue
#                 if perm not in {k[0].perm for k in ding.keys()}:
#                     continue

#                 drd2 = drd * ASx(perm, length).coproduct()
#                 if drd2 != DualTensorModule():
#                     print(f"{perm1.trimcode=}, {perm2.trimcode=}, {perm.trimcode=}")
#                     print("initial")
#                     print(drd)
#                     print("after")
#                     print(drd2)

    # for perm in perms:
    #     if perm.inv < 2:
    #         continue
    #     if len(perm.trimcode) > length:
    #         continue
    #     vm = drd * ASx(perm, length).coproduct()
    #     print(f"{perm.trimcode=}")
    #     for sudabaker, vvv in vm.items():
    #         print(f"{sudabaker=}")
    #         print(f"{vvv=}")
    #         print()
    # ring1 = ASx @ ASx
    # print(rd)
    # print()
    # print(rd*ASx(uncode([1,3,2])))



# if __name__ == "_pamko__":
#     from itertools import combinations_with_replacement

#     from schubmult import Sx
#     from schubmult.abc import x
#     from schubmult.rings.free_algebra_basis import *
#     FA = FreeAlgebra(basis=SchubertBasis)
#     FAW = FreeAlgebra(basis=WordBasis)

#     #spug = TensorModule({RCGraphTensor(RCGraph(()), RCGraph(())): 1})
#     spug = RCGraphModule({RCGraph(()): 1})
    
#     # spug1 = FA(uncode([0,1]), 2).coproduct() * spug
#     # print(spug1)

#     # spug2 = FA(uncode([1,0]), 2).coproduct() * spug
#     # print(spug2)
#     #resspug = TensorModule({RCGraphTensor(RCGraph(()), RCGraphTensor(RCGraph(), RCGraph(()))): 1})

#     res = RCGraphModule({RCGraphTensor(RCGraph(()), RCGraph(())): 1})
#     #spurg = (FA @ FA).zero
#     bungalo = {}
#     fratboy = {}
#     ring = FAW @ (FAW @ FAW)
#     mx = 3
#     length = 3

#     for c in combinations_with_replacement(range(mx + 1), length):
#         #oil =  FA(uncode([i]),1).coproduct() * FA(uncode([j]),1).coproduct()
#         #oil = ring.ext_multiply(FAW(i,j),FAW(i,j).coproduct())
#         oil = FAW(*c)
#         #spud_dog = oil.change_basis(SchubertBasis)
#         fratboy[tuple(c)] = fratboy.get(tuple(c), RCGraphModule()) + oil * res
#         # for k, v in spud_dog.items():
#         #     fratboy[k] = fratboy.get(k, RCGraphModule()) + v*spug2
#         #print(f"{FreeAlgebraBasis.change_tensor_basis(oil,WordBasis,WordBasis)=}")
#         # print(spug2)
#         # print(f"{dict(spug2)=}")
#         #res += spug2
#         #spurg += oil
#     foingle = {}
#     for k, v in fratboy.items():
#         foingle2 = FAW(*k).change_basis(SchubertBasis)
#         for kk, vv in foingle2.items():
#             foingle[kk] = foingle.get(kk, RCGraphModule()) + vv * v
#     # print(res)
#     # for bacon, v in fratboy.items():
#     #     print(f"{bacon=}")
#     #     print(v)
#     Permutation.print_as_code = True
#     for k, v in foingle.items():
#         if any(a > mx for a in k[0].trimcode):
#             continue
#         print(k[0].trimcode)
#         # print(v)
#         # print(v.polyvalue(x))
#         sm0 = S.Zero
#         sm1 = S.Zero
#         side0 = {}
#         side1 = {}
#         plathbucket = {}
#         graphd = RCGraphModule()

#         for kk, vv in v.items():
#             if (Sx(kk[0].perm)*Sx(kk[1].perm)).get(k[0], 0) == 0:
#                 continue
#             #     print("EXCLUDE")
#             #     print(kk)
#             #     print("END EXCLUDE")
#             #     continue
#             graphd += RCGraphModule({kk: vv})

#             key = (kk[0].perm, kk[1].perm)
#             side0[key[0]] = side0.get(key[0], 0) + vv*kk[0].polyvalue(x)
#             side1[key[1]] = side1.get(key[1], 0) + vv*kk[1].polyvalue(x)
#             plathbucket[key] = plathbucket.get(key, 0) + vv*kk[0].polyvalue(x)*kk[1].polyvalue(x)
#             # if kk[0].perm.inv ==0:
#             #     sm1 += vv * kk[1].polyvalue(x)
#             # elif kk[1].perm.inv ==0:
#             #     sm0 += vv * kk[0].polyvalue(x)

#             realval = (Sx(kk[0].perm)*Sx(kk[1].perm)).get(k[0], 0)
#             if realval > vv:
#                 print(f"Less  {vv} < {realval} {kk[0].perm.trimcode}, {kk[1].perm.trimcode}")
#             elif realval < vv:
#                 print(f"More  {vv} > {realval} {kk[0].perm.trimcode}, {kk[1].perm.trimcode}")
#             else:
#                 print(f"Equal {vv} = {realval} {kk[0].perm.trimcode}, {kk[1].perm.trimcode}")
#         # for perm, v in side0.items():
#         #     print(f"  L {perm.trimcode} : {v} {Sx(v)}")
#         # for perm, v in side1.items():
#         #     print(f"  R {perm.trimcode} : {v} {Sx(v)}")
#         # for (p0, p1), v in plathbucket.items():
#         #     print(f"  P {p0.trimcode}, {p1.trimcode} : {v} {Sx(v)}")
#         print(graphd)
        
#         print("-----")


# if __name__ == "__mainbubbles__":
#     from itertools import product

#     from schubmult import Sx
#     from schubmult.abc import x
#     from schubmult.rings.free_algebra_basis import *
#     FA = FreeAlgebra(basis=SchubertBasis)
#     FAW = FreeAlgebra(basis=WordBasis)

#     #spug = TensorModule({RCGraphTensor(RCGraph(()), RCGraph(())): 1})
#     spug = RCGraphModule({RCGraph(()): 1})
    
#     # spug1 = FA(uncode([0,1]), 2).coproduct() * spug
#     # print(spug1)

#     # spug2 = FA(uncode([1,0]), 2).coproduct() * spug
#     # print(spug2)
#     #resspug = TensorModule({RCGraphTensor(RCGraph(()), RCGraphTensor(RCGraph(), RCGraph(()))): 1})

#     res = RCGraphModule({RCGraph(()): 1})
#     #spurg = (FA @ FA).zero
#     fratboy = {}
#     ring = FAW
#     mx = 3
#     length = 3
#     lc = product(range(mx + 1), repeat=length)
#     for c in lc:
#         #oil =  FA(uncode([i]),1).coproduct() * FA(uncode([j]),1).coproduct()
#         #oil = ring.ext_multiply(FAW(i,j),FAW(i,j).coproduct())
#         soing = tuple(c)
#         print(c)
#         oil = FAW(*soing)
#         #spud_dog = oil.change_basis(SchubertBasis)
#         # dual module
#         modulebigot = oil * res
#         #tinseltown = fratboy.get(soing, ((), RCGraphModule()))
#         #fratboy[soing] = (soing, tinseltown[1] + oil * res)
#         # for k, v in spud_dog.items():
#         #     fratboy[k] = fratboy.get(k, RCGraphModule()) + v*spug2
#         #print(f"{FreeAlgebraBasis.change_tensor_basis(oil,WordBasis,WordBasis)=}")
#         # print(spug2)
#         # print(f"{dict(spug2)=}")
#         #res += spug2
#         #spurg += oil
#     # foingle = fratboy
#     # for k, v in fratboy.items():
#     #     foingle2 = FAW(*k).change_basis(SchubertBasis)
#     #     for kk, vv in foingle2.items():
#     #         if vv
#     #         foingle[kk] = foingle.get(kk, RCGraphModule()) + vv * v
#     # print(res)
#     # for bacon, v in fratboy.items():
#     #     print(f"{bacon=}")
#     #     print(v)

#     foingle = {}
    
#     for k, v in fratboy.items():
#         for graph, val in v.items():
#             foingle[(graph.perm,length)] = foingle.get((graph.perm,length), RCGraphModule()) + RCGraphModule({graph: val})

#     Permutation.print_as_code = True
#     for k, v in foingle.items():
#         sumbag = v.polyvalue(x)
#         for kk, vv in v.items():
#             assert k[0] == kk.perm, f"Failed {k} got {kk}, expected only {k[0]}"
#         # for kk, vv in v.items():
#         #     assert (vv == S.One and k[0].inv == 0) or k[0] == vv.perm
#         #     if vv != S.One:
#         #         sumbag += vv.polyvalue(x)
#         #     else:
#         #         sumbag += vv
#         try:
#             assert expand(sympify(sumbag) - Sx(k[0]).expand()) == S.Zero, f"Failed {k[0]}"
#         except Exception as e:
#             print(f"Failed {k[0].trimcode}")
#             continue
#         print(k[0].trimcode)
#         print(v)
#         print("-----")

# # if __name__ == "goodstuff":
# #     from itertools import product

# #     from schubmult import Sx
# #     from schubmult.abc import x
# #     from schubmult.rings.free_algebra_basis import *
# #     FA = FreeAlgebra(basis=SchubertBasis)
# #     FAW = FreeAlgebra(basis=WordBasis)

# #     #spug = TensorModule({RCGraphTensor(RCGraph(()), RCGraph(())): 1})
# #     spug = RCGraphModule({RCGraph(()): 1})
    
# #     # spug1 = FA(uncode([0,1]), 2).coproduct() * spug
# #     # print(spug1)

# #     # spug2 = FA(uncode([1,0]), 2).coproduct() * spug
# #     # print(spug2)
# #     #resspug = TensorModule({RCGraphTensor(RCGraph(()), RCGraphTensor(RCGraph(), RCGraph(()))): 1})

# #     res = TensorModule({RCGraphTensor(RCGraph(()), RCGraphTensor(RCGraph(), RCGraph())): 1})
# #     res0 = RCGraphModule({RCGraph(()): 1})
# #     res1 = TensorModule({RCGraphTensor(RCGraph(()), RCGraph(())): 1})
# #     #spurg = (FA @ FA).zero
# #     fratboy = {}
# #     ring = FAW @ (FAW @ FAW)
# #     mx = 3
# #     length = 3
# #     lc = product(range(mx + 1), repeat=length)
# #     for c in lc:
# #         #oil =  FA(uncode([i]),1).coproduct() * FA(uncode([j]),1).coproduct()
# #         #oil = ring.ext_multiply(FAW(i,j),FAW(i,j).coproduct())
# #         soing = tuple(c)
# #         print(c)
# #         oil1 = FAW(*soing)
# #         oil2 = FAW(*soing).coproduct()
# #         #spud_dog = oil.change_basis(SchubertBasis)


# #         fratboy[soing] = fratboy.get(soing, RCGraphModule()) + TensorModule.ext_multiply(oil1 * res0, oil2 * res1)
# #         # for k, v in spud_dog.items():
# #         #     fratboy[k] = fratboy.get(k, RCGraphModule()) + v*spug2
# #         #print(f"{FreeAlgebraBasis.change_tensor_basis(oil,WordBasis,WordBasis)=}")
# #         # print(spug2)
# #         # print(f"{dict(spug2)=}")
# #         #res += spug2
# #         #spurg += oil
# #     # foingle = fratboy
# #     # for k, v in fratboy.items():
# #     #     foingle2 = FAW(*k).change_basis(SchubertBasis)
# #     #     for kk, vv in foingle2.items():
# #     #         if vv
# #     #         foingle[kk] = foingle.get(kk, RCGraphModule()) + vv * v
# #     # print(res)
# #     # for bacon, v in fratboy.items():
# #     #     print(f"{bacon=}")
# #     #     print(v)

# #     foingle = {}
    
# #     for k, modbob in fratboy.items():
# #         badpairdct = {}
# #         for (graph0, (graph1, graph2)), val in modbob.items():
# #             if (graph1.perm, graph2.perm) in badpairdct.get(graph0.perm, set()):
# #                 continue
# #             #foingle[(graph0.perm,graph1.perm,graph2.perm)] = foingle.get((graph0.perm,graph1.perm,graph2.perm), TensorModule()) + TensorModule({RCGraphTensor(graph0, RCGraphTensor(graph1, graph2)): val})
# #             if expand(graph0.polyvalue(x) - graph1.polyvalue(x)*graph2.polyvalue(x)) != S.Zero:
# #                 continue
# #             good = True
# #             for i in range(length - 1, -1, -1):
# #                 g0 = graph0.rowrange(i, length)
# #                 g1 = graph1.rowrange(i, length)
# #                 g2 = graph2.rowrange(i, length)
# #                 if not g1.perm.bruhat_leq(g0.perm) or not g2.perm.bruhat_leq(g0.perm):
# #                     good = False
# #                     badpairdct[graph0.perm] = badpairdct.get(graph0.perm, set())
# #                     badpairdct[graph0.perm].add((graph1.perm, graph2.perm))
# #                     break
# #                 g0 = graph0.rowrange(i, i + 1)
# #                 g1 = graph1.rowrange(i, i + 1)
# #                 g2 = graph2.rowrange(i, i + 1)
# #                 if not g1.perm.bruhat_leq(g0.perm) or not g2.perm.bruhat_leq(g0.perm):
# #                     good = False
# #                     badpairdct[graph0.perm] = badpairdct.get(graph0.perm, set())
# #                     badpairdct[graph0.perm].add((graph1.perm, graph2.perm))
# #                     break

# #                 # g0 = RCGraph(graph0[i:i+1])
# #                 # g1 = RCGraph(graph1[i:i+1])
# #                 # g2 = RCGraph(graph2[i:i+1])
# #                 # if not g1.perm.bruhat_leq(g0.perm) or not g2.perm.bruhat_leq(g0.perm) or expand(g0.polyvalue(x) - g1.polyvalue(x)*g2.polyvalue(x)) != S.Zero:
# #                 #     good = False
# #                 #     break
# #             if good:
# #                 foingle[graph0] = foingle.get(graph0,TensorModule()) + TensorModule({RCGraphTensor(graph1, graph2): val})

# #     Permutation.print_as_code = True
# #     finglestick = foingle
# #     # for (k0, k1, k2), v in foingle.items():
        
# #     #     # for (kk0, kk1, kk2), vv in v.items():
# #     #     #     assert k0 == kk0.perm and k1 == kk1.perm and k2 == kk2.perm, f"Failed {k} got {kk}, expected only {k[0]}"
# #     #     # for kk, vv in v.items():
# #     #     #     assert (vv == S.One and k[0].inv == 0) or k[0] == vv.perm
# #     #     #     if vv != S.One:
# #     #     #         sumbag += vv.polyvalue(x)
# #     #     #     else:
# #     #     #         sumbag += vv
# #     #     # try:
# #     #     #     assert expand(sympify(sumbag) - Sx(k[0]).expand()) == S.Zero, f"Failed {k[0]}"
# #     #     # except Exception as e:
# #     #     #     print(f"Failed {k[0].trimcode}")
# #     #     #     continue
# #     #     # print(kk0.trimcode, kk1.trimcode, kk2.trimcode)
# #     #     for (kk0, (kk1, kk2)), vv in v.items():
# #     #         finglestick[kk0.perm] = finglestick.get(kk0.perm, TensorModule()) + TensorModule({RCGraphTensor(kk1, kk2): vv})

# #     for graph0, v in finglestick.items():
# #         print(f"Graph toinka: {graph0.perm.trimcode}")
# #         print(graph0)
# #         print(v)
# #         for (graph1, graph2), val in v.items():
# #             if (graph1.perm, graph2.perm) in badpairdct.get(graph0.perm, set()):
# #                 continue
# #             truth = (Sx(graph1.perm)*Sx(graph2.perm)).get(graph0.perm, 0) == val
# #             print(truth)
# #             if not truth:
# #                 print(RCGraphTensor(graph1, graph2))
# #                 gval = (Sx(graph1.perm)*Sx(graph2.perm)).get(graph0.perm, 0)
# #                 print(f"{val} != {gval}")
# #                 if gval != S.Zero:
# #                     raise Exception("Oh no")
# #         # if graph0.perm in permset:
# #         #     continue
# #         # permset.add(graph0.perm)
# #         # print(graph0.perm.trimcode)
# #         # print(v)
