# to encourage development

from bisect import bisect_left
from functools import cache

import sympy
from symengine import S, sympify
from sympy import Basic

import schubmult.rings._schubert_polynomial_ring as spr
import schubmult.rings._utils as utils
import schubmult.schub_lib.quantum as py
import schubmult.schub_lib.quantum_double as yz
from schubmult.perm_lib import Permutation, longest_element
from schubmult.poly_lib.poly_lib import elem_sym_poly, elem_sym_poly_q, xreplace_genvars
from schubmult.poly_lib.schub_poly import schubpoly_from_elems
from schubmult.poly_lib.variables import GeneratingSet, GeneratingSet_base
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import is_parabolic

## EMULATE POLYTOOLS

q_var = GeneratingSet("q")
# _def_printer = StrPrinter({"order": "none"})

logger = get_logger(__name__)


class QuantumDoubleSchubertAlgebraElement(spr.BasisSchubertAlgebraElement):
    def __new__(cls, _dict, basis):
        return spr.BasisSchubertAlgebraElement.__new__(cls, _dict, basis)

    def subs(self, old, new):
        logger.debug("ferefef")
        elb = self.as_classical().subs(old, new).as_quantum()
        logger.debug(f"{elb=}")
        return elb


# # TODO: not a noncommutative symbol, something else
# # Atomic Schubert polynomial
class QDSchubPoly(QuantumDoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, basis):
        return QDSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj = QuantumDoubleSchubertAlgebraElement.__new__(_class, sympy.Dict({(Permutation(k[0]), k[1]): 1}), basis)
        obj._perm = k[0]
        obj._coeff_var = k[1]
        # obj._base_var = base_var
        return obj

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset):
        return QDSchubPoly.__xnew__(_class, k, genset)

    def _sympystr(self, printer):
        if self._coeff_var == 0 or self._coeff_var == utils.NoneVar:
            return printer.doprint(f"QS{self.genset.label}({printer.doprint(self._perm)})")
        return printer.doprint(f"QDS{self.genset.label}({printer.doprint(self._perm)}, {spr._varstr(self._coeff_var)})")


class ParabolicQuantumDoubleSchubertAlgebraElement(spr.BasisSchubertAlgebraElement):
    def __new__(cls, _dict, basis):
        return spr.BasisSchubertAlgebraElement.__new__(cls, _dict, basis)

    @property
    def index_comp(self):
        return self.basis.index_comp

    def kill_ideal(self):
        length = sum(self.index_comp)
        new_dict = {}
        for k, v in self.coeff_dict.items():
            if len(k[0]) <= length:
                new_dict[k] = v
        return self.basis._from_dict(new_dict)


class PQDSchubPoly(ParabolicQuantumDoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, basis):
        return PQDSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj = ParabolicQuantumDoubleSchubertAlgebraElement.__new__(_class, sympy.Dict({(Permutation(k[0]), k[1]): 1}), basis)
        obj._perm = k[0]
        # obj._perm._print_as_code = True
        obj._coeff_var = k[1]
        # obj._base_var = base_var
        return obj

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset):
        return PQDSchubPoly.__xnew__(_class, k, genset)

    def _sympystr(self, printer):
        if self._coeff_var == 0 or self._coeff_var == utils.NoneVar:
            return printer.doprint(f"QPS{self.genset.label}{(tuple(self.index_comp))}({printer.doprint(self._perm)})")
        return printer.doprint(f"QPDS{self.genset.label}{tuple(self.index_comp)}({printer.doprint(self._perm)}, {spr._varstr(self._coeff_var)})")


class QuantumDoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset):
        return Basic.__new__(cls, genset)

    def _from_dict(self, _dict):
        return QuantumDoubleSchubertAlgebraElement(_dict, self)

    @property
    def genset(self):
        return self.args[0]

    @cache
    def cached_product(self, u, v, va, vb):
        return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}

    def in_quantum_basis(self, elem):
        return elem

    def in_classical_basis(self, elem):
        result = S.Zero
        for k, v in elem.coeff_dict.items():
            result += v * self.quantum_as_classical_schubpoly(k[0], k[1])
        return result

    def classical_elem_func(self, coeff_var):
        basis = spr.DoubleSchubertAlgebraElement_basis(self.genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis([], coeff_var)
            if p < 0 or p > k:
                return basis(0, coeff_var)
            return (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2) + elem_func(p, k - 1, varl1, varl2) + q_var[k - 1] * elem_func(p - 2, k - 2, varl1, varl2)

        return elem_func

    @property
    def single_element_class(self):
        return QDSchubPoly

    @cache
    def quantum_as_classical_schubpoly(self, perm, coeff_var="y"):
        return schubpoly_from_elems(perm, self.genset, utils.poly_ring(coeff_var), self.classical_elem_func(coeff_var))

    @cache
    def cached_schubpoly(self, k):
        return schubpoly_from_elems(k[0], self.genset, utils.poly_ring(k[1]), elem_func=elem_sym_poly_q)  # yz.schubpoly_quantum(k[0], self.genset, utils.poly_ring(k[1]))

    @cache
    def cached_positive_product(self, u, v, va, vb):
        return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()}

    @property
    def double_mul(self):
        return yz.schubmult_q_double_fast

    @property
    def single_mul(self):
        return py.schubmult_q_fast

    @property
    def mult_poly_single(self):
        return py.mult_poly_q

    @property
    def mult_poly_double(self):
        return yz.mult_poly_q_double

    def __call__(self, x, cv=None):
        genset = self.genset
        # logger.debug(f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            if cv is None:
                cv = "y"
            elem = self._from_dict({(Permutation(x), cv): 1})
        elif isinstance(x, Permutation):
            if cv is None:
                cv = "y"
            elem = self._from_dict({(x, cv): 1})
        # elif isinstance(x, spr.SchubertPolynomial):
        #     if x._parent._base_var == self._base_var:
        #         elem_dict = {(x, utils.NoneVar): v for k, v in x.coeff_dict.items()}
        #         elem = QuantumDoubleSchubertAlgebraElement(elem_dict, self)
        #         if cv is not None:
        #             elem = self([1, 2], cv) * elem
        #     else:
        #         return self(x.expand(), cv)
        elif isinstance(x, QuantumDoubleSchubertAlgebraElement):
            if x.is_Add or x.is_Mul:
                return x
            if x.genset == genset:
                elem = QuantumDoubleSchubertAlgebraElement(x.coeff_dict, self)  # , self)
            else:
                return self(x.expand(), cv, genset)
        elif isinstance(x, spr.DoubleSchubertAlgebraElement):
            if x.genset == self.genset:
                return self(x.expand(), cv, genset)
        else:
            logger.debug("bagelflap")
            x = sympify(x)
            if cv is None or cv == utils.NoneVar:
                cv = utils.NoneVar
                logger.debug(f"{x=} {list(genset)=}")
                result = py.mult_poly_q({Permutation([]): 1}, x, genset)
                logger.debug(f"{result=}")
            else:
                result = yz.mult_poly_q_double({Permutation([]): 1}, x, genset, utils.poly_ring(cv))
            elem = QuantumDoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()}, self)
        return elem


QDSx = QuantumDoubleSchubertAlgebraElement_basis(GeneratingSet("x"))

t = GeneratingSet("t")
a = GeneratingSet("a")

spunky_basis = spr.SchubertAlgebraElement_basis(t)


class QuantumSchubertAlgebraElement_basis(QuantumDoubleSchubertAlgebraElement_basis):
    def __new__(cls, genset):
        return QuantumDoubleSchubertAlgebraElement_basis.__new__(cls, genset)

    def _from_single_dict(self, _dict):
        return QuantumDoubleSchubertAlgebraElement({(k, utils.NoneVar): v for k, v in _dict.items()}, self)

    def __call__(self, x):
        genset = self.genset
        # logger.debug(f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self._from_single_dict({Permutation(x): 1})
        elif isinstance(x, Permutation):
            elem = self._from_single_dict({x: 1})
        # elif isinstance(x, spr.SchubertPolynomial):
        #     if x._parent._base_var == self._base_var:
        #         elem_dict = {(x, utils.NoneVar): v for k, v in x.coeff_dict.items()}
        #         elem = QuantumDoubleSchubertAlgebraElement(elem_dict, self)
        #         if cv is not None:
        #             elem = self([1, 2], cv) * elem
        #     else:
        #         return self(x.expand(), cv)
        elif isinstance(x, QuantumDoubleSchubertAlgebraElement):
            if x.is_Add or x.is_Mul:
                return x
            if x.genset == genset:
                elem = QuantumDoubleSchubertAlgebraElement(x.coeff_dict, self)  # , self)
            else:
                return self(x.expand())
        elif isinstance(x, spr.DoubleSchubertAlgebraElement):
            if x.genset == self.genset:
                return x.as_quantum()
        elif isinstance(x, ParabolicQuantumDoubleSchubertAlgebraElement):
            return x.as_quantum()
        else:
            x = sympify(x)
            result = py.mult_poly_q({Permutation([]): 1}, x, genset)
            elem = self._from_single_dict(result)
        return elem


class ParabolicQuantumDoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset, index_comp):
        obj = Basic.__new__(cls, genset, tuple(index_comp))
        obj._quantum_basis = QuantumDoubleSchubertAlgebraElement_basis(genset)
        obj._classical_basis = spr.DoubleSchubertAlgebraElement_basis(genset)
        obj._n = list(index_comp)
        obj._N = [sum(obj._n[:i]) for i in range(len(obj._n) + 1)]
        # print(f"{obj._N=}")
        # obj._D = []
        # obj._E = {}
        # from symengine import Matrix
        # for j in range(1, len(obj._N)):
        #     m_arr = [[0 for i in range(obj._N[j])] for p in range(obj._N[j])]
        #     for i in range(obj._N[j]):
        #         m_arr[i][i] = a[i+1] - t[1] #genset[i+1] - t[1]
        #         if i < obj._N[j] - 1:
        #             m_arr[i][i+1] = -1
        #     for b in range(1, j):
        #         njm1 = obj._N[b + 1] - 1
        #         njp1 = obj._N[b - 1]
        #         # print(f"{b=}")
        #         # print(f"{njm1=} {njp1=}")
        #         if njp1 < obj._N[j] and njm1 < obj._N[j]:
        #             # print(f"{b=} {obj._n[b]=}")
        #             m_arr[njm1][njp1] = -(-1)**(obj._n[b])*q_var[b]
        #     # print(Matrix(m_arr))
        #     poly = Matrix(m_arr).det().simplify()
        #     # print(f"{poly=}")
        #     # def dongle(v):
        #     #     return poly.subs(t[1], v)
        #     obj._D += [spunky_basis(poly)]
        #     obj._E[obj._N[j]] = {obj._N[j]: obj._D[-1]}
        #     for i in range(1,obj._N[j]):
        #         obj._E[obj._N[j]][obj._N[j] - i] = -obj._E[obj._N[j]][obj._N[j] - i + 1].divdiff(i)
        # # print(obj._E)
        # add_am = 6
        # index_comp += [add_am]
        # obj._N += [obj._N[-1] + add_am]
        parabolic_index = []
        start = 0
        # 1, 2 | 3
        for i in range(len(index_comp)):
            end = start + index_comp[i]
            parabolic_index += list(range(start + 1, end))
            # start += int(args.parabolic[i])
            start = end
        obj._parabolic_index = parabolic_index
        obj._otherlong = Permutation(list(range(obj._N[-1], 0, -1)))
        obj._longest = obj._otherlong * longest_element(parabolic_index)
        # obj._E[0] = obj._from_dict({(Permutation([]),utils.NoneVar): S.One})
        return obj

    @property
    def parabolic_index(self):
        return self._parabolic_index

    @property
    def quantum_basis(self):
        return self._quantum_basis

    @property
    def classical_basis(self):
        return self._classical_basis

    # def elem_sym_poly(self, p, k, varl1, varl2, xstart=0, ystart=0):
    #     # print(f"{p=} {k=} {xstart=} {ystart=} {len(varl1)=} {len(varl2)=}")
    #     if p > k:
    #         return zero
    #     if p == 0:
    #         return one
    #     if p == 1:
    #         res = varl1[xstart] - varl2[ystart]
    #         for i in range(1, k):
    #             res += varl1[xstart + i] - varl2[ystart + i]
    #         return res
    #     if p == k:
    #         res = (varl1[xstart] - varl2[ystart]) * (varl1[xstart + 1] - varl2[ystart])
    #         for i in range(2, k):
    #             res *= varl1[i + xstart] - varl2[ystart]
    #         return res
    #     mid = k // 2
    #     xsm = xstart + mid
    #     ysm = ystart + mid
    #     kmm = k - mid
    #     res = elem_sym_poly(p, mid, varl1, varl2, xstart, ystart) + elem_sym_poly(
    #         p,
    #         kmm,
    #         varl1,
    #         varl2,
    #         xsm,
    #         ysm,
    #     )
    #     for p2 in range(max(1, p - kmm), min(p, mid + 1)):
    #         res += elem_sym_poly(p2, mid, varl1, varl2, xstart, ystart) * elem_sym_poly(
    #             p - p2,
    #             kmm,
    #             varl1,
    #             varl2,
    #             xsm,
    #             ysm - p2,
    #         )
    #     return res

    def elem_sym(self):
        def elem_func(p, k, varl1, varl2):
            # print(f"{p=} {k=} {self._N=}")
            if p < 0 or p > k:
                return 0
            if p == 0 and k >= 0:
                return 1
            if k <= self._N[1]:
                return elem_sym_poly(p, k, varl1, varl2)
            ret = 0
            j = bisect_left(self._N, k)
            if j < len(self._N) and k == self._N[j]:
                ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            ret += elem_func(p, k - 1, varl1, varl2) + (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2)
            return ret

        return elem_func

    def boingle_elem_sym(self):
        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return S.One
            if p < 0 or p > k:
                return S.Zero
            # print(f"{p=} {k=}")
            spoink = self._E[k][p]
            # print(f"{spoink=}")
            # for i in range(1,k - p + 1):
            #     spoink = -spoink.divdiff(i)
            return sympify(sympify(spoink.as_polynomial()).xreplace({t[i]: varl2[i - 1] for i in range(1, len(varl2) + 1)})).xreplace({a[i]: varl1[i - 1] for i in range(1, len(varl1) + 1)})
            # TEMP
            # varl2 = utils.poly_ring(0)
            # if p < 0 or p > k:
            #     return 0
            # if p == 0 and k >= 0:
            #     return 1
            # if k == self._N[1]:
            #     return elem_sym_poly(p, k, varl1, varl2)
            # print(f"{p=} {k=} {self._N=}")
            # ret = 0
            # j = bisect_left(self._N, k)
            # if j < len(self._N) and k == self._N[j]:
            #     ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            #     #ret += elem_func(p, self._N[j-1], varl1, varl2)
            #     # print(f"{self._n[j-1]=}")
            #     for i in range(min(p+1,self._n[j-1]+1)):
            #         # print(f"{p=} {p - i=} {i=} {self._N[j-1]-1=} {self._N[j-1]-1-i=} {len(varl2)=}")
            #         ret += varl1[self._N[j] - 1 - i] * elem_func(p - i, self._N[j-1], varl1, varl2)
            # else:
            #     # print("Bob jones")
            # return ret

        return elem_func

    # def classical_elem(self, k, coeff_var):
    #     if k <= self._N[1]:
    #         return self(uncode([1 for i in range(k)]), coeff_var)

    def _from_dict(self, _dict):
        return ParabolicQuantumDoubleSchubertAlgebraElement(_dict, self)

    @property
    def genset(self):
        return self.args[0]

    @property
    def index_comp(self):
        return self.args[1]

    def process_coeff_dict(self, coeff_dict):
        max_len = max(len(w) for w in coeff_dict)
        parabolic_index = [*self._parabolic_index]
        # parabolic_index += list(range(parabolic_index[-1] + 2, max_len + 1))
        if max_len > len(self._longest):
            parabolic_index = []
            start = 0
            # 1, 2 | 3
            index_comp = [*self._n, max_len + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                # start += int(args.parabolic[i])
                start = end
        return yz.apply_peterson_woodward(coeff_dict, parabolic_index)

    @cache
    def cached_product(self, u, v, va, vb):
        initial_dict = {k: xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}
        return {(k, va): v for k, v in self.process_coeff_dict(initial_dict).items()}

    def in_quantum_basis(self, elem):
        result = S.Zero
        for k, v in elem.coeff_dict.items():
            result += v * schubpoly_from_elems(k[0], self.genset, utils.poly_ring(k[1]), self.quantum_elem_func(k[1]))
            # print(f"{result=}")
        return result

    def in_classical_basis(self, elem):
        result = S.Zero
        for k, v in elem.coeff_dict.items():
            result += v * self.quantum_as_classical_schubpoly(k[0], k[1])
            # print(f"{result=}")
        return result

    @cache
    def classical_in_basis(self, k):
        from symengine import expand

        a = self.classical_basis(*k)
        b = self(*k)
        if expand(a.as_polynomial() - b.as_polynomial()) == S.Zero:
            return b
        cd = dict(b.as_classical().coeff_dict)
        for k2, v in cd.items():
            if k != k2:
                b -= v * self.classical_in_basis(k2)
        return b

    def classical_elem_func(self, coeff_var):
        basis = spr.DoubleSchubertAlgebraElement_basis(self.genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            # print(f"{p=} {k=} {varl1=} {varl2=}")
            if p == 0 and k >= 0:
                return basis([], coeff_var)
            if p < 0 or p > k:
                return basis(0, coeff_var)
            if k <= self._N[1]:
                return basis(elem_sym_poly(p, k, varl1, varl2), coeff_var)
            ret = basis(0, coeff_var)
            j = bisect_left(self._N, k)
            if j < len(self._N) and k == self._N[j]:
                ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            ret += elem_func(p, k - 1, varl1, varl2) + (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2)
            # print(f"{ret=}")
            return ret

        return elem_func

    def quantum_elem_func(self, coeff_var):
        basis = QuantumDoubleSchubertAlgebraElement_basis(self.genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            # print(f"{p=} {k=} {varl1=} {varl2=}")
            if p == 0 and k >= 0:
                return basis([], coeff_var)
            if p < 0 or p > k:
                return basis(0, coeff_var)
            if k <= self._N[1]:
                return basis(elem_sym_poly(p, k, varl1, varl2), coeff_var)
            ret = basis(0, coeff_var)
            j = bisect_left(self._N, k)
            if j < len(self._N) and k == self._N[j]:
                ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            ret += elem_func(p, k - 1, varl1, varl2) + (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2)
            # print(f"{ret=}")
            return ret

        return elem_func

    @property
    def single_element_class(self):
        return PQDSchubPoly

    @cache
    def quantum_as_classical_schubpoly(self, perm, coeff_var="y"):
        k = (perm, coeff_var)
        # print(f"{k=}")
        if len(k[0]) > len(self._longest):
            parabolic_index = []
            start = 0
            # 1, 2 | 3
            index_comp = [*self._n, len(k[0]) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                # start += int(args.parabolic[i])
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation(longest_element(parabolic_index))
            # print(f"{longpar=} {parabolic_index=}")
            longest = otherlong * longpar
            # print(f"new longest = {longest=}")
        else:
            longest = self._longest
        return schubpoly_from_elems(perm, self.genset, utils.poly_ring(coeff_var), elem_func=self.classical_elem_func(coeff_var), mumu=~longest)

    @cache
    def cached_schubpoly(self, k):
        if len(k[0]) > len(self._longest):
            parabolic_index = []
            start = 0
            # 1, 2 | 3
            index_comp = [*self._n, len(k[0]) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                # start += int(args.parabolic[i])
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation(longest_element(parabolic_index))
            # print(f"{longpar=} {parabolic_index=}")
            longest = otherlong * longpar
            # print(f"new longest = {longest=}")
        else:
            longest = self._longest
        return schubpoly_from_elems(k[0], self.genset, utils.poly_ring(k[1]), elem_func=self.elem_sym(), mumu=~longest)  # yz.schubpoly_quantum(k[0], self.genset, utils.poly_ring(k[1]))

    @cache
    def cached_positive_product(self, u, v, va, vb):
        initial_dict = {k: xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()}
        return {(k, va): v for k, v in self.process_coeff_dict(initial_dict).items()}

    @property
    def double_mul(self):
        return yz.schubmult_q_double_fast

    @property
    def single_mul(self):
        return py.schubmult_q_fast

    @property
    def mult_poly_single(self):
        return py.mult_poly_q

    @property
    def mult_poly_double(self):
        return yz.mult_poly_q_double

    def __call__(self, x, cv=None):
        genset = self.genset
        # logger.debug(f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            if cv is None:
                cv = "y"
            perm = Permutation(x)
            if not is_parabolic(perm, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {perm} is not")
            elem = self._from_dict({(perm, cv): 1})
        elif isinstance(x, Permutation):
            if cv is None:
                cv = "y"
            if not is_parabolic(x, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {x} is not")
            elem = self._from_dict({(x, cv): 1})
        elif isinstance(x, ParabolicQuantumDoubleSchubertAlgebraElement):
            return x
        else:
            dct = self.classical_basis(x, cv)
            elem = 0
            if not isinstance(dct, spr.BasisSchubertAlgebraElement):
                return dct
            try:
                for k, v in dct.coeff_dict.items():
                    if elem == 0:
                        elem = v * self.classical_in_basis(k)
                    else:
                        elem += v * self.classical_in_basis(k)
            except ValueError:
                raise ValueError(f"Could not convert {x=} to quantum parabolic")
        return elem


QSx = QuantumSchubertAlgebraElement_basis(GeneratingSet("x"))

QuantumDoubleSchubertPolynomial = QuantumDoubleSchubertAlgebraElement


def make_parabolic_quantum_basis(index_comp):
    return ParabolicQuantumDoubleSchubertAlgebraElement_basis(GeneratingSet("x"), index_comp)


@cache
def QPDSx(*args):
    return make_parabolic_quantum_basis(args)


# is_Add = True
# is_Mul = True
# is_Add
# is_AlgebraicNumber
# is_Atom
# is_Boolean
# is_Derivative
# is_Dummy
# is_Equality
# is_Float
# is_Function
# is_Indexed
# is_Integer
# is_MatAdd
# is_MatMul
# is_Matrix
# is_Mul
# is_Not
# is_Number
# is_NumberSymbol
# is_Order
# is_Piecewise
# is_Point
# is_Poly
# is_Pow
# is_Rational
# is_Relational
# is_Symbol
# is_Vector
# is_Wild
# is_algebraic
# is_algebraic_expr
# is_antihermitian
# is_commutative
# is_comparable
# is_complex
# is_composite
# is_constant
# is_even
# is_extended_negative
# is_extended_nonnegative
# is_extended_nonpositive
# is_extended_nonzero
# is_extended_positive
# is_extended_real
# is_finite
# is_hermitian
# is_hypergeometric
# is_imaginary
# is_infinite
# is_integer
# is_irrational
# is_meromorphic
# is_negative
# is_noninteger
# is_nonnegative
# is_nonpositive
# is_nonzero
# is_number
# is_odd
# is_polar
# is_polynomial
# is_positive
# is_prime
# is_rational
# is_rational_function
# is_real
# is_scalar
# is_symbol
# is_transcendental
# is_zero
# is_polynomial = True
# is_Symbol = True


# class SchubAdd(QuantumDoubleSchubertAlgebraElement, Add):
#     is_Add = True

#     def __new__(cls, *args, evaluate=True, _sympify=True):
#         obj = Add.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)
#         if evaluate:
#             return obj.doit()
#         return obj

#     def doit(self):
#         ret = self.args[0]
#         for arg in self.args[1:]:
#             if arg.is_Add or arg.is_Mul:
#                 arg = arg.doit()
#             ret = _do_schub_add(ret, arg)
#         return ret

#     # def _sympystr(self, printer):
#     #     return _def_printer._print(f"SchubAdd({self.args}")


# class SchubMul(QuantumDoubleSchubertAlgebraElement, Mul):
#     is_Mul = True

#     def __new__(cls, *args, evaluate=True, _sympify=True):
#         if len(args) == 0:
#             return 1
#         # args, a, b = Mul.flatten(list(args))
#         # if len(args) == 0:
#         #     return 1
#         obj = Mul.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)

#         if evaluate:
#             return obj.doit()
#         return obj

#     def doit(self):
#         ret = self.args[0]
#         for arg in self.args[1:]:
#             if arg.is_Add or arg.is_Mul:
#                 arg = arg.doit()
#             ret = _do_schub_mul(ret, arg)
#         return ret


# Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
#     "Mul": [get_postprocessor(Mul)],
#     "Add": [get_postprocessor(Add)],
# }
