from bisect import bisect_left
from functools import cache

import sympy
from symengine import Add, Mul, Pow, S, Symbol, expand, sympify

import schubmult.utils.ring_utils as utils
import schubmult.rings.schubert_ring as spr
import schubmult.schub_lib.quantum as py
import schubmult.schub_lib.quantum_double as yz
from schubmult.perm_lib import Permutation, longest_element
from schubmult.poly_lib.poly_lib import complete_sym_poly, elem_sym_poly, elem_sym_poly_q, xreplace_genvars
from schubmult.poly_lib.schub_poly import schubpoly_from_elems
from schubmult.poly_lib.variables import GeneratingSet, GeneratingSet_base
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import is_parabolic

from .abstract_schub_poly import AbstractSchubPoly
from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing

q_var = GeneratingSet("q")

logger = get_logger(__name__)


class QuantumDoubleSchubertElement(BaseSchubertElement):
    def subs(self, old, new):
        return self.as_classical().subs(old, new).as_quantum()


_pretty_schub_char = "𝕼𝔖"  # noqa: RUF001


class QDSchubPoly(AbstractSchubPoly):
    is_Atom = True

    def __new__(cls, k, ring):
        return QDSchubPoly.__xnew_cached__(cls, k, ring)

    @staticmethod
    def __xnew__(_class, k, ring):
        return AbstractSchubPoly.__new__(_class, k, ring)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset):
        return QDSchubPoly.__xnew__(_class, k, genset)

    def _sympystr(self, printer):
        if self.ring.coeff_genset.label is None:
            return printer.doprint(f"QS{self.genset.label}({printer.doprint(self._perm)})")
        return printer.doprint(f"QDS{self.genset.label}({printer.doprint(self._perm)}, {self.ring.coeff_genset.label})")

    def _pretty(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        subscript = printer.doprint(int("".join([str(i) for i in self._key])))
        if self.ring.coeff_genset.label is None:
            return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(self.genset.label)))
        return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(f"{self.genset.label}; {self.ring.coeff_genset.label}")))

    def _latex(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        subscript = sympy.sstr(self._key)
        if self.ring.coeff_genset.label is None:
            return printer._print_Function(sympy.Function("\\widetilde{\\mathfrak{S}}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(self.genset.label)))
        return printer._print_Function(sympy.Function("\\widetilde{\\mathfrak{S}}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(f"{self.genset.label}; {self.ring.coeff_genset.label}")))


class ParabolicQuantumDoubleSchubertElement(BaseSchubertElement):
    @property
    def index_comp(self):
        return self.ring.index_comp

    def kill_ideal(self):
        length = sum(self.index_comp)
        new_dict = {}
        for k, v in self.items():
            if len(k) <= length:
                new_dict[k] = v
        return self.ring.from_dict(new_dict)


class PQDSchubPoly(AbstractSchubPoly):
    is_Atom = True

    def __new__(cls, k, basis, index_comp):
        return PQDSchubPoly.__xnew_cached__(cls, k, basis, index_comp)

    @staticmethod
    def __xnew__(_class, k, basis, index_comp):
        obj = AbstractSchubPoly.__new__(_class, k, basis)
        obj._perm = k
        obj._key = k
        obj._index_comp = index_comp
        return obj

    @property
    def index_comp(self):
        return self._index_comp

    @property
    def args(self):
        return (sympy.Tuple(*self._key), self._basis, sympy.Tuple(*self._index_comp))

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, basis, index_comp):
        return PQDSchubPoly.__xnew__(_class, k, basis, index_comp)

    def _sympystr(self, printer):
        if self.ring.coeff_genset.label is None:
            return printer.doprint(f"QPS{self.genset.label}{(tuple(self.index_comp))}({printer.doprint(self._perm)})")
        return printer.doprint(f"QPDS{self.genset.label}{tuple(self.index_comp)}({printer.doprint(self._perm)}, {self.ring.coeff_genset.label})")

    def _pretty(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(int("".join([str(i) for i in self._key])))
        if self.ring.coeff_genset.label is None:
            return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(f"{self.genset.label} | {self.ring.index_comp}")))
        return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(f"{self.genset.label}; {self.ring.coeff_genset.label} | {self.ring.index_comp}")))

    def _latex(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(self._key)
        supscript = printer._print(sympy.Tuple(*self.index_comp))
        if self.ring.coeff_genset.label is None:
            return printer._print_Function(sympy.Function("\\widetilde{\\mathfrak{S}}" + f"^{'{'}{supscript}{'}'}_{'{' + subscript + '}'}")(sympy.Symbol(self.genset.label)))
        return printer._print_Function(
            sympy.Function("\\widetilde{\\mathfrak{S}}" + f"^{'{'}{supscript}{'}'}_{'{' + subscript + '}'}")(sympy.Symbol(f"{self.genset.label}; {self.ring.coeff_genset.label}")),
        )


class QuantumDoubleSchubertRing(BaseSchubertRing):
    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "QDBS"))

    def __init__(self, genset, coeff_genset):
        super().__init__(genset, coeff_genset)
        self.dtype = type("QuantumDoubleSchubertElement", (QuantumDoubleSchubertElement,), {"ring": self})

    def __str__(self):
        return f"Quantum Double Schubert polynomial ring in {self.genset.label} and {self.coeff_genset.label}"

    def printing_term(self, k):
        return QDSchubPoly(k, self)

    def _coerce_mul(self, other):
        if isinstance(other, BaseSchubertElement):
            if type(other.ring) is type(self):
                if self.genset == other.ring.genset:
                    return other
            if type(other.ring) is QuantumSingleSchubertRing:
                if self.genset == other.ring.genset:
                    newbasis = QuantumDoubleSchubertRing(self.genset, utils.poly_ring(0))
                    return newbasis.from_dict(other)
        return None

    def _coerce_add(self, other):
        if isinstance(other, BaseSchubertElement):
            if type(other.ring) is type(self):
                if self.genset == other.ring.genset and self.coeff_genset == other.ring.coeff_genset:
                    return other
        return None

    @property
    def symbol_elem_func(self):
        def elem_func(p, k, varl1, varl2):  # noqa: ARG001
            if p == 0 and k >= 0:
                return 1
            if p < 0 or p > k:
                return 0
            return sympy.Add(*[(Symbol(f"e_{p - i}_{k}") if p - i > 0 else 1) * complete_sym_poly(i, k + 1 - p, [-v for v in varl2]) for i in range(p + 1)])

        return elem_func

    def elem_sym_subs(self, kk):
        elems = []
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(sympy.Symbol(f"e_{p}_{k}"), elem_sym_poly_q(p, k, self.genset[1:], utils.poly_ring(0)))]
        return dict(elems)

    @cache
    def cached_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}

    def in_quantum_basis(self, elem):
        return elem

    def in_classical_basis(self, elem):
        result = None
        for k, v in elem.items():
            if not result:
                result = v * self.quantum_as_classical_schubpoly(k)
            else:
                result += v * self.quantum_as_classical_schubpoly(k)
        return result if result else S.Zero

    @property
    def classical_elem_func(self):
        basis = spr.DoubleSchubertRing(self.genset, self.coeff_genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis([])
            if p < 0 or p > k:
                return basis(0)
            return (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2) + elem_func(p, k - 1, varl1, varl2) + q_var[k - 1] * elem_func(p - 2, k - 2, varl1, varl2)

        return elem_func

    @cache
    def quantum_as_classical_schubpoly(self, perm):
        return schubpoly_from_elems(perm, self.genset, self.coeff_genset, self.classical_elem_func)

    @cache
    def cached_schubpoly(self, k):
        return schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=elem_sym_poly_q)

    @cache
    def cached_positive_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()}

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

    def from_sympy(self, x):
        if isinstance(x, BaseSchubertElement):
            if x.ring == self:
                return x
        x = sympify(x)
        ind = self.genset.index(x)
        if ind != -1:
            return self.from_dict(yz.mult_poly_q_double({Permutation([]): 1}, x, self.genset, self.coeff_genset))
        if isinstance(x, Add):
            return self.sum([self.from_sympy(arg) for arg in x.args])
        if isinstance(x, Mul):
            res = self.one
            for arg in x.args:
                res *= self.from_sympy(arg)
            return res
        if isinstance(x, Pow):
            return self.from_sympy(x.args[0]) ** int(x.args[1])
        return self.from_dict({Permutation([]): x})

    # def from_sympy(self, x):
    #     x = sympify(x)
    #     result = yz.mult_poly_q_double({Permutation([]): 1}, x, self.genset, self.coeff_genset)
    #     return self.from_dict(result)

    def new(self, x):
        genset = self.genset
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): 1})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: 1})
        else:
            elem = self.from_sympy(x)
        return elem


def QDSx(x, genset=GeneratingSet("y")):
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    return QuantumDoubleSchubertRing(GeneratingSet("x"), genset)(x)


class QuantumSingleSchubertRing(QuantumDoubleSchubertRing):
    def __init__(cls, genset):
        super().__init__(genset, utils.poly_ring(utils.NoneVar))

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "QBS"))

    def _coerce_mul(self, other):
        """Coerce a basis schubert algebra element so it can be multiplied

        Args:
            other (_type_): _description_

        Returns:
            _type_: _description_
        """
        if type(other.ring) is type(self):
            if self.genset == other.ring.genset:
                return other
        if type(other.ring) is QuantumDoubleSchubertRing:
            if self.genset == other.ring.genset:
                return other
        return None

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            return py.schubmult_q_fast({u: S.One}, v)
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def from_sympy(self, x):
        if isinstance(x, BaseSchubertElement):
            if x.ring == self:
                return x
        x = sympify(x)
        ind = self.genset.index(x)
        if ind != -1:
            return self.from_dict(py.mult_poly_q({Permutation([]): 1}, x, self.genset))
        if isinstance(x, Add):
            return self.sum([self.from_sympy(arg) for arg in x.args])
        if isinstance(x, Mul):
            res = self.one
            for arg in x.args:
                res *= self.from_sympy(arg)
            return res
        if isinstance(x, Pow):
            return self.from_sympy(x.args[0]) ** int(x.args[1])
        return self.from_dict({Permutation([]): x})

    # def from_sympy(self, x):
    #     x = sympify(x)
    #     result = py.mult_poly_q({Permutation([]): 1}, x, self.genset)
    #     return self.from_dict(result)

    def new(self, x):
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): 1})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: 1})
        elif isinstance(x, spr.DoubleSchubertElement):
            if x.genset == self.genset:
                return x.as_quantum()
        elif isinstance(x, ParabolicQuantumDoubleSchubertElement):
            return x.as_quantum()
        else:
            elem = self.from_sympy(x)
        return elem

    def in_SEM_basis(self):
        result = S.Zero
        for k, v in self.items():
            if len(k) > len(self._longest):
                parabolic_index = []
                start = 0
                index_comp = [*self._n, len(k) + 1 - self._N[-1]]
                for i in range(len(index_comp)):
                    end = start + index_comp[i]
                    parabolic_index += list(range(start + 1, end))
                    start = end
                otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
                longpar = Permutation(longest_element(parabolic_index))
                longest = otherlong * longpar
            else:
                longest = self._longest
            result += v * schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=self.ring.symbol_elem_func, mumu=~longest)
        return result


class ParabolicQuantumDoubleSchubertRing(BaseSchubertRing):
    def __hash__(self):
        return hash((self.genset, self.coeff_genset, self.index_comp, "PBGBG"))

    def __init__(self, genset, coeff_genset, index_comp):
        super().__init__(genset, coeff_genset)
        self._quantum_basis = QuantumDoubleSchubertRing(genset, coeff_genset)
        self._classical_basis = spr.DoubleSchubertRing(genset, coeff_genset)
        self._index_comp = index_comp
        self._n = list(index_comp)
        self._N = [sum(self._n[:i]) for i in range(len(self._n) + 1)]
        parabolic_index = []
        start = 0
        for i in range(len(index_comp)):
            end = start + index_comp[i]
            parabolic_index += list(range(start + 1, end))
            start = end
        self._parabolic_index = parabolic_index
        self._otherlong = Permutation(list(range(self._N[-1], 0, -1)))
        self._longest = self._otherlong * longest_element(parabolic_index)
        self.dtype = type("ParabolicQuantumDoubleSchubertElement", (ParabolicQuantumDoubleSchubertElement,), {"ring": self})

    def _coerce_mul(self, other):
        if isinstance(other, BaseSchubertElement):
            if type(other.ring) is type(self):
                if self.genset == other.ring.genset:
                    return other
            if type(other.ring) is ParabolicQuantumSingleSchubertRing:
                if self.genset == other.ring.genset and self.index_comp == other.ring.index_comp:
                    newbasis = ParabolicQuantumDoubleSchubertRing(self.genset, utils.poly_ring(0), self.index_comp)
                    return newbasis.from_dict(other)
        return None

    def _coerce_add(self, other):
        if isinstance(other, BaseSchubertElement):
            if type(other.ring) is type(self):
                if self.genset == other.ring.genset and self.coeff_genset == other.ring.coeff_genset and self.index_comp == other.ring.index_comp:
                    return other
        return None

    @property
    def symbol_elem_func(self):
        def elem_func(p, k, varl1, varl2):  # noqa: ARG001
            if p == 0 and k >= 0:
                return 1
            if p < 0 or p > k:
                return 0
            return sympy.Add(*[(Symbol(f"e_{p - i}_{k}") if p - i > 0 else 1) * complete_sym_poly(i, k + 1 - p, [-v for v in varl2]) for i in range(p + 1)])

        return elem_func

    def elem_sym_subs(self, kk):
        elems = []
        elem_func = self.elem_sym
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(sympy.Symbol(f"e_{p}_{k}"), elem_func(p, k, self.genset[1:], utils.poly_ring(0)))]
        return dict(elems)

    @property
    def parabolic_index(self):
        return self._parabolic_index

    @property
    def quantum_basis(self):
        return self._quantum_basis

    @property
    def classical_basis(self):
        return self._classical_basis

    def elem_sym(self, p, k, varl1, varl2):
        if p < 0 or p > k:
            return 0
        if p == 0 and k >= 0:
            return 1
        if k <= self._N[1]:
            return elem_sym_poly(p, k, varl1, varl2)
        ret = 0
        j = bisect_left(self._N, k)
        if j < len(self._N) and k == self._N[j]:
            ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * self.elem_sym(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
        ret += self.elem_sym(p, k - 1, varl1, varl2) + (varl1[k - 1] - varl2[k - p]) * self.elem_sym(p - 1, k - 1, varl1, varl2)
        return ret

    @property
    def index_comp(self):
        return self._index_comp

    def process_coeff_dict(self, coeff_dict):
        max_len = max(len(w) for w in coeff_dict)
        parabolic_index = [*self._parabolic_index]
        if max_len > len(self._longest):
            parabolic_index = []
            start = 0
            index_comp = [*self._n, max_len + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                start = end
        return yz.apply_peterson_woodward(coeff_dict, parabolic_index)

    @cache
    def cached_product(self, u, v, basis2):
        initial_dict = {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}
        return self.process_coeff_dict(initial_dict)

    def in_quantum_basis(self, elem):
        result = S.Zero
        for k, v in elem.items():
            result += v * schubpoly_from_elems(k, self.genset, self.coeff_genset, self.quantum_elem_func)
        return result

    def in_classical_basis(self, elem):
        result = S.Zero
        for k, v in elem.items():
            result += v * self.quantum_as_classical_schubpoly(k)
        return result

    @cache
    def classical_in_basis(self, k):
        a = self.classical_basis(k)
        b = self(k)
        if expand(a.as_polynomial() - b.as_polynomial()) == S.Zero:
            return b

        cd = b.as_classical()
        for k2, v in cd.items():
            if k != k2:
                b -= v * self.classical_in_basis(k2)
        return b

    @property
    def classical_elem_func(self):
        basis = spr.DoubleSchubertRing(self.genset, self.coeff_genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis([])
            if p < 0 or p > k:
                return basis.zero
            if k <= self._N[1]:
                return basis(elem_sym_poly(p, k, varl1, varl2))
            ret = basis.zero
            j = bisect_left(self._N, k)
            if j < len(self._N) and k == self._N[j]:
                ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            ret += elem_func(p, k - 1, varl1, varl2) + (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2)
            return ret

        return elem_func

    @property
    def quantum_elem_func(self):
        basis = QuantumDoubleSchubertRing(self.genset, self.coeff_genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis([])
            if p < 0 or p > k:
                return basis(0)
            if k <= self._N[1]:
                return basis(elem_sym_poly(p, k, varl1, varl2))
            ret = basis(0)
            j = bisect_left(self._N, k)
            if j < len(self._N) and k == self._N[j]:
                ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            ret += elem_func(p, k - 1, varl1, varl2) + (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2)
            return ret

        return elem_func

    def printing_term(self, k):
        return PQDSchubPoly(k, self, self.index_comp)

    @cache
    def quantum_as_classical_schubpoly(self, perm):
        k = perm
        if len(k) > len(self._longest):
            parabolic_index = []
            start = 0
            index_comp = [*self._n, len(k) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation(longest_element(parabolic_index))
            longest = otherlong * longpar
        else:
            longest = self._longest
        return schubpoly_from_elems(perm, self.genset, self.coeff_genset, elem_func=self.classical_elem_func, mumu=~longest)

    @cache
    def cached_schubpoly(self, k):
        if len(k) > len(self._longest):
            parabolic_index = []
            start = 0
            index_comp = [*self._n, len(k) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation(longest_element(parabolic_index))
            longest = otherlong * longpar
        else:
            longest = self._longest
        return schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=self.elem_sym, mumu=~longest)

    @cache
    def cached_positive_product(self, u, v, basis2):
        initial_dict = {
            k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset if basis2.coeff_genset else utils.poly_ring(utils.NoneVar)) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()
        }
        return self.process_coeff_dict(initial_dict)

    @property
    def double_mul(self):
        from schubmult.schub_lib.quantum_double import _vars

        def do_double_mul(perm_dict, v, var2=_vars.var2, var3=_vars.var3, q_var=_vars.q_var):
            coeff_dict = yz.schubmult_q_double_fast(perm_dict, v, var2, var3, q_var)
            return self.process_coeff_dict(coeff_dict)

        return do_double_mul

    @property
    def single_mul(self):
        from schubmult.schub_lib.quantum_double import _vars

        def do_single_mul(perm_dict, v, q_var=_vars.q_var):
            coeff_dict = py.schubmult_q_fast(perm_dict, v, q_var)
            return self.process_coeff_dict(coeff_dict)

        return do_single_mul

    @property
    def mult_poly_single(self):
        return py.mult_poly_q

    @property
    def mult_poly_double(self):
        return yz.mult_poly_q_double

    def from_sympy(self, x):
        dct = self.classical_basis.from_sympy(x)
        elem = self.zero
        if not isinstance(dct, BaseSchubertElement):
            return dct
        for k, v in dct.items():
            if elem == self.zero:
                elem = v * self.classical_in_basis(k)
            else:
                elem += v * self.classical_in_basis(k)
        return elem

    def __call__(self, x):
        genset = self.genset
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            perm = Permutation(x)
            if not is_parabolic(perm, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {perm} is not")
            elem = self.from_dict({perm: 1})
        elif isinstance(x, Permutation):
            if not is_parabolic(x, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {x} is not")
            elem = self.from_dict({x: 1})
        elif isinstance(x, ParabolicQuantumDoubleSchubertElement):
            return x
        else:
            elem = self.from_sympy(sympify(x))
        return elem


QSx = QuantumSingleSchubertRing(GeneratingSet("x"))

QuantumDoubleSchubertPolynomial = QuantumDoubleSchubertElement


def make_parabolic_quantum_basis(index_comp, coeff_genset):
    return ParabolicQuantumDoubleSchubertRing(GeneratingSet("x"), coeff_genset, index_comp)


def QPDSx_index(*args):
    def this_QPDSx(x, coeff_genset=GeneratingSet("y")):
        return make_parabolic_quantum_basis(args, utils.poly_ring(coeff_genset) if isinstance(coeff_genset, str) else coeff_genset)(x)

    return this_QPDSx


@cache
def QPDSx(*args):
    return QPDSx_index(*args)


class ParabolicQuantumSingleSchubertRing(ParabolicQuantumDoubleSchubertRing):
    def __init__(self, genset, index_comp):
        super().__init__(genset, utils.poly_ring(utils.NoneVar), index_comp)

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, self.index_comp, "PQSB"))

    def _coerce_mul(self, other):
        if isinstance(other, BaseSchubertElement):
            if type(other.ring) is type(self):
                if self.genset == other.ring.genset:
                    return other
            if type(other.ring) is ParabolicQuantumDoubleSchubertRing:
                if self.genset == other.ring.genset and self.index_comp == other.ring.index_comp:
                    newbasis = ParabolicQuantumDoubleSchubertRing(self.genset, utils.poly_ring(0), self.index_comp)
                    return newbasis.from_dict(other)
        return None

    @property
    def coeff_genset(self):
        return utils.poly_ring(utils.NoneVar)

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            initial_dict = py.schubmult_q_fast({u: S.One}, v)
        else:
            initial_dict = {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}
        return self.process_coeff_dict(initial_dict)

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def __call__(self, x):
        genset = self.genset
        if not genset:
            genset = self.genset
        if isinstance(x, list) or isinstance(x, tuple):
            perm = Permutation(x)
            if not is_parabolic(perm, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {perm} is not")
            elem = self.from_dict({perm: 1})
        elif isinstance(x, Permutation):
            if not is_parabolic(x, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {x} is not")
            elem = self.from_dict({x: 1})
        elif isinstance(x, ParabolicQuantumDoubleSchubertElement):
            return x
        else:
            elem = self.from_sympy(x)
        return elem


def make_single_parabolic_quantum_basis(index_comp):
    return ParabolicQuantumSingleSchubertRing(GeneratingSet("x"), index_comp)


@cache
def QPSx(*args):
    return make_single_parabolic_quantum_basis(args)
