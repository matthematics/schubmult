from functools import cache, cached_property

import sympy
from symengine import Add, S, expand, sympify
from sympy import CoercionFailed
from sympy.polys.domains import EXRAW
from sympy.polys.domains.compositedomain import CompositeDomain
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.ring import Ring
from sympy.printing.defaults import DefaultPrinting
from sympy.printing.str import StrPrinter

import schubmult.rings._quantum_schubert_polynomial_ring as qsr
import schubmult.rings._utils as utils
import schubmult.schub_lib.double as yz
import schubmult.schub_lib.schub_lib as schub_lib
import schubmult.schub_lib.single as py
from schubmult.perm_lib import Permutation, uncode
from schubmult.poly_lib.poly_lib import elem_sym_poly, xreplace_genvars
from schubmult.poly_lib.schub_poly import schubpoly_classical_from_elems, schubpoly_from_elems
from schubmult.poly_lib.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

_def_printer = StrPrinter({"order": "none"})

logger = get_logger(__name__)


class NotEnoughGeneratorsError(ValueError):
    pass


def _mul_schub_dicts(dict1, dict2, basis1, basis2, best_effort_positive=True):
    # print(f"{dict1.items()=} {dict2.items()=}")
    this_dict = {}
    for k, v in dict2.items():
        for kd, vd in dict1.items():
            did_positive = False
            to_mul = v * vd
            if best_effort_positive:
                try:
                    this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis1.cached_positive_product(kd, k, basis2).items()})
                    did_positive = True
                except Exception:
                    did_positive = False
            if not did_positive:
                this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis1.cached_product(kd, k, basis2).items()})
    return this_dict


class BaseSchubertElement(DomainElement, DefaultPrinting, dict):
    _op_priority = 1e200
    precedence = 40

    def parent(self):
        return self.ring

    def has_free(self, *args):
        return any(s in args for s in self.free_symbols)

    def mult_poly(self, poly):
        res_dict2 = {}
        for k, v in self.items():
            if self.ring.coeff_genset.label is None:
                dict2 = self.ring.mult_poly_single({k: v}, poly, self.ring.genset)
            else:
                dict2 = self.ring.mult_poly_double({k: v}, poly, self.ring.genset, self.ring.coeff_genset)
            res_dict2 = add_perm_dict(res_dict2, dict2)
        return self.ring.from_dict(res_dict2)

    def in_SEM_basis(self):
        result = sympy.S.Zero
        for k, v in self.items():
            result += sympy.sympify(v) * schubpoly_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=self.ring.symbol_elem_func)
        return result

    def in_CEM_basis(self):
        result = sympy.S.Zero
        for k, v in self.items():
            result += sympy.sympify(v) * schubpoly_classical_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=self.ring.symbol_elem_func)
        return result

    def _sympystr(self, printer):
        if len(self.keys()) == 0:
            return printer._print(sympy.S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(self)

    def _pretty(self, printer):
        if len(self.keys()) == 0:
            return printer._print(sympy.S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(self)

    def _latex(self, printer):
        if len(self.keys()) == 0:
            return printer._print(sympy.S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(self)

    def as_terms(self):
        if len(self.keys()) == 0:
            return [sympy.sympify(S.Zero)]
        return [(self.ring.domain.to_sympy(self[k]) if k == Permutation([]) else sympy.Mul(self.ring.domain.to_sympy(self[k]), self.ring.printing_term(k))) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        return self.as_terms()

    def __add__(self, other):
        if isinstance(other, BaseSchubertElement):
            if self.ring == other.ring:
                return self.ring.add(self, other)
            return other.__radd__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({Permutation([]): other})
            return self.ring.add(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__add_(new_other)
        except CoercionFailed:
            return other.__radd__(self)

    def __sympy__(self):
        return self.ring.to_sympy(self)

    def __radd__(self, other):
        if isinstance(other, BaseSchubertElement):
            new_ring = TensorRing(other.ring, self.ring, tensor_symbol=None)
            return new_ring.add(new_ring.from_comp_ring(other), new_ring.from_comp_ring(self))
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({Permutation([]): other})
            return self.ring.add(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__add_(self)
        except CoercionFailed:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, BaseSchubertElement):
            if self.ring == other.ring:
                return self.ring.sub(self, other)
            return other.__rsub__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({Permutation([]): other})
            return self.ring.sub(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__sub__(new_other)
        except CoercionFailed:
            return other.__rsub__(self)

    def __rsub__(self, other):
        if isinstance(other, BaseSchubertElement):
            new_ring = TensorRing(other.ring, self.ring, tensor_symbol=None)
            return new_ring.sub(new_ring.from_comp_ring(other), new_ring.from_comp_ring(self))
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({Permutation([]): other})
            return self.ring.sub(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__sub_(self)
        except CoercionFailed:
            return NotImplemented

    def __neg__(self):
        return self.ring.neg(self)

    def __mul__(self, other):
        if isinstance(other, BaseSchubertElement):
            if isinstance(other.ring, type(self.ring)):
                return self.ring.from_dict(_mul_schub_dicts(self, other, self.ring, other.ring))
            if isinstance(other.ring, TensorRing) and other.ring.is_implicit:
                terms = [v * self * other.ring.combine(other.ring(k)) for k, v in other.items()]
                if self.ring in other.ring.rings:
                    return other.ring.sum([other.from_comp_ring(t) for t in terms])
                new_ring = TensorRing(self.ring, *other.ring.rings, tensor_symbol=None)
                return new_ring.sum([new_ring.from_comp_ring(t) for t in terms])
            new_other = self.ring._coerce_mul(other)
            if new_other:
                return self.ring.mul(self, new_other)
            return other.__rmul__(self)
        try:
            other = self.ring.domain_new(other)
            return self.ring.from_dict({k: sympify(other) * v for k, v in self.items()})
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__mul__(new_other)
        except CoercionFailed:
            return other.__rmul__(self)

    def __rmul__(self, other):
        if isinstance(other, BaseSchubertElement):
            if isinstance(other.ring, TensorRing) and other.ring.is_implicit:
                terms = [v * other.ring.combine(other.ring(k)) * self for k, v in other.items()]
                if self.ring in other.ring.rings:
                    return other.ring.sum([other.ring.from_comp_ring(t) for t in terms])
                new_ring = TensorRing(self.ring, *other.ring.rings, tensor_symbol=None)
                return new_ring.sum([new_ring.from_comp_ring(t) for t in terms])
            new_other = self.ring._coerce_mul(other)
            if new_other:
                return self.ring.mul(new_other, self)
        try:
            other = self.ring.domain_new(other)
            return self.ring.from_dict({k: sympify(other) * v for k, v in self.items()})
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__mul__(self)
        except CoercionFailed:
            return NotImplemented

    def as_coefficients_dict(self):
        return sympy.Dict({self.ring.printing_term(k, self.ring): sympy.sympify(v) for k, v in self.items()})

    def _eval_expand_basic(self, *args, **kwargs):  # noqa: ARG002
        return self.as_polynomial()

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        if not deep:
            return self.ring.from_dict({k: expand(v) for k, v in self.items()})
        return sympy.sympify(expand(sympify(self.as_polynomial())))

    def as_expr(self):
        return sympy.Add(*self.as_terms())

    def as_polynomial(self):
        return sympy.sympify(Add(*[v * self.ring.cached_schubpoly(k) for k, v in self.items()]))

    def as_classical(self):
        return self.ring.in_classical_basis(self)

    def as_quantum(self):
        return self.ring.in_quantum_basis(self)

    def __eq__(self, other):
        return (type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)) or self.almosteq(other)

    # unify base ring
    def almosteq(self, other):
        if isinstance(other, BaseSchubertElement):
            elem1 = self
            elem2 = other
            if isinstance(self, TensorRingElement):
                elem1 = elem1.ring.combine(elem1)
            if isinstance(other, TensorRingElement):
                elem2 = elem2.ring.combine(elem2)
            if elem1.ring == elem2.ring:
                return dict.__eq__(elem1, elem2)
            return elem1 == elem1.ring.one * elem2
        return self == self.ring.from_sympy(other)


class DoubleSchubertElement(BaseSchubertElement):
    """Algebra with sympy coefficients
    and a dict basis
    """

    def divdiff(self, i):
        return self.ring.from_dict({k.swap(i - 1, i): v for k, v in self.items() if i - 1 in k.descents()})

    def simpleref(self, i):
        return self + self.divdiff(i).mult_poly(self.ring.genset[i + 1] - self.ring.genset[i])

    def act(self, perm):
        perm = Permutation(perm)
        dset = perm.descents()
        if len(dset) == 0:
            return self
        i = next(iter(dset))
        return self.simpleref(i + 1).act(perm.swap(i, i + 1))

    def max_index(self):
        return max([max([0, *list(k.descents(zero_indexed=False))]) for k in self.keys()])

    def subs(self, old, new):
        result = 0
        if self.ring.genset.index(old) != -1:
            result = 0
            index = self.ring.genset.index(old)
            mindex = self.max_index()
            if mindex < index:
                return self
            perm = Permutation([]).swap(index - 1, mindex)  # index to max index + 1
            transf = self.act(perm)
            for k, v in transf.items():
                perm = k
                coeff_gens = self.ring.coeff_genset
                L = schub_lib.pull_out_var(mindex + 1, perm)
                for index_list, new_perm in L:
                    result += self.ring.from_dict({new_perm: v}).mult_poly(sympy.prod([(new - coeff_gens[index2]) for index2 in index_list]))
            return result

        for k, v in self.items():
            if self.ring.coeff_genset.label is None:
                add_dict = {k: v.subs(old, new)}
            else:
                coeff_genset = self.ring.coeff_genset
                if coeff_genset.index(old) != -1:
                    genset_list = [coeff_genset[i] for i in range(len(coeff_genset))]
                    genset_list[coeff_genset.index(old)] = 0
                    custom_genset = CustomGeneratingSet(genset_list)
                    new_add_dict = {k2: sympify(v2).subs(old, new) for k2, v2 in yz.schubmult_double({(): v}, k, custom_genset, coeff_genset).items()}  # remove the variable
                    add_dict = {}
                    for k3, v3 in new_add_dict.items():
                        to_add_dict = yz.schubmult_double({(): v3}, k3, coeff_genset, custom_genset)
                        add_dict = add_perm_dict(add_dict, to_add_dict)
                else:
                    add_dict = {k: sympify(v).subs(old, new)}
            for k5, v5 in add_dict.items():
                if any(self.ring.genset.index(s) != -1 for s in sympify(v5).free_symbols):
                    result += self.ring.from_dict({k5: 1}).mult_poly(v5)
                else:
                    result += self.ring.from_dict({k5: v5})
        return result

    @property
    def free_symbols(self):
        ret = set()
        for k, v in self.items():
            ret.update(v.free_symbols)
            perm = k
            if len(perm.descents()) > 0:
                ret.update([self.ring.genset[i] for i in range(1, max(perm.descents()) + 2)])
            if self.ring.coeff_genset:
                genset2 = self.ring.coeff_genset
                perm2 = ~perm
                if len(perm2.descents()) > 0:
                    ret.update([genset2[i] for i in range(1, max(perm2.descents()) + 2)])
        return ret

    def pull_out_gen(self, gen):
        ind = self.genset.index(gen)
        gens2 = MaskedGeneratingSet(self.genset, [ind])
        gens2.set_label(f"({self.ring.genset.label}\\{gen})")
        new_basis = DoubleSchubertRing(gens2)
        ret = new_basis(0)
        for (perm, cv), val in self.items():
            L = schub_lib.pull_out_var(ind, perm)
            for index_list, new_perm in L:
                toadd = S.One
                for index2 in index_list:
                    toadd *= gen - utils.poly_ring(cv)[index2]
                ret += toadd * val * new_basis(new_perm, cv)
        return ret

    def coproduct(self, *indices, alt_coeff_genset=None, on_coeff_gens=False, gname1=None, gname2=None):
        result_dict = {}
        genset = self.ring.genset
        if on_coeff_gens:
            genset = self.ring.coeff_genset
        if gname1 is None:
            gname1 = f"{genset.label}_A"
        if gname2 is None:
            gname2 = f"{genset.label}_B"
        gens2 = MaskedGeneratingSet(genset, indices)
        gens1 = gens2.complement()
        gens1.set_label(gname1)
        gens2.set_label(gname2)
        for k, v in self.items():
            key = k
            if isinstance(self.ring, SingleSchubertRing) and not alt_coeff_genset:
                coprod_dict = py.schub_coprod_py(key, indices)
            else:
                if on_coeff_gens:
                    coprod_dict = yz.schub_coprod_double(~key, indices, self.ring.genset, alt_coeff_genset if alt_coeff_genset else self.ring.genset)
                else:
                    coprod_dict = yz.schub_coprod_double(key, indices, self.ring.coeff_genset, alt_coeff_genset if alt_coeff_genset else self.ring.coeff_genset)
            if on_coeff_gens:
                result_dict = add_perm_dict(result_dict, {(~k1, ~k2): v * v2 for (k1, k2), v2 in coprod_dict.items()})
            else:
                result_dict = add_perm_dict(result_dict, {k: v * v2 for k, v2 in coprod_dict.items()})
        if on_coeff_gens:
            basis = TensorRing(
                DoubleSchubertRing(self.ring.genset, gens1),
                DoubleSchubertRing(alt_coeff_genset if alt_coeff_genset else self.ring.genset, gens2),
            )
        else:
            basis = TensorRing(
                DoubleSchubertRing(gens1, self.ring.coeff_genset),
                DoubleSchubertRing(gens2, alt_coeff_genset if alt_coeff_genset else self.ring.coeff_genset),
            )
        return basis.from_dict(result_dict)

    @cached_property
    def max_gens(self):
        return max([max(k.descents()) for k in self.keys()])


_pretty_schub_char = "𝔖"  # noqa: RUF001


# Atomic Schubert polynomial
class AbstractSchubPoly(sympy.Expr):
    is_Atom = True

    def __new__(cls, k, basis):
        obj = sympy.Expr.__new__(cls, k, basis)
        obj._key = k
        obj._genset = basis.genset
        obj._basis = basis
        obj._perm = k
        return obj

    def _eval_expand_basic(self, *args, **kwargs):  # noqa: ARG002
        return self.as_polynomial()

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        if not deep:
            return self.ring.from_dict({k: expand(v) for k, v in self.items()})
        return sympy.sympify(expand(sympify(self.as_polynomial())))

    def as_polynomial(self):
        return self.ring.cached_schubpoly(self._key)

    @property
    def ring(self):
        return self._basis

    @property
    def genset(self):
        return self._genset

    @property
    def coeff_genset(self):
        return self.ring.coeff_genset

    @property
    def args(self):
        return (sympy.Tuple(*self._key), self._basis)


class DSchubPoly(AbstractSchubPoly):
    is_Atom = True

    def __new__(cls, k, basis):
        return DSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        return AbstractSchubPoly.__new__(_class, k, basis)

    def _sympystr(self, printer):
        key = self._key
        if self._key == Permutation([]):
            return printer.doprint(1)
        if self.ring.coeff_genset.label is None:
            return printer.doprint(f"S{self.ring.genset.label}({printer.doprint(key)})")
        return printer.doprint(f"DS{self.ring.genset.label}({printer.doprint(key)}, {self.ring.coeff_genset.label})")

    def _pretty(self, printer):
        key = self._key
        gl = self.ring.genset.label
        if key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(int("".join([str(i) for i in key])))
        if self.ring.coeff_genset.label is None:
            return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(gl)))
        return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(f"{self.ring.genset.label}; {self.ring.coeff_genset.label}")))

    def _latex(self, printer):
        key = self._key
        gl = self.ring.genset.label
        if key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(key)
        if self.ring.coeff_genset.label is None:
            return printer._print_Function(sympy.Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(gl)))
        return printer._print_Function(sympy.Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(f"{self.ring.genset.label}; {self.ring.coeff_genset.label}")))

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, basis):
        return DSchubPoly.__xnew__(_class, k, basis)


class BaseSchubertRing(Ring, CompositeDomain):

    def __str__(self):
        return self.__class__.__name__

    def __eq__(self, other):
        return type(self) is type(other) and self.genset == other.genset and self.coeff_genset == other.coeff_genset

    def to_sympy(self, elem):
        return elem.as_expr()

    def __init__(self, genset, coeff_genset):
        self._genset = genset
        self._coeff_genset = coeff_genset
        self.symbols = list(genset)
        self.domain = EXRAW
        self.dom = self.domain
        self.zero_monom = Permutation([])

    def add(self, elem, other):
        return self.from_dict(add_perm_dict(elem, other))

    def sub(self, elem, other):
        return self.from_dict(add_perm_dict(elem, {k: -v for k, v in other.items()}))

    def neg(self, elem):
        return self.from_dict({k: -v for k, v in elem.items()})

    def mul(self, elem, other):
        if self.of_type(elem):
            if isinstance(other, BaseSchubertElement):
                if isinstance(other.ring, type(self)):
                    return self.from_dict(_mul_schub_dicts(elem, other, elem.ring, other.ring))
        try:
            other = self.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        return self.from_dict({k: other * v for k, v in elem.items()})

    def to_domain(self):
        return self

    def new(self, x): ...

    def printing_term(self, k): ...

    def _coerce_mul(self, other): ...

    @property
    def one(self):
        return self.from_dict({Permutation([]): S.One})

    def _coerce_add(self, other): ...

    def from_dict(self, element, orig_domain=None):
        domain_new = self.domain_new
        poly = self.zero

        for monom, coeff in element.items():
            coeff = domain_new(coeff, orig_domain)
            if expand(coeff) != S.Zero:
                poly[monom] = coeff
        return poly

    @property
    def zero(self):
        return self.dtype()

    @property
    def elem_sym(self): ...

    @property
    def symbol_elem_func(self): ...

    def elem_sym_subs(self, kk): ...

    def domain_new(self, element, orig_domain=None):
        if not sympy.sympify(element).has_free(*self.symbols):
            return self.domain.convert(element, orig_domain)
        raise CoercionFailed(f"{element} contains an element of the set of generators")

    @property
    def genset(self):
        return self._genset

    @property
    def coeff_genset(self):
        return self._coeff_genset

    def in_quantum_basis(self, elem): ...

    def in_classical_basis(self, elem): ...

    def quantum_schubpoly(self, perm): ...

    def cached_product(self, u, v, basis2): ...

    def cached_positive_product(self, u, v, basis2): ...

    @property
    def double_mul(self): ...

    @property
    def single_mul(self): ...

    @property
    def mult_poly_single(self): ...

    @property
    def mult_poly_double(self): ...

    @property
    def quantum_elem_func(self): ...

    def cached_schubpoly(self, k): ...


class DoubleSchubertRing(BaseSchubertRing):
    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "DBS"))

    def __init__(self, genset, coeff_genset):
        super().__init__(genset, coeff_genset)
        self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

    def __str__(self):
        return f"Double Schubert polynomial ring in {self.genset.label} and {self.coeff_genset.label}"

    def printing_term(self, k):
        return DSchubPoly(k, self)

    def _coerce_mul(self, other):
        if isinstance(other, BaseSchubertElement):
            if isinstance(other.ring, qsr.QuantumDoubleSchubertRing):
                return other.as_classical()
        return None

    def _coerce_add(self, other):  # noqa: ARG002
        # if isinstance(other, BaseSchubertElement):
        #     if isinstance(other.ring, DoubleSchubertRing) and other.ring.genset == self.genset:
        #         return self.one * other
        return None

    @property
    def elem_sym(self):
        genset = self.genset

        class esf(sympy.Function):
            @classmethod
            def eval(cls, *x):
                pass

            def _eval_expand_func(self, **_):
                if len(self.args) == 2:
                    return elem_sym_poly(int(self.args[0]), int(self.args[1]), genset[1:], utils.poly_ring(0))
                return elem_sym_poly(int(self.args[0]), int(self.args[1]), genset[1:], self.args[2:])

        return esf

    @property
    def symbol_elem_func(self):
        def elem_func(p, k, varl1, varl2):  # noqa: ARG001
            if p == 0 and k >= 0:
                return 1
            if p < 0 or p > k:
                return 0
            if self.coeff_genset.label:
                return self.elem_sym(p, k, *varl2)
            return self.elem_sym(p, k)

        return elem_func

    def elem_sym_subs(self, kk):
        elems = []
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(sympy.Symbol(f"e_{p}_{k}"), elem_sym_poly(p, k, self.genset[1:], utils.poly_ring(0)))]
        return dict(elems)

    def in_quantum_basis(self, elem):
        result = S.Zero
        for k, v in elem.items():
            result += v * self.quantum_schubpoly(k)
        return result

    def in_classical_basis(self, elem):
        return elem

    @cache
    def quantum_schubpoly(self, perm):
        return schubpoly_classical_from_elems(perm, self.genset, self.coeff_genset, self.quantum_elem_func)

    @cache
    def cached_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset if basis2.coeff_genset else utils.poly_ring(0)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset if basis2.coeff_genset else utils.poly_ring(0)) for k, x in yz.schubmult_generic_partial_posify(u, v).items()}

    @property
    def double_mul(self):
        return yz.schubmult_double

    @property
    def single_mul(self):
        return py.schubmult_py

    @property
    def mult_poly_single(self):
        return py.mult_poly_py

    @property
    def mult_poly_double(self):
        return yz.mult_poly_double

    @property
    def quantum_elem_func(self):
        basis = qsr.QuantumDoubleSchubertRing(self.genset, self.coeff_genset)

        def elem_func(p, k, varl1, varl2, xstart=0, ystart=0):
            if p > k:
                return basis(0)
            if p == 0:
                return basis([])
            if p == 1:
                res = basis(varl1[xstart] - varl2[ystart])
                for i in range(1, k):
                    res += basis(varl1[xstart + i] - varl2[ystart + i])
                return res
            if p == k:
                res = basis((varl1[xstart] - varl2[ystart]) * (varl1[xstart + 1] - varl2[ystart]))
                for i in range(2, k):
                    res *= basis(varl1[i + xstart] - varl2[ystart])
                return res
            mid = k // 2
            xsm = xstart + mid
            ysm = ystart + mid
            kmm = k - mid
            res = elem_func(p, mid, varl1, varl2, xstart, ystart) + elem_func(
                p,
                kmm,
                varl1,
                varl2,
                xsm,
                ysm,
            )
            for p2 in range(max(1, p - kmm), min(p, mid + 1)):
                res += elem_func(p2, mid, varl1, varl2, xstart, ystart) * elem_func(
                    p - p2,
                    kmm,
                    varl1,
                    varl2,
                    xsm,
                    ysm - p2,
                )
            return res

        return elem_func

    def monomial_schub(self, monom):
        monom = [*monom]
        while len(monom) > 0 and monom[-1] == 0:
            monom.pop()
        return self._monomial_schub_cache(tuple(monom))

    @cache
    def _monomial_schub_cache(self, monom):
        srt_perm = Permutation.sorting_perm([-i for i in monom])
        schub_perm = uncode(sorted(monom, reverse=True))
        return self.from_dict({(schub_perm, utils.NoneVar): S.One}).act(srt_perm)

    @cache
    def cached_schubpoly(self, k):
        return schubpoly_classical_from_elems(k, self.genset, self.coeff_genset, elem_func=elem_sym_poly)

    def from_sympy(self, x):
        x = sympify(x)
        result = yz.mult_poly_double({Permutation([]): 1}, x, self.genset, self.coeff_genset)
        return self.from_dict(result)

    def new(self, x):
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({p_x: 1})
        elif isinstance(x, Permutation):
            if max([0, *list(x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({x: 1})

        elif isinstance(x, DoubleSchubertElement):
            if x.is_Add or x.is_Mul:
                return x.doit()
            if x.genset == genset:
                return x
            raise ValueError("Different generating set")
        # poly
        # elif isinstance(x, sympy.Poly):
        #     # eject generators we don't want
        #     if not x.has_only_gens():
        #         x, _ = sympy.poly_from_expr(x.as_expr())
        #     new_gens = [g for g in x.gens if self.genset.index(g) != -1]
        #     # end_gens = [g for g in x.gens if self.genset.index(g) == -1]
        #     if len(new_gens) == 0:
        #         # # # logger.debug((f"Didn't find any gens in {x=}")
        #         return self(x.as_expr())
        #     new_gens.sort(key=lambda g: self.genset.index(g))
        #     # expand_gens = [self.genset[i] for i in range(self.genset.index(new_gens[-1])+1)] + end_gens
        #     x = sympy.poly(x, gens=tuple(new_gens))
        #     dct = x.as_dict()
        #     result = 0
        #     for monom, coeff in dct.items():
        #         srt_perm = Permutation.sorting_perm([-i for i in monom])
        #         schub_perm = uncode(sorted(monom, reverse=True))
        #         result += self.from_dict({(schub_perm, utils.NoneVar): coeff}).act(srt_perm)
        #     return result
        else:
            elem = self.from_sympy(x)
        return elem


def DSx(x, genset=GeneratingSet("y")):
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    return DoubleSchubertRing(GeneratingSet("x"), genset)(x)


"""DSx: Double Schubert polynomial generator
DSx is an alias for a DoubleSchubertAlgebraElement_basis object with
GeneratingSet being variables with name x_i for i an integer up to 99.
It is a callable object, and the signature is

DSx(x, cv=None, genset=None)

x is either a tuple, a list, a schubmult.Permutation, or a sympy
or symengine object that you are trying to express in terms of
double Schubert polynomials. cv is a string that is the name of
the base GeneratingSet for the coefficient variable (defaults to
"y"), and genset is the "x" variable generating set by default,
but can be subsituted with a custom GeneratingSet_base object.
"""


class SingleSchubertRing(DoubleSchubertRing):
    def __init__(self, genset):
        super().__init__(genset, utils.poly_ring(utils.NoneVar))

    def __str__(self):
        return f"Schubert polynomial ring in {self.genset.label}"

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "SBS"))

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
        if type(other.ring) is DoubleSchubertRing:
            if self.genset == other.ring.genset:
                return other
        if isinstance(other.ring, qsr.QuantumDoubleSchubertRing):
            return self._coerce_mul(other.as_classical())
        return None

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            return py.schubmult_py({u: S.One}, v)
        return {k: xreplace_genvars(x, utils.poly_ring(0), basis2.coeff_genset if basis2.coeff_genset else utils.poly_ring(0)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def from_sympy(self, x):
        x = sympify(x)
        result = py.mult_poly_py({Permutation([]): 1}, x, self.genset)
        return self.from_dict(result)

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
        elif isinstance(x, DoubleSchubertElement):
            if x.is_Add or x.is_Mul:
                return x
            if x.genset == genset:
                elem = DoubleSchubertElement(x, self)  # , self)
            else:
                return self(x.expand())
        else:
            elem = self.from_sympy(x)
        return elem


Sx = SingleSchubertRing(GeneratingSet("x"))


def _tensor_product_of_dicts(d1, d2):
    ret_dict = {}
    for k1, v1 in d1.items():
        this_dict = {}
        for k2, v2 in d2.items():
            if isinstance(k1, tuple):
                this_dict[(*k1, k2)] = v1 * v2
            else:
                this_dict[(k1, k2)] = v1 * v2
        ret_dict = add_perm_dict(ret_dict, this_dict)
    return ret_dict


class TensorBasisElement(AbstractSchubPoly):
    is_commutative = False

    def __new__(cls, k, basis):
        return TensorBasisElement.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj = AbstractSchubPoly.__new__(_class, k, basis)
        if not basis.is_implicit:
            obj.precedence = 50
        else:
            obj.precedence = 1000
        return obj

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, basis):
        return TensorBasisElement.__xnew__(_class, k, basis)

    def _sympystr(self, printer):
        if self.ring.is_implicit:
            # print(f"{self._key=}")
            if self._key == self.ring.zero_monom:
                return printer._print(S.One)
            for i in range(len(self._key)):
                if self._key[i] != Permutation([]):
                    # print(f"Found {self._key[i]=} {i=} {self.ring.rings[i]=}")
                    return printer._print(self.ring.rings[i].printing_term(self._key[i]))
            # print(f"This shouldn't happen {self._key=}")
            # return printer._print_Mul(sympy.Mul(self._elem1, self._elem2))
        return self.ring.tensor_symbol.join([printer._print(self.ring.rings[i].printing_term(self._key[i])) for i in range(len(self._key))])


class TensorRingElement(BaseSchubertElement):
    # tensor ring
    # def __new__(cls, _dict, basis):
    #     return TensorAlgebraElement.__xnew_cached__(cls, Dict(_dict), basis)
    # def as_terms(self):
    #     if len(self.keys()) == 0:
    #         return [sympy.sympify(S.Zero)]
    #     return [
    #         (self.ring.domain.to_sympy(self[k]) if k == self.ring.zero_monom else sympy.Mul(self.ring.domain.to_sympy(self[k]), *[self.ring.rings[i].printing_term(k[i]) for i in range(len(k))]))
    #         for k in self.keys()
    #     ]

    @property
    def free_symbols(self):
        ret = set()
        for k, v in self.items():
            ret.update(v.free_symbols)
            for i in range(len(k)):
                ret.update(self.rings.rings[i](k[i]).free_symbols)
        return ret

    # @cahe
    # def _sympystr(self, printer):
    #     ret_list = [sympy.Mul(v, TensorBasisElement(k[0], k[1], self.ring)) for k, v in self.items()]
    #     return printer.doprint(Add(*ret_list))

    # def _sympystr(self, printer):
    #     return self.cached_sympystr(printer)

    # def expand(self, **_):
    #     return sympify(
    #         Add(
    #             *[v * sympify(self.ring.ring1.from_dict({k[0]: 1}).expand()) * sympify(self.ring.ring2.from_dict({k[1]: 1}).expand()) for k, v in self.items()],
    #         ),
    #     )


class TensorRing(BaseSchubertRing):
    # tensor ring
    def __eq__(self, other):
        return type(self) is type(other) and self.rings == other.rings

    def __init__(self, *rings, tensor_symbol=" # "):
        self._rings = list(rings)
        new_rings = []
        unpacked = True
        while unpacked:
            unpacked = False
            for r in self._rings:
                if isinstance(r, TensorRing):
                    unpacked = True
                    for r2 in r.rings:
                        if r2 not in new_rings:
                            new_rings += [r2]
                else:
                    if r not in new_rings:
                        new_rings += [r]
            self._rings = new_rings
        self._rings = tuple(new_rings)
        genset = set()
        for r in self._rings:
            genset.update(set(r.genset))
        genset = tuple(genset)
        coeff_genset = set()
        for r in self._rings:
            coeff_genset.update(set(r.coeff_genset))
        coeff_genset = tuple(coeff_genset)
        super().__init__(list(genset), list(coeff_genset))
        # self._rings = tuple(sorted(set(rings)))
        self.zero_monom = tuple([self.rings[i].zero_monom for i in range(len(self.rings))])
        self.tensor_symbol = tensor_symbol
        self.dtype = type("TensorRingElement", (TensorRingElement,), {"ring": self})

    def __hash__(self):
        return hash(self._rings)

    @property
    def rings(self):
        return self._rings

    def mul(self, elem1, elem2):
        ret_dict = {}
        for k1, v1 in elem1.items():
            for k2, v2 in elem2.items():
                dct = self.rings[0].from_dict({k1[0]: v1 * v2}) * self.rings[0].from_dict({k2[0]: 1})
                for i in range(1, len(self.rings)):
                    dct2 = self.rings[i].from_dict({k1[i]: 1}) * self.rings[i].from_dict({k2[i]: 1})
                    dct = _tensor_product_of_dicts(dct, dct2)
                ret_dict = add_perm_dict(ret_dict, dct)
        elem = self.from_dict(ret_dict)
        if self.is_implicit:
            return self.combine(elem)
        return elem

    def combine(self, elem):
        if self.is_implicit:
            if len(elem.keys()) == 1:
                k = next(iter(elem.keys()))
                if k == self.zero_monom:
                    return self.zero
                for i in range(len(k)):
                    if k[i] != self.rings[i].zero_monom:
                        return elem[k] * self.rings[i](k[i])
                return self.domain_new(next(iter(elem.values())))
        res = self.rings[0].zero
        for k, v in elem.items():
            to_add = v
            for i in range(len(self.rings)):
                to_add *= self.rings[i](k[i])
            res += to_add
        return res

    @property
    def is_implicit(self):
        return self.tensor_symbol is None

    def _coerce_add(self, x):  # noqa: ARG002
        return None

    def _coerce_mul(self, x):  # noqa: ARG002
        return None

    @property
    def one(self):
        return self.from_dict({self.zero_monom: S.One})

    @cache
    def cached_schubpoly(self, k):
        return sympy.Mul(*[self.rings[i].cached_schubpoly(k[i]) for i in range(len(self.rings))])

    def printing_term(self, k):
        # print(f"beffle {k=}")
        return TensorBasisElement(k, self)

    def from_comp_ring(self, t):
        dct = {}
        for k, v in t.items():
            new_k = list(self.zero_monom)
            if isinstance(t.ring, TensorRing):
                for i in range(len(t.ring.rings)):
                    new_k[self.rings.index(t.ring.rings[i])] = k[i]
            else:
                new_k[self.rings.index(t.ring)] = k
            dct[tuple(new_k)] = v
        return self.from_dict(dct)

    def from_sympy(self, x):
        elem1 = self.rings[0].from_sympy(x)
        for i in range(1, len(self.rings)):
            res = self.zero
            for k, v in elem1.items():
                res += self.from_dict(_tensor_product_of_dicts(elem1, self.rings[i].from_sympy(v)))
            elem1 = res
        return elem1

    def __call__(self, x):
        if isinstance(x, tuple):
            return self.from_dict({x: 1})
        return self.from_sympy(x)


# def TensorAlgebra_basis(Basic):
#     pass
