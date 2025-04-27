from functools import cache

import sympy
from symengine import Add, S, expand, sympify
from sympy import CoercionFailed
from sympy.polys.domains import EXRAW
from sympy.polys.domains.compositedomain import CompositeDomain
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.ring import Ring
from sympy.printing.defaults import DefaultPrinting

import schubmult.rings._utils as utils
from schubmult.perm_lib import Permutation
from schubmult.poly_lib.schub_poly import schubpoly_classical_from_elems, schubpoly_from_elems
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from .abstract_schub_poly import AbstractSchubPoly

logger = get_logger(__name__)


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
        if len(self.keys()) == 0:
            return [sympy.sympify(S.Zero)]
        return [
            (self.ring.domain.to_sympy(self[k]) if k == Permutation([]) else sympy.Mul(self.ring.domain.to_sympy(self[k]), self.ring.printing_term(k)))
            for k in sorted(self.keys(), key=lambda kk: (kk.inv, tuple(kk)))
        ]

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
            return self.__add__(new_other)
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
                return self.ring.from_dict(utils._mul_schub_dicts(self, other, self.ring, other.ring))
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
                    return self.from_dict(utils._mul_schub_dicts(elem, other, elem.ring, other.ring))
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


class TensorRing(BaseSchubertRing):
    # tensor ring
    def __eq__(self, other):
        return type(self) is type(other) and self.rings == other.rings

    @property
    def args(self):
        return ()

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
        self.latex_tensor_symbol = "\\otimes"
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
                    dct = utils._tensor_product_of_dicts(dct, dct2)
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
                res += self.from_dict(utils._tensor_product_of_dicts(elem1, self.rings[i].from_sympy(v)))
            elem1 = res
        return elem1

    def __call__(self, x):
        if isinstance(x, tuple):
            return self.from_dict({x: 1})
        return self.from_sympy(x)


# def TensorAlgebra_basis(Basic):
#     pass


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
            if self._key == self.ring.zero_monom:
                return printer._print(S.One)
            for i in range(len(self._key)):
                if self._key[i] != Permutation([]):
                    return printer._print(self.ring.rings[i].printing_term(self._key[i]))
        return self.ring.tensor_symbol.join([printer._print(self.ring.rings[i].printing_term(self._key[i])) for i in range(len(self._key))])

    def _pretty(self, printer):
        if self.ring.is_implicit:
            if self._key == self.ring.zero_monom:
                return printer._print(S.One)
            for i in range(len(self._key)):
                if self._key[i] != Permutation([]):
                    return printer._print(self.ring.rings[i].printing_term(self._key[i]))
        return printer._print_TensorProduct(sympy.Mul(*[self.ring.rings[i].printing_term(self._key[i]) for i in range(len(self._key))]))

    def _latex(self, printer):
        if self.ring.is_implicit:
            if self._key == self.ring.zero_monom:
                return printer._print(S.One)
            for i in range(len(self._key)):
                if self._key[i] != Permutation([]):
                    return printer._print(self.ring.rings[i].printing_term(self._key[i]))
        return printer._print_TensorProduct(sympy.Mul(*[self.ring.rings[i].printing_term(self._key[i]) for i in range(len(self._key))]))


class TensorRingElement(BaseSchubertElement):
    @property
    def free_symbols(self):
        ret = set()
        for k, v in self.items():
            ret.update(v.free_symbols)
            for i in range(len(k)):
                ret.update(self.rings.rings[i](k[i]).free_symbols)
        return ret

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympy.sympify(S.Zero)]
        return [
            self.ring.domain.to_sympy(self[k]) if k == self.ring.zero_monom else sympy.Mul(self.ring.domain.to_sympy(self[k]), self.ring.printing_term(k))
            for k in sorted(self.keys(), key=lambda kkt: [(kk.inv, tuple(kk)) for kk in kkt])
        ]
