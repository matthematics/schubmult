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
            return MixedSchubertElement(other, self)
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
            return MixedSchubertElement(other, -self)
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
            if isinstance(self, MixedSchubertElement):
                if isinstance(other, MixedSchubertElement):
                    return dict.__eq__(self, other)
                return (self - other).almosteq(sympy.S.Zero)
            if isinstance(other, MixedSchubertElement):
                return (self - other).almosteq(sympy.S.Zero)
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


class MixedSchubertElement(BaseSchubertElement, dict):
    def __new__(cls, *elems):
        obj = dict.__new__(cls)
        for elem in elems:
            obj[elem.ring] = elem.ring.add(obj.get(elem.ring, elem.ring.zero), elem)
        if len(obj.keys()) == 1:
            return next(iter(obj.values()))
        return obj

    def __init__(self, *elems):  # noqa: ARG002
        self._ring = BaseSchubertRing([], [])
        self._ring.dtype = type("MixedSchubertElement", (MixedSchubertElement,), {"ring": self})

    def __hash__(self):
        return hash(tuple(self.values()))

    @property
    def ring(self):
        return self._ring

    def __add__(self, other):
        if isinstance(other, MixedSchubertElement):
            return MixedSchubertElement(*list(add_perm_dict(self, other).values()))
        if isinstance(other, BaseSchubertElement):
            elem = MixedSchubertElement(*list(self.values()))
            elem[other.ring] = other.ring.add(elem.get(other.ring, other.ring.zero), other)
            return elem
        elem = MixedSchubertElement(*list(self.values()))

        for ring in self:
            try:
                new_other = ring.domain_new(other)
                new_other = ring.from_dict({Permutation([]): new_other})
                elem[ring] = ring.add(elem.get(ring, ring.zero), new_other)
                if elem[ring] == ring.zero:
                    del elem[ring]
                return MixedSchubertElement(*list(elem.values()))
            except CoercionFailed:
                new_other = ring.from_sympy(other)
                elem[ring] = ring.add(elem.get(ring, ring.zero), new_other)
                if elem[ring] == ring.zero:
                    del elem[ring]
                return MixedSchubertElement(*list(elem.values()))
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, MixedSchubertElement):
            return MixedSchubertElement(*list(add_perm_dict(self, -other).values()))
        if isinstance(other, BaseSchubertElement):
            elem = MixedSchubertElement(*list(self.values()))
            elem[other.ring] = other.ring.sub(elem.get(other.ring, other.ring.zero), other)
            return elem
        elem = MixedSchubertElement(*list(self.values()))
        for ring in self:
            try:
                new_other = ring.domain_new(other)
                new_other = ring.from_dict({Permutation([]): new_other})
                elem[ring] = ring.sub(elem.get(ring, ring.zero), new_other)
                if elem[ring] == ring.zero:
                    del elem[ring]
                return MixedSchubertElement(*list(elem.values()))
            except CoercionFailed:
                new_other = ring.from_sympy(other)
                elem[ring] = ring.sub(elem.get(ring, ring.zero), new_other)
                if elem[ring] == ring.zero:
                    del elem[ring]
                return MixedSchubertElement(*list(elem.values()))
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, BaseSchubertElement):
            elem = MixedSchubertElement(*list(self.values()))
            elem[other.ring] = other.ring.add(other, elem.get(other.ring, other.ring.zero))
            return elem
        elem = MixedSchubertElement(*list(self.values()))
        for ring in self:
            try:
                new_other = ring.domain_new(other)
                new_other = ring.from_dict({Permutation([]): new_other})
                elem[ring] = ring.add(new_other, elem.get(ring, ring.zero))
                if elem[ring] == ring.zero:
                    del elem[ring]
                return MixedSchubertElement(*list(elem.values()))
            except CoercionFailed:
                new_other = ring.from_sympy(other)
                elem[ring] = ring.add(new_other, elem.get(ring, ring.zero))
                if elem[ring] == ring.zero:
                    del elem[ring]
                return MixedSchubertElement(*list(elem.values()))
        return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, BaseSchubertElement):
            elem = MixedSchubertElement(*list(self.values()))
            elem[other.ring] = other.ring.sub(other, elem.get(other.ring, other.ring.zero))
            return elem
        elem = MixedSchubertElement(*list(self.values()))
        for ring in self:
            try:
                new_other = ring.domain_new(other)
                new_other = ring.from_dict({Permutation([]): new_other})
                elem[ring] = ring.sub(new_other, elem.get(ring, ring.zero))
                if elem[ring] == ring.zero:
                    del elem[ring]
                return MixedSchubertElement(*list(elem.values()))
            except CoercionFailed:
                new_other = ring.from_sympy(other)
                elem[ring] = ring.sub(new_other, elem.get(ring, ring.zero))
                if elem[ring] == ring.zero:
                    del elem[ring]
                return MixedSchubertElement(*list(elem.values()))
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, MixedSchubertElement):
            elem = MixedSchubertElement()
            for v1 in self.values():
                for v2 in other.values():
                    elem[v1.ring] = elem.get(v1.ring, v1.ring.zero) + (v1 * v2)
            return elem
        if isinstance(other, BaseSchubertElement):
            return MixedSchubertElement(*[a * other for a in self.values()])
        elem = MixedSchubertElement(*list(self.values()))
        for ring in elem:
            try:
                new_other = ring.domain_new(other)
                new_other = ring.from_dict({Permutation([]): new_other})
                elem[ring] = ring.mul(elem[ring], new_other)
            except CoercionFailed:
                new_other = ring.from_sympy(other)
                elem[ring] = ring.mul(elem[ring], new_other)
        return elem

    def __rmul__(self, other):
        if isinstance(other, BaseSchubertElement):
            return MixedSchubertElement(*[other * a for a in self.values()])
        elem = MixedSchubertElement(*list(self.values()))
        for ring in elem:
            try:
                new_other = ring.domain_new(other)
                new_other = ring.from_dict({Permutation([]): new_other})
                elem[ring] = ring.mul(new_other, elem[ring])
                return elem
            except CoercionFailed:
                new_other = ring.from_sympy(other)
                elem[ring] = ring.mul(new_other, elem[ring])
        return elem

    @property
    def free_symbols(self):
        ret = set()
        for v in self.values():
            ret.update(v.free_symbols)
        return ret

    def as_ordered_terms(self, *_, **__):
        tms = []
        for v in self.values():
            tms += v.as_ordered_terms()
        return tms
