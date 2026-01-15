from sympy import pretty

from schubmult.schub_lib.perm_lib import Permutation
from schubmult.schub_lib.schub_poly import schubpoly_from_elems
from schubmult.symbolic import EXRAW, Add, CoercionFailed, CompositeDomain, DefaultPrinting, DomainElement, Ring, S, expand, sympify, sympify_sympy, sympy_Add, sympy_Mul
from schubmult.utils._mul_utils import _mul_schub_dicts
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

logger = get_logger(__name__)


class BaseSchubertElement(DomainElement, DefaultPrinting, dict):
    _op_priority = 1e200
    precedence = 40

    __sympy__ = True

    def __reduce__(self):
        return (self.__class__, self.items())

    @property
    def is_zero(self):
        return all(v == S.Zero for v in self.values())

    def parent(self):
        return self.ring

    def has_free(self, *args):
        return any(s in args for s in self.free_symbols)

    def eval(self, *args):
        pass

    def mult_poly(self, poly):
        res_dict2 = {}
        for k, v in self.items():
            if self.ring.coeff_genset.label is None:
                dict2 = self.ring.mult_poly_single({k: v}, poly, self.ring.genset)
            else:
                dict2 = self.ring.mult_poly_double({k: v}, poly, self.ring.genset, self.ring.coeff_genset)
            res_dict2 = add_perm_dict(res_dict2, dict2)
        return self.ring.from_dict(res_dict2)

    def in_schubert_schur_basis(self, numvars):
        res = (self.ring @ self.ring).zero
        for perm, v in self.items():
            res += v * self.ring.in_schubert_schur_basis(perm, numvars)
        return res

    def in_SEM_basis(self):
        result = S.Zero
        for k, v in self.items():
            result += sympify(v) * schubpoly_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=self.ring.symbol_elem_func)
        return result

    def _sympystr(self, printer):
        return printer._print(pretty(self, use_unicode=False))

    def _pretty(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _latex(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def as_terms(self):
        if len(self.keys()) == 0:
            return [sympify_sympy(S.Zero)]
        return [((self[k]) if k == self.ring.zero_monom else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == self.ring.zero_monom else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys(), key=lambda kk: (kk.inv, tuple(kk)))]

    def __add__(self, other):
        if isinstance(other, BaseSchubertElement):
            # if self.ring == other.ring:
            return self.ring.add(self, other)
            # return other.__radd__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self.ring.zero_monom: other})
            return self.ring.add(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__add__(new_other)
        except CoercionFailed:
            return other.__radd__(self)

    def __radd__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self.ring.zero_monom: other})
            return self.ring.add(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__add__(self)
        except CoercionFailed:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, BaseSchubertElement):
            return self.ring.sub(self, other)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self.ring.zero_monom: other})
            return self.ring.sub(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__sub__(new_other)
        except CoercionFailed:
            return other.__rsub__(self)

    def __rsub__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self.ring.zero_monom: other})
            return self.ring.sub(other, self)
        except CoercionFailed:
            return NotImplemented

    def __neg__(self):
        return self.ring.neg(self)

    def __pow__(self, val):
        try:
            val = int(val)
        except Exception:
            return NotImplemented
        if val == 0:
            return self.ring.one
        if val < 0:
            return NotImplemented
        if val == 1:
            return self
        return (self ** (val - 1)) * self

    def __mul__(self, other):
        return self.ring.mul(self, other)

    def __rmul__(self, other):
        try:
            return self.ring.rmul(self, other)
        except CoercionFailed:
            return NotImplemented

    def coproduct(self):
        result = self.ring.zero
        for k, v in self.items():
            result += v * self.ring.coproduct_on_basis(k)
        return result

    def as_coefficients_dict(self):
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    def _eval_expand_basic(self, *args, **kwargs):  # noqa: ARG002
        return self.as_polynomial()

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        if not deep:
            return self.ring.from_dict({k: expand(v, **kwargs) for k, v in self.items()})
        return sympify(expand(self.as_polynomial()))

    def as_expr(self):
        return Add(*self.as_terms())

    def as_polynomial(self):
        # print(f"{self=}")
        # try:
        return Add(*[v * self.ring.cached_schubpoly(k) for k, v in self.items()])
        # except SympifyError:
        #     return Add(*[sympify(v) * self.ring.cached_schubpoly(k) for k, v in self.items()])

    def as_classical(self):
        return self.ring.in_classical_basis(self)

    def as_quantum(self):
        return self.ring.in_quantum_basis(self)

    def __eq__(self, other):
        return type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)

    # unify base ring
    def almosteq(self, other):
        if isinstance(other, BaseSchubertElement):
            elem1 = self
            elem2 = other
            if elem1.ring == elem2.ring:
                return (self - other).expand(deep=False).almosteq(S.Zero)
            return elem1.almosteq(elem1.ring.one * elem2)
        return (self - self.ring.from_expr(other)).expand(deep=False) == self.ring.zero


class BaseSchubertRing(Ring, CompositeDomain):
    def __str__(self):
        return self.__class__.__name__

    def __matmul__(self, other):
        from .tensor_ring import TensorRing

        return TensorRing(self, other)
        # return NotImplemented

    def __eq__(self, other):
        return type(self) is type(other) and self.genset == other.genset and self.coeff_genset == other.coeff_genset

    def to_sympy(self, elem):
        return elem.as_expr()

    def __init__(self, genset, coeff_genset, domain=None):
        self._genset = genset
        self._coeff_genset = coeff_genset
        self.symbols = list(genset)
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain
        self.zero_monom = Permutation([])

    def add(self, elem, other):
        res = self.from_dict(add_perm_dict(elem, other))
        return self.from_dict({k: v for k, v in res.items() if expand(v) != S.Zero})

    def sub(self, elem, other):
        res = self.from_dict(add_perm_dict(elem, {k: -v for k, v in other.items()}))
        return self.from_dict({k: v for k, v in res.items() if expand(v) != S.Zero})

    def neg(self, elem):
        return self.from_dict({k: -v for k, v in elem.items()})

    def rmul(self, elem, other):
        try:
            other = self.domain_new(other)
            return self.from_dict({k: v * other for k, v in elem.items()})
        except Exception:
            return self.mul_expr(elem, other)

    def mul(self, elem, other):
        try:
            other = self.domain_new(other)
            # print(f"{other=} {type(other)=}")
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            if isinstance(other, BaseSchubertElement):
                other = self._coerce_mul(other)
                if not other:
                    raise CoercionFailed(f"Could not coerce {other} of type {type(other)} to {type(elem)}")
                return self.from_dict(_mul_schub_dicts(elem, other, elem.ring, other.ring))
            return self.mul_expr(elem, other)
        return NotImplemented

    def to_domain(self):
        return self

    def from_sympy(self, expr):
        return self.from_expr(expr)

    def new(self, x): ...

    def printing_term(self, k): ...

    def coproduct_on_basis(self, k): ...

    def _coerce_mul(self, other): ...

    @property
    def one(self):
        return self.from_dict({Permutation([]): S.One})

    def is_elem_mul_type(self, elem): ...

    def elem_mul(self, ring_elem, elem): ...

    def _coerce_add(self, other): ...

    def from_dict(self, element, orig_domain=None):
        # print(f"{element=}")
        domain_new = self.domain_new
        poly = self.zero

        for monom, coeff in element.items():
            coeff = domain_new(coeff, orig_domain)
            if coeff != self.domain.zero:
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

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        try:
            if isinstance(element, BaseSchubertElement) or isinstance(element, DomainElement):
                raise CoercionFailed("Not a domain element")
            if not any(x in self.genset for x in sympify_sympy(element).free_symbols):
                return sympify(element)
            raise CoercionFailed(f"{element} contains an element of the set of generators")
        except Exception:
            raise CoercionFailed(f"Could not coerce {element} of type {type(element)} to {self.__class__.__name__}")

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

    def from_expr(self, x):
        return self.mul_expr(self.one, x)

    def mul_expr(self, elem, x): ...

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
