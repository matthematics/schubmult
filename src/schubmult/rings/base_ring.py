from schubmult.symbolic import EXRAW, Add, CoercionFailed, CompositeDomain, DefaultPrinting, DomainElement, Ring, S, expand, sympify, sympify_sympy, sympy_Add, sympy_Mul
from schubmult.utils.perm_utils import add_perm_dict


class BaseRingElement(DomainElement, DefaultPrinting, dict):
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

    def _sympystr(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):
            return printer._print_Add(sympy_Add(*self.as_ordered_terms()), order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _pretty(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _latex(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def as_terms(self):
        if len(self.keys()) == 0:
            return [sympify_sympy(S.Zero)]
        return [((self[k]) if k == self.ring.zero_monom else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == self.ring.zero_monom else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys())]

    def __add__(self, other):
        if isinstance(other, BaseRingElement):
            return self.ring.add(self, other)
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
        if isinstance(other, BaseRingElement):
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
        return self.as_expr()

    def __eq__(self, other):
        return type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)

    def almosteq(self, other):
        if isinstance(other, BaseRingElement):
            elem1 = self
            elem2 = other
            if elem1.ring == elem2.ring:
                return (self - other).expand(deep=False).almosteq(S.Zero)
            return elem1.almosteq(elem1.ring.one * elem2)
        return (self - self.ring.from_expr(other)).expand(deep=False) == self.ring.zero

    def __matmul__(self, other):
        return (self.ring @ self.ring).ext_multiply(self, other)


class BaseRing(Ring, CompositeDomain):
    def __str__(self):
        return self.__class__.__name__

    def __matmul__(self, other):
        from .tensor_ring import TensorRing

        return TensorRing(self, other)

    def __eq__(self, other):
        return type(self) is type(other)

    def to_sympy(self, elem):
        return elem.as_expr()

    def __init__(self, domain=None):
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain
        self.zero_monom = ()

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
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            if isinstance(other, BaseRingElement):
                coerced = self._coerce_mul(other)
                if not coerced:
                    raise CoercionFailed(f"Could not coerce {other} of type {type(other)} to {type(elem)}")
                return self._mul_elements(elem, coerced)
            return self.mul_expr(elem, other)

    def _mul_elements(self, elem, other):
        raise NotImplementedError

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
        return self.from_dict({self.zero_monom: S.One})

    def is_elem_mul_type(self, elem): ...

    def elem_mul(self, ring_elem, elem): ...

    def _coerce_add(self, other): ...

    def from_dict(self, element, orig_domain=None):
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

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        try:
            if isinstance(element, BaseRingElement) or isinstance(element, DomainElement):
                raise CoercionFailed("Not a domain element")
            return sympify(element)
        except Exception:
            raise CoercionFailed(f"Could not coerce {element} of type {type(element)} to {self.__class__.__name__}")

    def from_expr(self, x):
        return self.mul_expr(self.one, x)

    def mul_expr(self, elem, x): ...
