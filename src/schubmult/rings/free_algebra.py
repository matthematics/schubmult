from functools import cache

from schubmult.perm_lib import uncode
from schubmult.symbolic import (
    EXRAW,
    Add,
    CoercionFailed,
    CompositeDomain,
    DefaultPrinting,
    DomainElement,
    Integer,
    Mul,
    Ring,
    S,
    Symbol,
    expand,
    is_of_func_type,
    sstr,
    sympify,
    sympify_sympy,
    sympy_Add,
    sympy_Mul,
)
from schubmult.symmetric_polynomials import FactorialElemSym
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from .abstract_schub_poly import GenericPrintingTerm
from .base_schubert_ring import BaseSchubertElement
from .schubert_ring import DSx, Sx
from .separated_descents import SeparatedDescentsRing

ASx = SeparatedDescentsRing(Sx([]).ring)
ADSx = SeparatedDescentsRing(DSx([]).ring)

logger = get_logger(__name__)


# keys are tuples of nonnegative integers
class FreeAlgebraElement(DomainElement, DefaultPrinting, dict):
    precedence = 40

    __sympy__ = True

    def parent(self):
        return self.ring

    def eval(self, *args):
        pass

    def _sympystr(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(sympy_Add(*self.as_ordered_terms()), order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

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
        return [self[k] if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k)) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys())]

    def __add__(self, other):
        if isinstance(other, FreeAlgebraElement):
            if self.ring == other.ring:
                return self.ring.add(self, other)
            return other.__radd__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
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
            other = self.ring.from_dict({(): other})
            return self.ring.add(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__add__(self)
        except CoercionFailed:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, FreeAlgebraElement):
            if self.ring == other.ring:
                return self.ring.sub(self, other)
            return other.__rsub__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
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
            other = self.ring.from_dict({(): other})
            return self.ring.sub(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__sub__(self)
        except CoercionFailed:
            return NotImplemented

    def __neg__(self):
        return self.ring.neg(self)

    def __mul__(self, other):
        try:
            return self.ring.mul(self, other)
        except CoercionFailed:
            return other.__rmul__(self)

    def __rmul__(self, other):
        try:
            return self.ring.rmul(self, other)
        except CoercionFailed:
            return NotImplemented

    def as_coefficients_dict(self):
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        return self.ring.from_dict({k: expand(v, **kwargs) for k, v in self.items()})

    @property
    def free_symbols(self):
        return set()

    @staticmethod
    @cache
    def tup_expand(tup):
        res = ASx([])
        if len(tup) == 0:
            return res
        tup2 = tup[:-1]
        return FreeAlgebraElement.tup_expand(tup2) * ASx(uncode([tup[-1]]), 1)

    @staticmethod
    @cache
    def tup_double_expand(tup):
        res = ADSx([])
        if len(tup) == 0:
            return res
        tup2 = tup[:-1]
        return FreeAlgebraElement.tup_double_expand(tup2) * ADSx(uncode([tup[-1]]), 1)

    def schub_expand(self):
        res = ASx([]).ring.zero
        for tup, val in self.items():
            res += val * self.__class__.tup_expand(tup)
        return res

    def schub_double_expand(self):
        res = ADSx([]).ring.zero
        for tup, val in self.items():
            res += val * self.__class__.tup_double_expand(tup)
        return res

    def coproduct(self):
        T = self.ring @ self.ring
        res = T.zero

        for key, val in self.items():
            res += val * self.ring.coproduct_on_basis(key)
        return res

    def as_expr(self):
        return Add(*self.as_terms())

    def __eq__(self, other):
        return type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)

    def __str__(self):
        return sstr(self)

    def _nsymtup(self, tup, R):
        if len(tup) == 0:
            return R.one
        if tup[-1] == 0:
            return self._nsymtup(tup[:-1], R)
        return self._nsymtup(tup[:-1], R) * R.from_dict({(tup[-1],): S.One}, R)

    def remove_zeros(self):
        new_elem = self.ring.zero
        for k, v in self.items():
            new_elem += self.ring.from_dict({tuple([a for a in k if a != 0]): v})
        return new_elem

    def nsymexpand(self):
        R = NSym()
        ret = R.zero
        for k, v in self.items():
            ret += v * self._nsymtup(k, R)
        return ret

    def to_schub(self):
        res = Sx([]).ring.zero
        for k, v in self.items():
            res += v * self.ring.tup_to_schub(k)
        return res


class FreeAlgebra(Ring, CompositeDomain):
    def __str__(self):
        return self.__class__.__name__

    def tensor_schub_expand(self, tensor):
        T = ASx([]).ring @ ASx([]).ring
        ret = T.zero
        for (key1, key2), val in tensor.items():
            ret += val * T.ext_multiply(self(key1).schub_expand(), self(key2).schub_expand())
        return ret

    def __hash__(self):
        return hash((self.domain, "whatabong"))

    def __eq__(self, other):
        return type(self) is type(other) and self.domain == other.domain

    def __init__(self, domain=None):
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain
        self.zero_monom = ()
        self.dtype = type("FreeAlgebraElement", (FreeAlgebraElement,), {"ring": self})

    @staticmethod
    def right_pad(tup, n):
        if len(tup) < n:
            return (*tup, *([0] * (n - len(tup))))
        return tup

    def __matmul__(self, other):
        from .tensor_ring import TensorRing

        return TensorRing(self, other)

    def coproduct_on_basis(self, key):
        T = self @ self
        if len(key) == 0:
            return T.one
        return self.__class__._single_coprod(key[0], T) * self.coproduct_on_basis(key[1:])

    @cache
    def _single_coprod(p, T):
        res = T.zero
        for i in range(p + 1):
            res += T.from_dict({((i,), (p - i,)): S.One})
        return res

    # @cache
    #     def coproduct(self, key):
    #         T = self @ self
    #         val = self(*key)

    #         cprd_val = T.zero

    #         while val != val.ring.zero:
    #             mx = [k[0].code for k in val.keys() if val[k] != S.Zero ]
    #             mx.sort(reverse=True)
    #             cd = mx[0]

    #             mx_key = next(iter([k for k in val.keys() if k[0].code == cd]))
    #             if len(cd) == 0:
    #                 return cprd_val + T.from_dict({((Permutation([]),mx_key[1]),(Permutation([]),mx_key[1])): val[mx_key] * S.One})
    #             cd = [*cd]
    #             fv = cd.pop(0)
    #             while len(cd) > 1 and cd[-1] == 0:
    #                 cd.pop()
    #             cf = val[mx_key]
    #             cprd_val += (T.from_dict({((Permutation([]),0),(Permutation([]),0)): cf*S.One}))*_single_coprod(fv, 1, T) * self.coproduct((uncode(cd), mx_key[1] - 1))
    #             val -= cf * self(uncode([fv]),1)*self(uncode(cd), mx_key[1] - 1)
    #         return cprd_val

    def schub_elem(self, perm, numvars):
        res = self.zero
        expr = Sx(perm * ~uncode(list(range(perm.inv + numvars, perm.inv, -1)))).in_SEM_basis().expand()
        args = expr.args
        if not isinstance(expr, Add):
            args = [expr]
        for arg in args:
            tup = list(range(perm.inv + numvars, perm.inv, -1))
            coeff = S.One
            if is_of_func_type(sympify(arg), FactorialElemSym):
                arg = sympify_sympy(arg)
                tup[perm.inv + numvars - arg.numvars] = arg.numvars - arg.degree
            elif isinstance(arg, Mul):
                for arg0 in arg.args:
                    if is_of_func_type(sympify(arg0), FactorialElemSym):
                        arg0 = sympify_sympy(arg0)
                        tup[perm.inv + numvars - arg0.numvars] = arg0.numvars - arg0.degree
                    else:
                        coeff = Integer(arg0)
            else:
                coeff = Integer(arg)
            res += coeff * self(tup)
        return res

    def tup_to_schub(self, tup):
        from schubmult.abc import e, x
        pinv = sum(tup)
        res = Sx([])
        for i in range(len(tup)):
            #numvars = len(tup) + 1 - i
            # i = len - numvars
            # tup[i] = numvars - degree
            # numvars = len - i
            # degree = numvars - tup[i]

            numvars = len(tup) - i + pinv
            degree = numvars - tup[i]

            if degree > 0:
                res *= e(degree, numvars, x[1:])
        return res

    def schub_elem_double(self, perm, n, N):
        ADSx = SeparatedDescentsRing(DSx([]).ring)
        val = ADSx(perm, n)
        res = self.zero

        while any(sympify(v) != S.Zero for v in val.values()):
            mx = [k[0].code for k in val.keys() if val[k] != S.Zero]
            mx.sort(key=lambda bob: (sum(bob), bob), reverse=True)
            cd = mx[0]

            mx_key = next(iter([k for k in val.keys() if k[0].code == cd]))
            if len(cd) == 0:
                return res + val[mx_key] * self((*([0] * n),))
            cd = [*cd]
            fv = cd.pop(0)
            while len(cd) > 1 and cd[-1] == 0:
                cd.pop()
            cf = val[mx_key]
            res += cf * self((fv,)) * self.schub_elem_double(uncode(cd), mx_key[1] - 1, N)
            val -= cf * val.ring(uncode([fv]), 1) * val.ring(uncode(cd), mx_key[1] - 1)
            keys = val.keys()
            for key in list(keys):
                if len(key[0]) > N or expand(val[key]) == S.Zero:
                    del val[key]
        return res

    def add(self, elem, other):
        # print(f"{elem=}")
        # print(f"{other=}")
        return self.from_dict(add_perm_dict(elem, other))

    def sub(self, elem, other):
        return self.from_dict(add_perm_dict(elem, {k: -v for k, v in other.items()}))

    def neg(self, elem):
        return self.from_dict({k: -v for k, v in elem.items()})

    def mul_perm(self, elem, perm):
        dct = {}
        for k, v in elem.items():
            newperm = k * perm
            if newperm.inv == k.inv + perm.inv:
                dct[newperm] = v
        return self.from_dict(dct)

    def rmul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        if isinstance(other, FreeAlgebraElement):
            raise NotImplementedError
        return self.from_dict({k: v * other for k, v in elem.items()})

    def mul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        try:
            other = self.domain_new(other)
            # print(f"{other=} {type(other)=}")
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        if isinstance(other, FreeAlgebraElement):
            ret = self.zero
            for k0, v0 in elem.items():
                for k, v in other.items():
                    new_key = (*k0, *k)
                    ret += self.from_dict({new_key: v * v0})
            return ret
        raise CoercionFailed

    def to_domain(self):
        return self

    def new(self, x):
        if isinstance(x, FreeAlgebraElement):
            return x
        if isinstance(x, list) or isinstance(x, tuple):
            return self.from_dict({tuple(x): S.One})
        return self.mul_scalar(self.one, x)

    def printing_term(self, k):
        if all(a < 10 for a in k):
            return Symbol("[" + "".join([str(a) for a in k]) + "]")
        return Symbol("[" + " ".join([str(a) for a in k]) + "]")

    def _coerce_mul(self, other): ...

    @property
    def one(self):
        return self.from_dict({(): S.One})

    def from_dict(self, element, orig_domain=None):  # noqa: ARG002
        # print(f"{element=}")
        poly = self.zero

        for monom, coeff in element.items():
            if coeff != self.domain.zero:
                poly[monom] = coeff
        return poly

    @property
    def zero(self):
        return self.dtype()

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        # print(f"{element=} {type(element)=} bagels {type(sympify(element))=} {sympify(element).has(*self.symbols)=}")
        if isinstance(element, FreeAlgebraElement) or isinstance(element, BaseSchubertElement):
            raise CoercionFailed("Not a domain element")
        # if not any(x in self.genset for x in sympify_sympy(element).free_symbols):
        return sympify(element)
        # raise CoercionFailed(f"{element} contains an element of the set of generators")

    # @property
    # def genset(self):
    #     return self._genset

    # def from_expr(self, x):
    #     return self.mul_scalar(self.one, x)


FA = FreeAlgebra()


class NSym(FreeAlgebra):
    def __hash__(self):
        return hash((self.domain, "whatabong2"))

    def __init__(self, domain=None):
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain
        self.zero_monom = ()
        self.dtype = type("NSymElement", (NSymElement,), {"ring": self})

    def printing_term(self, k):
        return GenericPrintingTerm(k, "R")

    def rmul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        if isinstance(other, NSymElement):
            raise NotImplementedError
        return self.from_dict({k: v * other for k, v in elem.items()})

    def mul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        try:
            other = self.domain_new(other)
            # print(f"{other=} {type(other)=}")
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        if isinstance(other, NSymElement):
            ret = self.zero
            for k0, v0 in elem.items():
                for k, v in other.items():
                    if len(k0) > 0 and len(k) > 0:
                        new_key0 = (*k0[:-1], k0[-1] + k[0], *k[1:])
                        ret += self.from_dict({new_key0: v * v0})
                    new_key1 = (*k0, *k)
                    ret += self.from_dict({new_key1: v * v0})
            return ret
        raise CoercionFailed


class NSymElement(FreeAlgebraElement):
    precedence = 40

    __sympy__ = True

    def parent(self):
        return self.ring

    def eval(self, *args):
        pass

    def _sympystr(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(sympy_Add(*self.as_ordered_terms()), order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

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
        return [self[k] if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k)) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys())]

    def __add__(self, other):
        if isinstance(other, NSymElement):
            if self.ring == other.ring:
                return self.ring.add(self, other)
            return other.__radd__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
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
            other = self.ring.from_dict({(): other})
            return self.ring.add(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__add__(self)
        except CoercionFailed:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, NSymElement):
            if self.ring == other.ring:
                return self.ring.sub(self, other)
            return other.__rsub__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
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
            other = self.ring.from_dict({(): other})
            return self.ring.sub(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__sub__(self)
        except CoercionFailed:
            return NotImplemented

    def __neg__(self):
        return self.ring.neg(self)

    def __mul__(self, other):
        try:
            return self.ring.mul(self, other)
        except CoercionFailed:
            return other.__rmul__(self)

    def __rmul__(self, other):
        try:
            return self.ring.rmul(self, other)
        except CoercionFailed:
            return NotImplemented

    def as_coefficients_dict(self):
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    @property
    def free_symbols(self):
        return set()

    def as_expr(self):
        return Add(*self.as_terms())

    def __eq__(self, other):
        return type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)

    def __str__(self):
        return sstr(self)
