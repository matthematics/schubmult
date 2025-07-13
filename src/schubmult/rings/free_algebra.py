from functools import cache

from schubmult.perm_lib import Permutation, trimcode, uncode
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

splugSx = SeparatedDescentsRing(Sx([]).ring)
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

    def __hash__(self):
        return hash(set(self.items()))

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
    def tup_double_expand(tup):
        res = ADSx([])
        if len(tup) == 0:
            return res
        tup2 = tup[:-1]
        return FreeAlgebraElement.tup_double_expand(tup2) * ADSx(uncode([tup[-1]]), 1)

    @staticmethod
    @cache
    def tup_expand(tup):
        res = splugSx([])
        if len(tup) == 0:
            return res
        if len(tup) == 1:
            return splugSx(uncode(tup), 1)
        mid = len(tup) // 2
        return FreeAlgebraElement.tup_expand(tup[:mid]) * FreeAlgebraElement.tup_expand(tup[mid:])

    def change_basis(self, other_basis):
        new_ring = FreeAlgebra(basis=other_basis)
        ret = new_ring.zero
        tfunc = self.ring._basis.transition(other_basis)
        for k, v in self.items():
            ret += v * new_ring.from_dict(tfunc(k))
        return ret

    def schub_expand(self):
        res = splugSx([]).ring.zero
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
        if isinstance(other, FreeAlgebraElement):
            if other.ring == self.ring:
                diff = self - other
                return all(v == S.Zero for v in diff.values())
        return False
        # return type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)

    def __str__(self):
        return sstr(self)

    def _nsymtup(self, tup, R):
        if len(tup) == 0:
            return R.one
        if 0 in tup:
            # return self._nsymtup(tup[:-1], R)
            return R.zero
        return self._nsymtup(tup[:-1], R) * R.from_dict({(tup[-1],): S.One}, R)

    def remove_zeros(self, inserter=S.One):
        new_elem = self.ring.zero
        for k, v in self.items():
            # if 0 in k:
            #     continue
            # add x?
            new_elem += self.ring.from_dict({tuple([a for a in k if a != 0]): v}) * (inserter ** len([a for a in k if a == 0]))
        return new_elem

    def nsymexpand(self):
        R = NSym()
        ret = R.zero
        for k, v in self.items():
            ret += v * self._nsymtup(k, R)
        return ret

    def split(self, p):
        T = self.ring @ self.ring
        ret = T.zero
        for tup, val in self.items():
            if len(tup) < p:
                ret += val * T((tup, tuple()))
            else:
                ret += val * T((tup[:p], tup[p:]))
        return ret

    def to_schub(self, sym=False):
        res = Sx([]).ring.zero
        for k, v in self.items():
            res += v * self.ring.tup_to_schub(k, sym=sym)
        return res


class FreeAlgebraBasis:
    def is_key(self, x): ...

    def as_key(self, x): ...

    def product(self, key1, key2, coeff=S.One): ...

    def coproduct(self, key, coeff=S.One): ...

    def transition(self, other_basis): ...

    @property
    def zero_monom(self): ...

    def printing_term(self, key): ...

    @classmethod
    def compose_transition(tkeyfunc, output):
        ret = {}
        for key, v in output.items():
            ret = add_perm_dict(ret, {k: v*v0 for k, v0 in tkeyfunc(key).items()})
        return ret


class WordBasis(FreeAlgebraBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def product(self, key1, key2, coeff=S.One):
        return {(*key1, *key2): coeff}

    @property
    def zero_monom(self):
        return ()

    @cache
    def coproduct(self, key, coeff=S.One):
        if len(key) == 0:
            return {((), ()): coeff}
        if len(key) == 1:
            key = key[0]
            dct = {}
            for i in range(key + 1):
                dct[((i,), (key - i,))] = coeff
            return dct
        mid = len(key) // 2
        cp1 = self.coproduct(key[:mid], coeff)
        cp2 = self.coproduct(key[mid:])
        ret = {}
        for k0, v0 in cp1.items():
            for k1, v1 in cp2.items():
                ret = add_perm_dict(ret, {((*k0[0], *k1[0]), (*k0[1], *k1[1])): v0 * v1})
        return ret

    def printing_term(self, k):
        if all(a < 10 for a in k):
            return Symbol("[" + "".join([str(a) for a in k]) + "]")
        return Symbol("[" + " ".join([str(a) for a in k]) + "]")

    @staticmethod
    @cache
    def tup_expand(tup):
        res = splugSx([])
        if len(tup) == 0:
            return res
        if len(tup) == 1:
            return splugSx(uncode(tup), 1)
        mid = len(tup) // 2
        return WordBasis.tup_expand(tup[:mid]) * WordBasis.tup_expand(tup[mid:])

    def transition_schubert(self, key):
        return dict(WordBasis.tup_expand(key))

    def transition(self, other_basis):
        if isinstance(other_basis, SchubertBasis):
            return self.transition_schubert
        if isinstance(other_basis, WordBasis):
            return lambda x: x
        if isinstance(other_basis, SchubertSchurBasis):
            return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis().transition(SchubertSchurBasis()), self.transition_schubert(x))
        return None


class SchubertBasis(FreeAlgebraBasis):
    def is_key(self, x):
        return (len(x) == 1 and isinstance(x[0], Permutation | list | tuple)) or (len(x) == 2 and isinstance(x[0], Permutation | list | tuple) and isinstance(x[1], int))

    def as_key(self, x):
        # print(f"{x=} {type(x)=}")
        if len(x) == 1:
            perm = Permutation(x[0])
            return (perm, 0) if len(perm.descents()) == 0 else (perm, max(perm.descents()) + 1)
        return (Permutation(x[0]), x[1])

    def product(self, key1, key2, coeff=S.One):
        return dict(coeff * splugSx(*self.as_key(key1)) * splugSx(*self.as_key(key2)))

    @property
    def zero_monom(self):
        return (Permutation([]), 0)

    @cache
    def coproduct(self, key):
        from ._mul_utils import _tensor_product_of_dicts_first

        dct = self.transition_word(*key)
        res = {}
        wbasis = WordBasis()
        for key_word, v in dct.items():
            dct2 = wbasis.coproduct(key_word, v)
            # print(f"{dct2=}")
            for (k1, k2), v2 in dct2.items():
                dct0 = wbasis.transition_schubert(k1)
                dct1 = wbasis.transition_schubert(k2)
                res = add_perm_dict(res, {k: v0 * v2 for k, v0 in _tensor_product_of_dicts_first(dct0, dct1).items()})
        return res

    def transition_schubert_schur(self, *x):
        perm, numvars = x
        # schur(n)schubert(n)schur(m,init)schubert(n,end)
        # expansion vmu^{-1} in m schur (complement) n schubert (times w0)
        extra = len(perm) - numvars

        if extra == 0:
            return {(tuple([0]*numvars),perm,numvars): 1}
        dom = uncode([numvars] * extra + list(range(numvars - 1, 0, -1)))
        tosplit = perm * dom
        dct = Sx(tosplit).coproduct(*list(range(1,extra+1)))
        w0 = uncode(list(range(numvars-1,0,-1)))
        w0s = uncode([numvars]*extra)
        dct2 = {}
        for (lambd, perm0), v in dct.items():
            perm1 = perm0 * w0
            lambd2 = tuple(trimcode(lambd * (~w0s)))
            dct2[(lambd2, perm1)] = v
        return {(tuple(FreeAlgebra.right_pad(k[0],numvars)),k[1]): v for k, v in dct2.items()}



    def transition(self, other_basis):
        if isinstance(other_basis, SchubertBasis):
            return lambda x: x
        if isinstance(other_basis, SchubertSchurBasis):
            return lambda x: self.transition_schubert_schur(*x)
        if isinstance(other_basis, WordBasis):
            return lambda x: self.transition_word(*x)
        return None

    def transition_word(self, perm, numvars):
        res = {}
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
            tup = tuple(tup)
            res[tup] = res.get(tup, S.Zero) + coeff
        return res

    def printing_term(self, k):
        return splugSx([]).ring.printing_term(k)

class SchubertSchurBasis(FreeAlgebraBasis):
    def is_key(self, x):
        return len(x) == 2 and isinstance(x[0], list | tuple) and isinstance(x[1], Permutation | list | tuple)

    def as_key(self, x):
        return (tuple(x[0]), Permutation(x[1]))

    def product(self, key1, key2, coeff=S.One):
        #return dict(coeff * splugSx(*self.as_key(key1)) * splugSx(*self.as_key(key2)))
        # from ._mul_utils import _tensor_product_of_dicts_first
        # part1 = key1[0]
        # part2 = key2[0]
        # sym_part = coeff * splugSx(uncode(part1), len(part1)) * splugSx(uncode(part2), len(part2))
        # perm1 = key1[1]
        # perm2 = key2[1]
        # bym_part = splugSx(perm1, len(part1)) * splugSx(perm2, len(part2))
        # kys = list(sym_part.keys())
        # for k in kys:
        #     if len(k[0].descents()) > 1:
        #         del sym_part[k]
        # sym_part = {tuple(FreeAlgebra.right_pad(trimcode(k[0]),k[1])): v for k, v in sym_part.items()}
        
        # bym_part = {k[0]: v for k, v in bym_part.items()}
        # # for k0, v0 in sym_part.items():
        # #     for k1, v1 in bym_part.items():
        # #         ret = add_perm_dict(ret, {(k0, k1)})
        # return _tensor_product_of_dicts_first(sym_part, bym_part)
        pass

    @property
    def zero_monom(self):
        return ((), Permutation([]))

    @cache
    def coproduct(self, key): ...
        # from ._mul_utils import _tensor_product_of_dicts_first

        # dct = self.transition_word(*key)
        # res = {}
        # wbasis = WordBasis()
        # for key_word, v in dct.items():
        #     dct2 = wbasis.coproduct(key_word, v)
        #     # print(f"{dct2=}")
        #     for (k1, k2), v2 in dct2.items():
        #         dct0 = wbasis.transition_schubert(k1)
        #         dct1 = wbasis.transition_schubert(k2)
        #         res = add_perm_dict(res, {k: v0 * v2 for k, v0 in _tensor_product_of_dicts_first(dct0, dct1).items()})
        # return res

    #def schubert_matrix(self, lambd, perm0, perm):
        # div schubert: perm * (~mu)
        # div lambd: lambd * (~mu0)
        # remaining: div(lambd * (~mu0), perm * (~mu)) on +n perm0 w9
        # multiply the +n perm0 w0

    def transition_schubert(self, lambd, perm):
        #pass
        from .schubert_ring import SingleSchubertRing
        from .variables import MaskedGeneratingSet
        if lambd[-1] == 0:
            return {(perm, len(lambd)): 1}
        numvars = len(lambd)
        extra = lambd[-1] + len(lambd) - 1
        dom = uncode([numvars] * extra)
        grass_perm = uncode(lambd) * dom
        w0 = uncode(list(range(numvars - 1, 0, -1)))
        lower_perm = perm * w0
        dom_perm = uncode(([numvars] * extra) + list(range(numvars - 1, 0, -1)))
        shifted_ring = SingleSchubertRing(MaskedGeneratingSet(Sx([]).ring.genset,list(range(1,extra+1))))
        start_schub = Sx(grass_perm)
        start_schub *= shifted_ring(lower_perm).in_SEM_basis()
        return {(k*~dom_perm, numvars): v for k, v in start_schub.items()}



    @classmethod
    def transition_word(cls, lambd, perm):
        return FreeAlgebraBasis.compose_transition(SchubertBasis.transition(WordBasis()),self.transition_schubert(lambd, perm))

    @classmethod
    def transition(cls, other_basis):
        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(*x)
        if other_basis == WordBasis:
            return lambda x: cls.transition_word(*x)
        if isinstance(other_basis, SchubertSchurBasis):
            return lambda x: x
        return None

    # def transition_word(self, perm, numvars):
    #     res = {}
    #     expr = Sx(perm * ~uncode(list(range(perm.inv + numvars, perm.inv, -1)))).in_SEM_basis().expand()
    #     args = expr.args
    #     if not isinstance(expr, Add):
    #         args = [expr]
    #     for arg in args:
    #         tup = list(range(perm.inv + numvars, perm.inv, -1))
    #         coeff = S.One
    #         if is_of_func_type(sympify(arg), FactorialElemSym):
    #             arg = sympify_sympy(arg)
    #             tup[perm.inv + numvars - arg.numvars] = arg.numvars - arg.degree
    #         elif isinstance(arg, Mul):
    #             for arg0 in arg.args:
    #                 if is_of_func_type(sympify(arg0), FactorialElemSym):
    #                     arg0 = sympify_sympy(arg0)
    #                     tup[perm.inv + numvars - arg0.numvars] = arg0.numvars - arg0.degree
    #                 else:
    #                     coeff = Integer(arg0)
    #         else:
    #             coeff = Integer(arg)
    #         tup = tuple(tup)
    #         res[tup] = res.get(tup, S.Zero) + coeff
    #     return res

    def printing_term(self, k):
        return Symbol(f"SS{sstr(k)}")

class FreeAlgebra(Ring, CompositeDomain):
    def __str__(self):
        return self.__class__.__name__

    def tensor_schub_expand(self, tensor):
        T = splugSx([]).ring @ splugSx([]).ring
        ret = T.zero
        for (key1, key2), val in tensor.items():
            ret += val * T.ext_multiply(self(key1).schub_expand(), self(key2).schub_expand())
        return ret

    def tensor_nsym_expand(self, tensor, inserter=S.Zero):
        T = NSym() @ NSym()
        ret = T.zero
        for (key1, key2), val in tensor.items():
            ret += val * T.ext_multiply(self(key1).remove_zeros(inserter=inserter).nsymexpand(), self(key2).remove_zeros(inserter=inserter).nsymexpand())
        return ret

    def __hash__(self):
        return hash((self.domain, "whatabong"))

    def __eq__(self, other):
        return type(self) is type(other) and self.domain == other.domain

    def __init__(self, basis=WordBasis(), domain=None):
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain
        self._basis = basis
        self.zero_monom = self._basis.zero_monom
        self.dtype = type("FreeAlgebraElement", (FreeAlgebraElement,), {"ring": self})

    @staticmethod
    def right_pad(tup, n):
        if len(tup) < n:
            return (*tup, *([0] * (n - len(tup))))
        return tup

    def __matmul__(self, other):
        from .tensor_ring import TensorRing

        return TensorRing(self, other)

    @cache
    def coproduct_on_basis(self, key):
        T = self @ self
        return T.from_dict(self._basis.coproduct(key))

    def add(self, elem, other):
        return self.from_dict(add_perm_dict(elem, other))

    def sub(self, elem, other):
        return self.from_dict(add_perm_dict(elem, {k: -v for k, v in other.items()}))

    def neg(self, elem):
        return self.from_dict({k: -v for k, v in elem.items()})

    def rmul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        if isinstance(other, FreeAlgebraElement):
            raise NotImplementedError
        return self.from_dict({k: v * other for k, v in elem.items()})

    def mul(self, elem, other):
        try:
            other = self.domain_new(other)
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        if isinstance(other, FreeAlgebraElement):
            ret = self.zero
            for k0, v0 in elem.items():
                for k, v in other.items():
                    ret += self.from_dict(self._basis.product(k0, k, v * v0))
            return ret
        raise CoercionFailed

    def to_domain(self):
        return self

    def new(self, *x):
        if len(x) == 0 and isinstance(x, FreeAlgebraElement):
            return x
        if self._basis.is_key(x):
            return self.from_dict({self._basis.as_key(x): S.One})
        return self.mul_scalar(self.one, x)

    def printing_term(self, k):
        # print(f"pingdunkit {k=}")
        return self._basis.printing_term(k)

    def _coerce_mul(self, other): ...

    @property
    def one(self):
        return self.from_dict({(): S.One})

    def from_dict(self, element):
        poly = self.zero
        for monom, coeff in element.items():
            if coeff != self.domain.zero:
                poly[monom] = coeff
        return poly

    @property
    def zero(self):
        return self.dtype()

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        if isinstance(element, FreeAlgebraElement) or isinstance(element, BaseSchubertElement):
            raise CoercionFailed("Not a domain element")
        return sympify(element)


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
        return GenericPrintingTerm(k, "N")

    def rmul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        if isinstance(other, NSymElement):
            raise NotImplementedError
        return self.from_dict({k: v * other for k, v in elem.items()})

    def sepify(self, elem):
        return splugSx([]).ring.from_dict({(uncode([a - 1 for a in k]), len(k)): v for k, v in elem.items()})

    def from_sep(self, elem):
        dct = {}
        for (k, n), v in elem.items():
            cd = k.code
            if len(cd) < n:
                cd += [0] * (n - len(cd))
            elif len(cd) > n:
                cd = cd[:n]
            cd = tuple([c + 1 for c in cd])
            dct[cd] = v
        return self.from_dict(dct)

    def mul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        try:
            other = self.domain_new(other)
            # print(f"{other=} {type(other)=}")
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        if isinstance(other, NSymElement):
            return self.from_sep(self.sepify(elem) * self.sepify(other))
            # ret = self.zero
            # for k0, v0 in elem.items():
            #     for k, v in other.items():
            #         if len(k0) > 0 and len(k) > 0:
            #             new_key0 = (*k0[:-1], k0[-1] + k[0], *k[1:])
            #             ret += self.from_dict({new_key0: v * v0})
            #         new_key1 = (*k0, *k)
            #         ret += self.from_dict({new_key1: v * v0})
            # return ret
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
