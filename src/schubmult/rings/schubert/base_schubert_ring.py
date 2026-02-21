from sympy import pretty

from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import Add, CoercionFailed, S, expand, sympify, sympify_sympy, sympy_Mul
from schubmult.symbolic.poly.schub_poly import schubpoly_from_elems
from schubmult.utils._mul_utils import _mul_schub_dicts
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from ..base_ring import BaseRing, BaseRingElement

logger = get_logger(__name__)


class BaseSchubertElement(BaseRingElement):

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

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == self.ring.zero_monom else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys(), key=(lambda kk: (kk.inv, tuple(kk)) if hasattr(kk, "inv") else kk))]

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

    def almosteq(self, other):
        if isinstance(other, BaseSchubertElement):
            elem1 = self
            elem2 = other
            if elem1.ring == elem2.ring:
                return (self - other).expand(deep=False).almosteq(S.Zero)
            return elem1.almosteq(elem1.ring.one * elem2)
        return (self - self.ring.from_expr(other)).expand(deep=False) == self.ring.zero


class BaseSchubertRing(BaseRing):

    def __eq__(self, other):
        return type(self) is type(other) and self.genset == other.genset and self.coeff_genset == other.coeff_genset

    def __init__(self, genset, coeff_genset, domain=None):
        super().__init__(domain=domain)
        self._genset = genset
        self._coeff_genset = coeff_genset
        self.symbols = list(genset)
        self.zero_monom = Permutation([])

    def _mul_elements(self, elem, other):
        return self.from_dict(_mul_schub_dicts(elem, other, elem.ring, other.ring))

    def new(self, x): ...

    def printing_term(self, k): ...

    def coproduct_on_basis(self, k): ...

    def _coerce_mul(self, other): ...

    def is_elem_mul_type(self, elem): ...

    def elem_mul(self, ring_elem, elem): ...

    def _coerce_add(self, other): ...

    @property
    def elem_sym(self): ...

    @property
    def symbol_elem_func(self): ...

    def elem_sym_subs(self, kk): ...

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        try:
            if isinstance(element, BaseRingElement):
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
