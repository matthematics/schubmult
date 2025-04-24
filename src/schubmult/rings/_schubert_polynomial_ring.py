from functools import cache, cached_property

import sympy
from symengine import Add, S, Symbol, SympifyError, expand, sympify
from sympy import Basic
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.printing.str import StrPrinter

import schubmult.rings._quantum_schubert_polynomial_ring as qsr
import schubmult.rings._tensor_schub_ring as tsr
import schubmult.rings._utils as utils
import schubmult.schub_lib.double as yz
import schubmult.schub_lib.schub_lib as schub_lib
import schubmult.schub_lib.single as py
from schubmult.perm_lib import Permutation, inv, uncode
from schubmult.poly_lib.poly_lib import complete_sym_poly, elem_sym_poly, xreplace_genvars
from schubmult.poly_lib.schub_poly import schubpoly_classical_from_elems, schubpoly_from_elems
from schubmult.poly_lib.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

# init_logging(True)
## EMULATE POLYTOOLS

_def_printer = StrPrinter({"order": "none"})
# _def_printer = StrPrinter()

logger = get_logger(__name__)

# numpy arrays
# sympy parsing
# quantum

# COPRODUCT


class NotEnoughGeneratorsError(ValueError):
    pass


def _mul_schub_dicts(dict1, dict2, basis1, basis2, best_effort_positive=True):
    this_dict = {}
    for k, v in dict2.items():
        for kd, vd in dict1.items():
            did_positive = False
            to_mul = v * vd
            if best_effort_positive:
                try:
                    # logger.critical(f"{to_mul=} {kd=} {k=}")
                    this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis1.cached_positive_product(kd, k, basis2).items()})
                    did_positive = True
                except Exception:
                    # # # logger.debug(("Failed to compute")
                    did_positive = False
            if not did_positive:
                this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis1.cached_product(kd, k, basis2).items()})
    # results = add_perm_dict(results, this_dict)

    return this_dict


class BasisSchubertAlgebraElement(Expr):
    _op_priority = 1e200
    _kind = NumberKind
    is_commutative = False
    # precedence = 40
    do_parallel = False

    def __new__(cls, _dict, basis):
        obj = Expr.__new__(cls)
        obj._dict = {k: sympify(v) for k, v in _dict.items() if expand(v) != S.Zero}
        if len(obj._dict.keys()) == 1 and next(iter(obj._dict.values())) == S.One:
            obj.precedence = 1000
        else:
            obj.precedence = 40
        obj._basis = basis
        return obj

    # 217 per night
    # 569 per night
    @property
    def args(self):
        return (sympy.Dict(self._dict), self._basis)

    @property
    def coeff_dict(self):
        return self._dict

    @property
    def genset(self):
        return self.basis.genset

    @property
    def coeff_genset(self):
        return self.basis.coeff_genset

    @property
    def basis(self):
        return self._basis

    def _hashable_content(self):
        return self.args

    # def prune(self):
    #     keys = list(self._dict.keys())
    #     for k in keys:
    #         if expand(self._dict[k]) == S.Zero:
    #             del self._dict[k]
    #     return self

    def mult_poly(self, poly):
        res_dict2 = {}
        # poly = self.genset[i + 1] - self.genset[i]
        for k, v in self.coeff_dict.items():
            if self.basis.coeff_genset.label is None:
                dict2 = self.basis.mult_poly_single({k: v}, poly, self.genset)
            else:
                dict2 = self.basis.mult_poly_double({k: v}, poly, self.genset, self.basis.coeff_genset)
            res_dict2 = add_perm_dict(res_dict2, dict2)
        # # # logger.debug((f"{res_dict2=}")
        return self.basis._from_dict(res_dict2)

    def in_SEM_basis(self):
        result = S.Zero
        for k, v in self.coeff_dict.items():
            result += v * schubpoly_from_elems(k, self.genset, self.basis.coeff_genset, elem_func=self.basis.symbol_elem_func)
        return result

    def in_CEM_basis(self):
        result = S.Zero
        for k, v in self.coeff_dict.items():
            result += v * schubpoly_classical_from_elems(k, self.genset, self.coeff_genset, elem_func=self.basis.symbol_elem_func)
        return result

    def _sympystr(self, printer):
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(self)

    def _pretty(self, printer):
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(self)

    def _latex(self, printer):
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(self)

    def as_terms(self):
        if len(self.coeff_dict.keys()) == 0:
            return [sympy.sympify(S.Zero)]
        return [
            (sympy.sympify(self.coeff_dict[k]) if k == Permutation([]) else sympy.Mul(self.coeff_dict[k], self.basis.single_element_class(k, self.basis)))
            for k in sorted(self.coeff_dict.keys(), key=lambda bob: (inv(bob), *bob))
        ]

    def as_ordered_terms(self, *_, **__):
        return self.as_terms()

    # def _eval_simplify(self, *args, measure, **kwargs):
    #     return self.basis._from_dict({k: sympify(sympy.simplify(v, *args, measure=measure, **kwargs)) for k, v in self.coeff_dict.items()})

    # def __iadd__(self, other):
    #     return self.__add__(other)

    def __add__(self, other):
        # if isinstance(self)
        # # # # logger.debug((f"{type(other)=} {self.genset=}")
        if isinstance(other, BasisSchubertAlgebraElement):
            new_other = self.basis._coerce_add(other)
            if new_other:
                return self.basis._from_dict(add_perm_dict(self.coeff_dict, new_other.coeff_dict))
            # new_self = other.basis._coerce_add(self)
            # if new_self:
            #     return other.basis._from_dict(add_perm_dict(new_self.coeff_dict, other.coeff_dict))
            return other.__radd__(self)
        # # logger.debug(f"{self=} {other=} {self.basis=}")
        return super().__add__(other)

    def __radd__(self, other):
        # # # logger.debug((f"{type(other)=}")
        if isinstance(other, BasisSchubertAlgebraElement):
            new_other = self.basis._coerce_add(other)
            if new_other:
                return self.basis._from_dict(add_perm_dict(new_other.coeff_dict, self.coeff_dict))
        # # logger.debug(f"{self=} {other=} {self.basis=}")
        return super().__radd__(other)

    def __sub__(self, other):
        # # # logger.debug((f"{type(other)=}")
        if isinstance(other, BasisSchubertAlgebraElement):
            new_other = self.basis._coerce_add(other)
            if new_other:
                return self.basis._from_dict(add_perm_dict(self.coeff_dict, {k: -v for k, v in new_other.coeff_dict.items()}))
            return other.__rsub__(self)
        return super().__sub__(other)

    def __rsub__(self, other):
        if isinstance(other, BasisSchubertAlgebraElement):
            new_other = self.basis._coerce_add(other)
            if new_other:
                return self.basis._from_dict(add_perm_dict(other.coeff_dict, {k: -v for k, v in new_other.coeff_dict.items()}))
        return super().__rsub__(other)

    def __neg__(self):
        return self.basis._from_dict({k: -v for k, v in self.coeff_dict.items()})

    def __mul__(self, other):
        if isinstance(other, BasisSchubertAlgebraElement):
            new_other = self.basis._coerce_mul(other)
            if new_other:
                return self.basis._from_dict(_mul_schub_dicts(self.coeff_dict, new_other.coeff_dict, self.basis, new_other.basis))
            return other.__rmul__(self)
        try:
            other = sympify(other)
            if not any(x in self.basis.genset for x in other.free_symbols):
                return self.basis._from_dict({k: other * v for k, v in self.coeff_dict.items()})
            new_other = self.basis(other)
            if new_other:
                return self.__mul__(new_other)
        except SympifyError:
            if isinstance(other, sympy.Expr):
                try:
                    return other.__mul__(sympify(self))
                except Exception:
                    pass
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, BasisSchubertAlgebraElement):
            # # logger.debug(f"{type(self.basis)=} {type(other.basis)=}")
            new_other = self.basis._coerce_mul(other)
            # # logger.debug(f"{new_other=}")
            if new_other:
                return self.basis._from_dict(_mul_schub_dicts(new_other.coeff_dict, self.coeff_dict, new_other.basis, self.basis))
            return NotImplemented
        try:
            other = sympify(other)
            if not any(x in self.basis.genset for x in other.free_symbols):
                return self.basis._from_dict({k: other * v for k, v in self.coeff_dict.items()})
            new_other = self.basis(other)
            if new_other:
                return self.__rmul__(new_other)
        except SympifyError:
            if isinstance(other, sympy.Expr):
                # # logger.debug(f"{self=} {other=}")
                # # logger.debug(f"{self.basis=} {other.basis=}")
                return self.__rmul__(sympify(other))
        return NotImplemented

    # def equals(self, other):
    #     return self.__eq__(other)

    # def test_equality(self, other, disp=False):
    #     elem1 = self
    #     elem2 = other
    #     done = set()
    #     import sys

    #     for k, v in elem1.coeff_dict.items():
    #         done.add(k)
    #         if expand(v - elem2.coeff_dict.get(k, 0)) != 0:
    #             if disp:
    #                 # print(f"{k=} {v=} {elem2.coeff_dict.get(k, 0)=} {expand(v - elem2.coeff_dict.get(k, 0))=}", file=sys.stderr)
    #             return False
    #     for k, v in elem2.coeff_dict.items():
    #         if k in done:
    #             continue
    #         if expand(v - elem1.coeff_dict.get(k, 0)) != 0:
    #             if disp:
    #                 # print(f"{k=} {v=} {expand(v - elem1.coeff_dict.get(k, 0))=}", file=sys.stderr)
    #             return False
    #     return True

    def as_coefficients_dict(self):
        return sympy.Dict({self.basis.single_element_class(k, self.basis): sympy.sympify(v) for k, v in self.coeff_dict.items()})

    def _eval_expand_basic(self, *args, **kwargs):  # noqa: ARG002
        return self.as_polynomial()

    # def _eval_expand_mul(self, *args, **kwargs):
    #     # print(f"wogboppy {self} {args=} {kwargs=}")
    #     return self

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        # print(f"Frable gaboopa {self=} {args=} {kwargs=}")
        if not deep:
            return self.basis._from_dict({k: expand(v) for k, v in self.coeff_dict.items()})
        return sympy.sympify(expand(sympify(self.as_polynomial())))

    def as_polynomial(self):
        return sympy.sympify(Add(*[v * self.basis.cached_schubpoly(k) for k, v in self.coeff_dict.items()]))

    def as_classical(self):
        return self.basis.in_classical_basis(self)

    def as_quantum(self):
        return self.basis.in_quantum_basis(self)


class DoubleSchubertAlgebraElement(BasisSchubertAlgebraElement):
    """Algebra with sympy coefficients
    and a dict basis
    """

    # __slots__ = ("_dict", "_parent")
    # is_polynomial = True

    def __new__(cls, _dict, basis):
        return BasisSchubertAlgebraElement.__new__(cls, _dict, basis)

    def divdiff(self, i):
        return self.basis._from_dict({k.swap(i - 1, i): v for k, v in self.coeff_dict.items() if i - 1 in k.descents()})

    def simpleref(self, i):
        return self + self.divdiff(i).mult_poly(self.genset[i + 1] - self.genset[i])

    def act(self, perm):
        perm = Permutation(perm)
        dset = perm.descents()
        if len(dset) == 0:
            return self
        i = next(iter(dset))
        return self.simpleref(i + 1).act(perm.swap(i, i + 1))

    def max_index(self):
        return max([max([0, *list(k.descents(zero_indexed=False))]) for k in self.coeff_dict.keys()])

    def subs(self, old, new):
        result = 0
        if self.genset.index(old) != -1:
            result = 0
            index = self.genset.index(old)
            mindex = self.max_index()
            if mindex < index:
                return self
            # if already equal to the max index, we don't want to move it over
            perm = Permutation([]).swap(index - 1, mindex)  # index to max index + 1
            # # # logger.debug((f"{mindex=}")
            # # # logger.debug((f"{perm=}")
            transf = self.act(perm)
            # # # # logger.debug((f"{transf=}")
            # # # # logger.debug((f"{self.expand()=}")
            # # # # logger.debug((f"{transf.expand().expand()=}")
            # # # # logger.debug((f"{transf2=}")
            # for (k1, k2), v in transf2.coeff_dict.items():
            #     result += self.basis._from_dict({k1: v}) * (new**k2[0].inv)
            # don't want to go nuts
            # res_dict = {}
            for k, v in transf.coeff_dict.items():
                perm = k
                coeff_gens = self.basis.coeff_genset
                # cached mul_poly
                L = schub_lib.pull_out_var(mindex + 1, perm)
                # # # # logger.debug((f"{perm=} {L=}")
                for index_list, new_perm in L:
                    result += self.basis._from_dict({new_perm: v}).mult_poly(sympy.prod([(new - coeff_gens[index2]) for index2 in index_list]))
            return result

        for k, v in self.coeff_dict.items():
            if self.coeff_genset.label is None:
                add_dict = {k: v.subs(old, new)}
            else:
                coeff_genset = self.coeff_genset
                if coeff_genset.index(old) != -1:
                    genset_list = [coeff_genset[i] for i in range(len(coeff_genset))]
                    # print(f"{genset_list=}")
                    genset_list[coeff_genset.index(old)] = 0
                    # print(f"{genset_list=}")
                    custom_genset = CustomGeneratingSet(genset_list)
                    # print(f"{custom_genset=}")
                    new_add_dict = {k2: sympify(v2).subs(old, new) for k2, v2 in yz.schubmult_double({(): v}, k, custom_genset, coeff_genset).items()}  # remove the variable
                    # print(f"{new_add_dict=}")
                    add_dict = {}
                    for k3, v3 in new_add_dict.items():
                        # convert back to coeff_genset
                        to_add_dict = yz.schubmult_double({(): v3}, k3, coeff_genset, custom_genset)
                        # print(f"{to_add_dict=}")
                        add_dict = add_perm_dict(add_dict, to_add_dict)
                else:
                    add_dict = {k: sympify(v).subs(old, new)}
            for k5, v5 in add_dict.items():
                if any(self.genset.index(s) != -1 for s in sympify(v5).free_symbols):
                    result += self.basis._from_dict({k5: 1}).mult_poly(v5)
                else:
                    result += self.basis._from_dict({k5: v5})
            # check correct, change vars to zeroed coeff var for coeff
        return result

    @property
    def free_symbols(self):
        ret = set()
        for k, v in self.coeff_dict.items():
            ret.update(v.free_symbols)
            perm = k
            if len(perm.descents()) > 0:
                ret.update([self.genset[i] for i in range(1, max(perm.descents()) + 2)])
            if self.coeff_genset:
                genset2 = self.basis.coeff_genset
                perm2 = ~perm
                if len(perm2.descents()) > 0:
                    ret.update([genset2[i] for i in range(1, max(perm2.descents()) + 2)])
        return ret

    def pull_out_gen(self, gen):
        ind = self.genset.index(gen)
        gens2 = MaskedGeneratingSet(self.genset, [ind])
        gens2.set_label(f"({self.genset.label}\\{gen})")
        new_basis = DoubleSchubertAlgebraElement_basis(gens2)
        ret = new_basis(0)
        for (perm, cv), val in self.coeff_dict.items():
            L = schub_lib.pull_out_var(ind, perm)
            for index_list, new_perm in L:
                toadd = S.One
                for index2 in index_list:
                    toadd *= gen - utils.poly_ring(cv)[index2]
                ret += toadd * val * new_basis(new_perm, cv)
        return ret

    def coproduct(self, *indices, alt_coeff_genset=None, on_coeff_gens=False, gname1=None, gname2=None):
        result_dict = {}
        genset = self.genset
        if on_coeff_gens:
            genset = self.coeff_genset
        # if not alt_coeff_genset:
        #     alt_coeff_genset = self.coeff_genset
        if gname1 is None:
            gname1 = f"{genset.label}_A"
        if gname2 is None:
            gname2 = f"{genset.label}_B"
        gens2 = MaskedGeneratingSet(genset, indices)
        # # # logger.debug((f"{indices=}")
        gens1 = gens2.complement()
        # # # logger.debug((f"{gens1.index_mask=}")
        # # # logger.debug((f"{list(gens1)=}")
        # # # logger.debug((f"{gens2.index_mask=}")
        # # # logger.debug((f"{list(gens2)=}")
        gens1.set_label(gname1)
        gens2.set_label(gname2)
        for k, v in self.coeff_dict.items():
            key = k
            if isinstance(self.basis, SchubertAlgebraElement_basis) and not alt_coeff_genset:
                coprod_dict = py.schub_coprod_py(key, indices)
            else:
                if on_coeff_gens:
                    coprod_dict = yz.schub_coprod_double(~key, indices, self.basis.genset, alt_coeff_genset if alt_coeff_genset else self.basis.genset)
                else:
                    coprod_dict = yz.schub_coprod_double(key, indices, self.basis.coeff_genset, alt_coeff_genset if alt_coeff_genset else self.basis.coeff_genset)
            # print(f"{coprod_dict=}")
            if on_coeff_gens:
                result_dict = add_perm_dict(result_dict, {(~k1, ~k2): v * v2 for (k1, k2), v2 in coprod_dict.items()})
            else:
                result_dict = add_perm_dict(result_dict, {k: v * v2 for k, v2 in coprod_dict.items()})
        if on_coeff_gens:
            basis = tsr.TensorAlgebraBasis(
                DoubleSchubertAlgebraElement_basis(self.basis.genset, gens1),
                DoubleSchubertAlgebraElement_basis(alt_coeff_genset if alt_coeff_genset else self.basis.genset, gens2),
            )
        else:
            basis = tsr.TensorAlgebraBasis(
                DoubleSchubertAlgebraElement_basis(gens1, self.basis.coeff_genset),
                DoubleSchubertAlgebraElement_basis(gens2, alt_coeff_genset if alt_coeff_genset else self.basis.coeff_genset),
            )
        return basis._from_dict(result_dict)

    @cached_property
    def max_gens(self):
        return max([max(k.descents()) for k in self.coeff_dict.keys()])


_pretty_schub_char = "𝔖"  # noqa: RUF001


# Atomic Schubert polynomial
class DSchubPoly(DoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, basis):
        return DSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        _coeff_dict = sympy.Dict({Permutation(k): 1})
        # if not isinstance(genset, GeneratingSet_base):
        #     raise TypeError
        obj = DoubleSchubertAlgebraElement.__new__(_class, _coeff_dict, basis)
        obj._key = k
        obj._genset = basis.genset
        obj._coeff_dict = _coeff_dict
        obj._basis = basis
        return obj

    # @property
    # def coeff_dict(self):
    #     return self._coeff_dict

    @property
    def perm(self):
        return self._key

    @property
    def args(self):
        return (sympy.Tuple(*self._key), self._basis)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, basis):
        return DSchubPoly.__xnew__(_class, k, basis)

    def _sympystr(self, printer):
        if self._key == Permutation([]):
            return printer.doprint(1)
        if self.coeff_genset.label is None:
            return printer.doprint(f"S{self.genset.label}({printer.doprint(self._key)})")
        return printer.doprint(f"DS{self.genset.label}({printer.doprint(self._key)}, {self.coeff_genset.label})")

    def _pretty(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(int("".join([str(i) for i in self._key])))
        if self.coeff_genset.label is None:
            return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(self.genset.label)))
        return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(f"{self.genset.label}; {self.coeff_genset.label}")))

    def _latex(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(self._key)
        if self.coeff_genset.label is None:
            return printer._print_Function(sympy.Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(self.genset.label)))
        return printer._print_Function(sympy.Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(f"{self.genset.label}; {self.coeff_genset.label}")))




# can coerce, otherwise unevaluated mul
class DoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset, coeff_genset):
        return DoubleSchubertAlgebraElement_basis.__xnew_cached__(cls, genset, coeff_genset)

    @cache
    def __xnew_cached__(_class, genset, coeff_genset):
        return DoubleSchubertAlgebraElement_basis.__xnew__(_class, genset, coeff_genset)

    def __xnew__(_class, genset, coeff_genset):
        return Basic.__new__(_class, genset, coeff_genset)

    def _coerce_mul(self, other):
        if isinstance(other, BasisSchubertAlgebraElement):
            if isinstance(other, DoubleSchubertAlgebraElement):
                if self.genset == other.basis.genset:
                    return other
            if isinstance(other.basis, qsr.QuantumDoubleSchubertAlgebraElement_basis):
                return self._coerce_mul(other.as_classical())
        return None

    def _coerce_add(self, other):
        if isinstance(other, BasisSchubertAlgebraElement):
            if type(other.basis) is type(self):
                if self.genset == other.basis.genset and self.coeff_genset == other.basis.coeff_genset:
                    return other
        return None

    def express(self, other):
        if other.basis == self:
            return other
        if type(other.basis) is DoubleSchubertAlgebraElement_basis:
            if other.genset == self.genset:
                return self([]) * other
        raise NotImplementedError

    # elem syms as functions
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
                elems += [(sympy.Symbol(f"e_{p}_{k}"), elem_sym_poly(p, k, self.genset[1:], utils.poly_ring(0)))]
        return dict(elems)

    # def in_SEM_basis(self, elem):
    #     return

    @property
    def genset(self):
        return self.args[0]

    @property
    def coeff_genset(self):
        return self.args[1]

    def _from_dict(self, _dict):
        return DoubleSchubertAlgebraElement(_dict, self)

    @property
    def single_element_class(self):
        return DSchubPoly

    def in_quantum_basis(self, elem):
        result = S.Zero
        for k, v in elem.coeff_dict.items():
            result += v * self.quantum_schubpoly(k)
            # print(f"{v=} {result=}")
        return result

    def in_classical_basis(self, elem):
        return elem

    @cache
    def quantum_schubpoly(self, perm):
        # print(f"quantum bucket schubpoly {self=} {self.genset=} {self.coeff_genset=}")
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
        basis = qsr.QuantumDoubleSchubertAlgebraElement_basis(self.genset, self.coeff_genset)

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
            # print(f"{res=}")
            # # # logger.debug((f"{res=}")
            return res

        return elem_func

    # @cache
    # def cached_schubpoly_oink(self, u):
    #     return yz.shubpoly(u)

    def monomial_schub(self, monom):
        monom = [*monom]
        while len(monom) > 0 and monom[-1] == 0:
            monom.pop()
        return self._monomial_schub_cache(tuple(monom))

    @cache
    def _monomial_schub_cache(self, monom):
        srt_perm = Permutation.sorting_perm([-i for i in monom])
        schub_perm = uncode(sorted(monom, reverse=True))
        return self._from_dict({(schub_perm, utils.NoneVar): S.One}).act(srt_perm)

    @cache
    def cached_schubpoly(self, k):
        # return yz.schubpoly(u)
        # print(f"Brablesink {self.genset=} {self.coeff_genset=}")
        return schubpoly_classical_from_elems(k, self.genset, self.coeff_genset, elem_func=elem_sym_poly)

    def __call__(self, x):
        # print(f"frivol {x=} {cv=}")
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        # # # logger.debug((f"{type(x)=}")
        # if isinstance(x, Mul) or isinstance(x, Add):
        #     raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self._from_dict({p_x: 1})
        elif isinstance(x, Permutation):
            if max([0, *list(x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self._from_dict({x: 1})

        elif isinstance(x, DoubleSchubertAlgebraElement):
            # # # logger.debug(("Line record")
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
        #         # srt_perm.reverse()
        #         # srt_perm = Permutation(srt_perm)
        #         # print(sorted(monom,reverse=True))
        #         schub_perm = uncode(sorted(monom, reverse=True))
        #         result += self._from_dict({(schub_perm, utils.NoneVar): coeff}).act(srt_perm)
        #     return result
        else:
            # # # logger.debug((f"{x=}")
            x = sympify(x)
            result = yz.mult_poly_double({Permutation([]): 1}, x, genset, self.coeff_genset)
            elem = self._from_dict(result)
            # # # logger.debug((f"Returning {elem=}")
        return elem




DoubleSchubertPolynomial = DoubleSchubertAlgebraElement


# class SchubAdd(Add):
#     is_Add = True

#     def __new__(cls, *args, evaluate=False, _sympify=True, **_):
#         obj = sympy.Add.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)
#         obj._args = args
#         if evaluate:
#             return obj.doit()
#         return obj

#     @property
#     def args(self):
#         return self._args

#     def doit(self):
#         ret = self.args[0]
#         # # # logger.debug((f"ADD {self.args=}")
#         for arg in self.args[1:]:
#             # # # logger.debug((f"{arg=} {type(arg)=}")
#             # # # logger.debug((f"{ret=} {type(ret)=}")
#             ret += sympy.expand(arg)
#         return ret

#     def _sympystr(self, printer):
#         return printer._print_Add(self)

#     def expand(self, deep=True, *_, **__):
#         return SchubAdd(*[sympy.expand(arg) for arg in self.args]).doit()


# class SchubMul(sympy.Mul):
#     is_Mul = True

#     def __new__(cls, *args, evaluate=False, _sympify=True, **_):
#         # args, a, b = Mul.flatten(list(args))
#         # if len(args) == 0:
#         #     return 1
#         obj = Mul.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)
#         obj._args = args
#         if evaluate:
#             return obj.doit()
#         return obj

#     @property
#     def args(self):
#         return self._args

#     def doit(self):
#         ret = self.args[0]
#         # # # logger.debug((f"MUL {self.args=}")
#         for arg in self.args[1:]:
#             # # # logger.debug((f"{arg=} {type(arg)=}")
#             ret *= sympy.expand(arg)
#         return ret

#     def _sympystr(self, printer):
#         return printer._print_Mul(self)

#     def __neg__(self):
#         return SchubMul(sympy.Integer(-1), self)

#     def _eval_expand_mul(self, *_, **__):
#         # # # logger.debug((f"Pringles {self.args=}")
#         return SchubMul(*[sympy.expand(arg) for arg in self.args]).doit()


# Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
#     "Mul": [get_postprocessor(Mul)],
#     "Add": [get_postprocessor(Add)],
# }


def DSx(x, genset=GeneratingSet("y")):
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    return DoubleSchubertAlgebraElement_basis(GeneratingSet("x"), genset)(x)


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


# def Sx(x):
#     return DSx(x, utils.NoneVar)


class SchubertAlgebraElement_basis(DoubleSchubertAlgebraElement_basis):
    def __new__(cls, genset):
        return DoubleSchubertAlgebraElement_basis.__new__(cls, genset, utils.poly_ring(utils.NoneVar))

    # def __hash__(self):
    #     return hash(*self.args)

    # def _from_single_dict(self, _dict):
    #     return DoubleSchubertAlgebraElement({(k, utils.NoneVar): v for k, v in _dict.items()}, self)
    # maybe don't reinvent the wheel
    # def _coerce_mul(self, other):
    #     """Coerce a basis schubert algebra element so it can be multiplied

    #     Args:
    #         other (_type_): _description_

    #     Returns:
    #         _type_: _description_
    #     """
    #     if type(other.basis) is type(self):
    #         if self.genset == other.basis.genset:
    #             return other
    #     if type(other.basis) is DoubleSchubertAlgebraElement_basis:
    #         if self.genset == other.basis.genset:
    #             return other
    #     if isinstance(other.basis, qsr.QuantumDoubleSchubertAlgebraElement_basis):
    #         return self._coerce_mul(other.as_classical())
    #     return None

    @property
    def coeff_genset(self):
        return utils.poly_ring(utils.NoneVar)

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            return py.schubmult_py({u: S.One}, v)
        return {k: xreplace_genvars(x, utils.poly_ring(0), basis2.coeff_genset if basis2.coeff_genset else utils.poly_ring(0)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def __call__(self, x):
        genset = self.genset
        # # # logger.debug((f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self._from_dict({Permutation(x): 1})
        elif isinstance(x, Permutation):
            elem = self._from_dict({x: 1})
        # elif isinstance(x, spr.SchubertPolynomial):
        #     if x._parent._base_var == self._base_var:
        #         elem_dict = {(x, utils.NoneVar): v for k, v in x.coeff_dict.items()}
        #         elem = QuantumDoubleSchubertAlgebraElement(elem_dict, self)
        #         if cv is not None:
        #             elem = self([1, 2], cv) * elem
        #     else:
        #         return self(x.expand(), cv)
        elif isinstance(x, DoubleSchubertAlgebraElement):
            if x.is_Add or x.is_Mul:
                return x
            if x.genset == genset:
                elem = DoubleSchubertAlgebraElement(x.coeff_dict, self)  # , self)
            else:
                return self(x.expand())
        # elif isinstance(x, spr.DoubleSchubertAlgebraElement):
        #     if x.genset == self.genset:
        #         return x.as_quantum()
        else:
            x = sympify(x)
            result = py.mult_poly_py({Permutation([]): 1}, x, genset)
            elem = self._from_dict(result)
        return elem


Sx = SchubertAlgebraElement_basis(GeneratingSet("x"))
