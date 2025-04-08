from functools import cache, cached_property

import sympy
from symengine import expand, sympify
from sympy import Add, Basic, Mul
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.printing.str import StrPrinter

import schubmult.rings._tensor_schub_ring as tsr
import schubmult.rings._utils as utils
import schubmult.schub_lib.double as yz

# from schubmult.poly_lib.schub_poly import pull_out_var
import schubmult.schub_lib.schub_lib as schub_lib
import schubmult.schub_lib.single as py
from schubmult.perm_lib import (
    Permutation,
    add_perm_dict,
    inv,
)
from schubmult.poly_lib.poly_lib import xreplace_genvars
from schubmult.poly_lib.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet
from schubmult.utils.logging import get_logger

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


def _varstr(v):
    if v == utils.NoneVar:
        return "NoneVar"
    if v == utils.ZeroVar:
        return "0"
    return f"'{v}'"


# def self._from_dict(coeff_dict):
#     return DoubleSchubertAlgebraElement(coeff_dict)


@cache
def cached_schubpoly(u):
    return yz.schubpoly(u)


@cache
def cached_product(u, v, va, vb):
    return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}


@cache
def cached_positive_product(u, v, va, vb):
    return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_generic_partial_posify(u, v).items()}


def _mul_schub_dicts(dict1, dict2, best_effort_positive=True):
    by_var = {}

    none_dict = {}
    for k, v in dict1.items():
        if k[1] == utils.NoneVar or k[1] == utils.ZeroVar:
            none_dict[k[0]] = v
        else:
            if k[1] not in by_var:
                by_var[k[1]] = {}
            by_var[k[1]][k[0]] = v

    results = {}
    # import sys

    for _vstr, _dict in by_var.items():
        this_dict = {}
        for k, v in dict2.items():
            for kd, vd in _dict.items():
                did_positive = False
                if best_effort_positive:
                    try:
                        this_dict = add_perm_dict(this_dict, {k1: v1 * v * vd for k1, v1 in cached_positive_product(kd, k[0], _vstr, k[1]).items()})
                        did_positive = True
                    except Exception:
                        logger.debug("Failed to compute")
                        did_positive = False
                if not did_positive:
                    this_dict = add_perm_dict(this_dict, {k1: v1 * v * vd for k1, v1 in cached_product(kd, k[0], _vstr, k[1]).items()})
        results = add_perm_dict(results, this_dict)

    by_var2 = {}
    none_dict2 = {}
    for k, v in dict2.items():
        if k[1] == utils.NoneVar or k[1] == utils.ZeroVar:
            none_dict2[k[0]] = v
        else:
            if k[1] not in by_var2:
                by_var2[k[1]] = {}
            by_var2[k[1]][k[0]] = v

    for _vstr, _dict in by_var2.items():
        this_dict = {}
        for k, v in none_dict.items():
            if not best_effort_positive:
                this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in yz.schubmult_double(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
            else:
                this_dict = add_perm_dict(this_dict, {(k1, _vstr): expand(v1) * v for k1, v1 in yz.schubmult_double(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
        results = add_perm_dict(results, this_dict)

    none_dict, none_dict2 = sorted([none_dict, none_dict2], key=lambda x: -len(x.keys()))
    for k, v in none_dict2.items():
        results = add_perm_dict(results, {(k1, utils.NoneVar): v1 * v for k1, v1 in py.schubmult_py(none_dict, k).items()})

    return results


class DoubleSchubertAlgebraElement(Expr):
    """Algebra with sympy coefficients
    and a dict basis
    """

    _op_priority = 1e200
    # __slots__ = ("_dict", "_parent")
    _kind = NumberKind
    is_commutative = False
    # is_polynomial = True

    # default_coeff_var = "y"

    def __new__(cls, _dict, genset):
        return DoubleSchubertAlgebraElement.__xnew_cached__(cls, sympy.Dict({k: v for k, v in _dict.items() if expand(v) != 0}), genset)

    @staticmethod
    def __xnew__(_class, _dict, genset):
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if max([0] + [max([d + 1 for d in k[0].descents()] + [0]) for k in _dict.keys()]) > len(genset):
            raise NotEnoughGeneratorsError("Not enough generators")
        return Expr.__new__(_class, _dict, genset)

    @staticmethod
    @cache
    def __xnew_cached__(_class, _dict, genset):
        return DoubleSchubertAlgebraElement.__xnew__(_class, _dict, genset)

    @property
    def coeff_dict(self):
        return self.args[0]

    @property
    def genset(self):
        return self.args[1]

    def __hash__(self):
        return hash(tuple(self.args))

    # @property
    # def _add_handler(self):
    #     return SchubAdd

    # @property
    # def _mul_handler(self):
    #     return SchubMul

    # confers existing generating set
    def _from_dict(self, _dict):
        return DoubleSchubertAlgebraElement(_dict, self.genset)

    def divdiff(self, i):
        return self._from_dict({(k[0].swap(i - 1, i), k[1]): v for k, v in self.coeff_dict.items() if i - 1 in k[0].descents()})

    def mult_poly(self, poly):
        res_dict2 = {}
        # poly = self.genset[i + 1] - self.genset[i]
        for k, v in self.coeff_dict.items():
            if k[1] == utils.ZeroVar or k[1] == utils.NoneVar:
                dict2 = py.mult_poly_py({k[0]: v}, poly, self.genset)
            else:
                dict2 = yz.mult_poly_double({k[0]: v}, poly, self.genset, utils.poly_ring(k[1]))
            res_dict2 = add_perm_dict(res_dict2, {(k2, k[1]): v for k2, v in dict2.items()})

        return self._from_dict(res_dict2)

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
        return max([max([0, *list(k[0].descents(zero_indexed=False))]) for k in self.coeff_dict.keys()])

    def _eval_subs(self, old, new):
        result = 0
        if self.genset.index(old) != -1:
            # coproduct might help here
            logger.debug(f"I is the found {old=} {self.genset.index(old)=}")
            result = 0
            index = self.genset.index(old)
            mindex = self.max_index()
            logger.debug(f"{mindex=}")
            if mindex < index:
                return self
            # if already equal to the max index, we don't want to move it over
            perm = Permutation([]).swap(index - 1, mindex)  # index to max index + 1
            logger.debug(f"{mindex=}")
            logger.debug(f"{perm=}")
            transf = self.act(perm)
            #logger.debug(f"{transf=}")
            # logger.debug(f"{self.expand()=}")
            # logger.debug(f"{transf.expand().expand()=}")
            # transf2 = transf.coproduct([i for i in range(1,self.max_index()+1)],coeff_var=utils.NoneVar)
            # logger.debug(f"{transf2=}")
            # for (k1, k2), v in transf2.coeff_dict.items():
            #     result += self._from_dict({k1: v}) * (new**k2[0].inv)
            # don't want to go nuts
            # res_dict = {}
            for k, v in transf.coeff_dict.items():
                perm = k[0]
                coeff_var = k[1]
                coeff_gens = utils.poly_ring(coeff_var)
                # cached mul_poly
                L = schub_lib.pull_out_var(mindex+1, perm)
                #logger.debug(f"{perm=} {L=}")
                for index_list, new_perm in L:
                    result += self._from_dict({(new_perm, k[1]): v}).mult_poly(sympy.prod([(new - coeff_gens[index2]) for index2 in index_list]))
            return result

        for k, v in self.coeff_dict.items():
            if k[1] == utils.ZeroVar or k[1] == utils.NoneVar:
                add_dict = {k: v.subs(old, new)}
            else:
                coeff_genset = utils.poly_ring(k[1])
                if coeff_genset.index(old) != -1:
                    genset_list = [coeff_genset[i] for i in range(len(coeff_genset))]
                    genset_list[coeff_genset.index(old)] = 0
                    custom_genset = CustomGeneratingSet(genset_list)
                    new_add_dict = {k2: sympify(v2).subs(old, new) for k2, v2 in yz.schubmult_double({(): v}, k[0], custom_genset, coeff_genset).items()}  # remove the variable
                    add_dict = {}
                    for k3, v3 in new_add_dict.items():
                        # convert back to coeff_genset
                        to_add_dict = {(k4, k[1]): v4 for k4, v4 in yz.schubmult_double({(): v3}, k3, coeff_genset, custom_genset).items()}
                        add_dict = add_perm_dict(add_dict, to_add_dict)
                else:
                    add_dict = {k: sympify(v).subs(old, new)}
            for k5, v5 in add_dict.items():
                if any(self.genset.index(s) != -1 for s in sympify(v5).free_symbols):
                    result += self._from_dict({k5: 1}).mult_poly(v5)
                else:
                    result += self._from_dict({k5: v5})
            # check correct, change vars to zeroed coeff var for coeff
        return result

        # for k, v in self.coeff_dict.items():
        #     # can permute it to the end and substitute
        #     perm = k[0]
        #     coeff_var = k[1]

    @property
    def free_symbols(self):
        ret = set()
        for k, v in self.coeff_dict.items():
            ret.update(v.free_symbols)
            perm = k[0]
            coeff_var = k[1]
            if len(perm.descents()) > 0:
                ret.update([self.genset[i] for i in range(1, max(perm.descents()) + 2)])
            if coeff_var != utils.NoneVar and coeff_var != utils.ZeroVar:
                genset2 = utils.poly_ring(coeff_var)
                perm2 = ~perm
                if len(perm2.descents()) > 0:
                    ret.update([genset2[i] for i in range(1, max(perm2.descents()) + 2)])
        return ret

    # def _eval_Eq(self, other):
    #     # this will prevent sympy from acting like an idiot
    #     return self.__eq__(other)

    # def _eval_subs(self, old, new):
    #     b_old = sympify(old)
    #     b_new = sympify(new)
    #     result = {}
    #     stuff_to_do = False
    #     lots_of_stuff_to_do = False
    #     if b_new in utils.poly_ring(self._base_var):
    #         stuff_to_do = True
    #     if b_old in utils.poly_ring(self._base_var):
    #         lots_of_stuff_to_do = True
    #     for k, v in self.coeff_dict.items():
    #         if lots_of_stuff_to_do:
    #             poley = sympify(self._from_dict({k: 1}).change_vars(0).expand() * v)
    #             if b_old in poley.free_symbols:
    #                 poley = poley.subs(b_old, b_new)
    #                 new_dict = yz.mult_poly_double({(1, 2): 1}, poley, utils.poly_ring(self._base_var), utils.poly_ring(k[1]))
    #                 new_p = {(koifle, k[1]): voifle for koifle, voifle in new_dict.items()}
    #                 result = add_perm_dict(result, new_p)
    #         elif stuff_to_do:
    #             this_p = self._from_dict({k: v}).change_vars(0)
    #             for kkk, vvv in this_p.coeff_dict.items():
    #                 vvvv = sympify(vvv).subs(b_old, b_new)
    #                 if b_new in sympify(vvvv).free_symbols:
    #                     s_dict = {kkk[0]: 1}
    #                     r_dict = py.mult_poly_py(s_dict, vvvv, utils.poly_ring(self._base_var))
    #                 else:
    #                     r_dict = {kkk[0]: vvvv}
    #                 r_dict = {(kk, 0): voif for kk, voif in r_dict.items()}
    #                 new_p = self._from_dict(r_dict).change_vars(k[1])
    #                 result = add_perm_dict(result, new_p.coeff_dict)
    #         else:
    #             result[k] = result.get(k, 0) + sympify(v).subs(b_old, b_new)
    #     return self._from_dict(result)

    @cache
    def _cached_sympystr(self, printer):
        return printer.doprint(
            sympy.Add(
                *[
                    self.coeff_dict[k] if k[0] == Permutation([]) else sympy.Mul(self.coeff_dict[k], DSchubPoly(k, self.genset))
                    for k in sorted(self.coeff_dict.keys(), key=lambda bob: (inv(bob[0]), str(bob[1]), *bob[0]))
                ],
            ),
        )

    def _sympystr(self, printer):
        return self._cached_sympystr(printer)

    # def _eval_simplify(self, *args, measure, **kwargs):
    #     return self._from_dict({k: sympify(sympy.simplify(v, *args, measure=measure, **kwargs)) for k, v in self.coeff_dict.items()})

    def __add__(self, other):
        # if isinstance(self)
        # logger.debug(f"{type(other)=} {self.genset=}")
        try:
            other = DSx(other, genset=self.genset)
        except Exception:
            logger.debug(f"{other=} {list(self.genset)=}")
            return sympify(other) + self.as_polynomial()
        return self._from_dict(add_perm_dict(self.coeff_dict, other.coeff_dict))

    def __radd__(self, other):
        logger.debug(f"{type(other)=}")
        try:
            other = DSx(other, genset=self.genset)
        except Exception:
            logger.debug(f"{other=} {list(self.genset)=}")
            return sympify(other) + self.as_polynomial()
        return self._from_dict(add_perm_dict(other.coeff_dict, self.coeff_dict))

    def __sub__(self, other):
        logger.debug(f"{type(other)=}")
        try:
            other = DSx(other, genset=self.genset)
        except Exception:
            logger.debug(f"{other=} {list(self.genset)=}")
            return self.as_polynomial() - sympify(other)
        double_dict = add_perm_dict(self.coeff_dict, {k: -v for k, v in other.coeff_dict.items()})
        return self._from_dict(double_dict)

    def __rsub__(self, other):
        logger.debug(f"{type(other)=}")
        try:
            other = DSx(other, genset=self.genset)
        except Exception:
            logger.debug(f"{other=} {list(self.genset)=}")
            return sympify(other) - self.as_polynomial()
        double_dict = add_perm_dict(other.coeff_dict, {k: -v for k, v in self.coeff_dict.items()})
        return self._from_dict(double_dict)

    def __neg__(self):
        elem = self
        if self.is_Add or self.is_Mul:
            elem = self.doit()
        double_dict = {k: -sympify(v) for k, v in elem.coeff_dict.items()}
        return self._from_dict(double_dict)

    def __mul__(self, other):
        logger.debug(f"{type(other)=}")
        try:
            other = DSx(other, genset=self.genset)
        except Exception:
            logger.debug(f"{other=} {list(self.genset)=}")
            return self.as_polynomial() * sympify(other)
        return self._from_dict(_mul_schub_dicts(self.coeff_dict, other.coeff_dict))

    def __rmul__(self, other):
        logger.debug(f"{type(other)=}")
        try:
            other = DSx(other, genset=self.genset)
        except Exception:
            logger.debug(f"{other=} {list(self.genset)=}")
            return sympify(other) * self.as_polynomial()
        return self._from_dict(_mul_schub_dicts(other.coeff_dict, self.coeff_dict))

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

    # def __eq__(self, other):
    #     if self.is_Add or self.is_Mul:
    #         return self.doit().equals(other)
    #     cv = "y"
    #     elem1 = self
    #     elem2 = other

    #     if not elem1.test_equality(elem2):
    #         elem1_o = elem1.change_vars(cv)
    #         elem2_o = elem2.change_vars(cv)
    #         return elem1_o.test_equality(elem2_o)
    #     return True
    # assert all([k[1] == cv for k in elem1.coeff_dict.keys()])
    # assert all([k[1] == cv for k in elem2.coeff_dict.keys()])

    # def __str__(self):
    #     pieces = []
    #     keys = list(self.coeff_dict.keys())
    #     for k in sorted(keys, key=lambda b: (inv(b[0]), b[1], *b[0])):
    #         v = self.coeff_dict[k]
    #         dvar = "D"
    #         if sympy.expand(v) != 0:
    #             pieces += [
    #                 sympy.Mul(
    #                     v,
    #                     sympy.Symbol(
    #                         f"{dvar}S{DSx._base_var}({list(k[0])}, {_varstr(k[1])})",
    #                         commutative=False,
    #                     )
    #                     if k[0] != Permutation([])
    #                     else 1,
    #                 ),
    #             ]
    #     return sympy.sstr(sympy.Add(*pieces, evaluate=False), order="none")

    # def __repr__(self):
    #     return str(self)
    # def _as_ordered_terms(self, *args, **kwargs):

    @cache
    def change_vars(self, cv):
        result = {}
        for k, v in self.coeff_dict.items():
            result = add_perm_dict(result, {k1: v1 * v for k1, v1 in cached_positive_product(Permutation([]), k[0], cv, k[1]).items()})
        return self._from_dict(result)

    def as_coefficients_dict(self):
        return self.coeff_dict

    def expand(self, *args, **kwargs):  # noqa: ARG002
        return sympy.sympify(expand(sympify(self.as_polynomial())))

    def coproduct(self, indices, coeff_var="y", gname1=None, gname2=None):
        result_dict = {}
        if gname1 is None:
            gname1 = f"{self.genset.label}_A"
        if gname2 is None:
            gname2 = f"{self.genset.label}_B"
        gens2 = MaskedGeneratingSet(self.genset, indices)
        logger.debug(f"{indices=}")
        gens1 = gens2.complement()
        logger.debug(f"{gens1.index_mask=}")
        logger.debug(f"{list(gens1)=}")
        logger.debug(f"{gens2.index_mask=}")
        logger.debug(f"{list(gens2)=}")
        gens1.set_label(gname1)
        gens2.set_label(gname2)
        for k, v in self.coeff_dict.items():
            key = k[0]
            var_str = k[1]
            # print(f"{var_str=}")
            # print(f"{coeff_var=}")
            if var_str in (utils.NoneVar, utils.ZeroVar) and coeff_var in (utils.NoneVar, utils.ZeroVar):
                coprod_dict = py.schub_coprod_py(key, indices)
            else:
                coprod_dict = yz.schub_coprod_double(key, indices, utils.poly_ring(var_str), utils.poly_ring(coeff_var))
            # print(f"{coprod_dict=}")
            result_dict = add_perm_dict(result_dict, {((k1, var_str), (k2, coeff_var)): v * v2 for (k1, k2), v2 in coprod_dict.items()})
        basis = tsr.TensorAlgebraBasis(DoubleSchubertAlgebraElement_basis(gens1), DoubleSchubertAlgebraElement_basis(gens2))
        return basis._from_dict(result_dict)

    def legacy_coproduct(self, indices, coeff_var="y", gname1=None, gname2=None):
        result_dict = {}
        if gname1 is None:
            gname1 = f"{self.genset.label}_A"
        if gname2 is None:
            gname2 = f"{self.genset.label}_B"
        gens2 = MaskedGeneratingSet(self.genset, indices)
        logger.debug(f"{indices=}")
        gens1 = gens2.complement()
        logger.debug(f"{gens1.index_mask=}")
        logger.debug(f"{list(gens1)=}")
        logger.debug(f"{gens2.index_mask=}")
        logger.debug(f"{list(gens2)=}")
        gens1.set_label(gname1)
        gens2.set_label(gname2)
        for k, v in self.coeff_dict.items():
            key = k[0]
            var_str = k[1]
            # print(f"{var_str=}")
            # print(f"{coeff_var=}")
            if var_str in (utils.NoneVar, utils.ZeroVar) and coeff_var in (utils.NoneVar, utils.ZeroVar):
                coprod_dict = py.schub_coprod_py(key, indices)
            else:
                coprod_dict = yz.schub_coprod_double(key, indices, utils.poly_ring(var_str), utils.poly_ring(coeff_var))
            # print(f"{coprod_dict=}")
            result_dict = add_perm_dict(result_dict, {((k1, var_str), (k2, coeff_var)): v for (k1, k2), v in coprod_dict.items()})
        result_list = []
        for ktuple, v in result_dict.items():
            logger.debug(f"{ktuple=}")
            A = DSchubPoly(ktuple[0], gens1)
            B = DSchubPoly(ktuple[1], gens2)
            # if A.perm == Permutation([]):
            #     if B.perm == Permutation([]):
            #         result_list += [v]
            #     result_list += [sympy.Mul(v, B)]
            # elif B.perm == Permutation([]):
            #     result_list += [sympy.Mul(v, A)]
            # else:
            #     logger.debug(f"{A=} {B=}")
            #     result_list += [sympy.Mul(v, A, B, evaluate=False)]
            logger.debug(f"{A=} {B=}")
            result_list += [sympy.Mul(sympy.sympify(v), A, B)]
        return sympy.Add(*result_list)

        # # will not allow zeros

        # return {k: v for k, v in self.coeff_dict.items() if expand(v) != 0}

    # def normalize_coefficients(self, coeff_var):
    #     return DSx([1, 2], coeff_var) * self

    # def expand(self, *_a, **_):
    #     if isinstance(self, SchubAdd):
    #         return self.doit().expand()
    #     if isinstance(self, SchubMul):
    #         return self.doit().expand()
    #     return expand(Add(*[v * schubpoly(k[0], self.genset, utils.poly_ring(k[1])) for k, v in self.coeff_dict.items()]))

    @cached_property
    def max_gens(self):
        return max([max(k[0].descents()) for k in self.coeff_dict.keys()])

    def as_polynomial(self):
        return sympy.sympify(Add(*[v * xreplace_genvars(cached_schubpoly(k[0]), self.genset, utils.poly_ring(k[1])) for k, v in self.coeff_dict.items()]))


# Atomic Schubert polynomial
class DSchubPoly(DoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, genset):
        return DSchubPoly.__xnew_cached__(cls, k, genset)

    @staticmethod
    def __xnew__(_class, k, genset):
        _coeff_dict = sympy.Dict({(Permutation(k[0]), k[1]): 1})
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        obj = DoubleSchubertAlgebraElement.__new__(_class, _coeff_dict, genset)
        obj._key = k
        obj._genset = genset
        obj._coeff_dict = _coeff_dict
        return obj

    @property
    def coeff_dict(self):
        return self._coeff_dict

    @property
    def perm(self):
        return self._key[0]

    @property
    def args(self):
        return (sympy.Tuple(*self._key), self._genset)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset):
        return DSchubPoly.__xnew__(_class, k, genset)

    def _sympystr(self, printer):
        if self._key[0] == Permutation([]):
            return printer.doprint(1)
        if self._key[1] == 0 or self._key[1] == utils.NoneVar:
            return printer.doprint(f"S{self.genset.label}({list(self._key[0])})")
        return printer.doprint(f"DS{self.genset.label}({list(self._key[0])}, {_varstr(self._key[1])})")


# None is faster to store
class DoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset):
        return Basic.__new__(cls, genset)

    @property
    def genset(self):
        return self.args[0]

    def _from_dict(self, _dict):
        return DoubleSchubertAlgebraElement(_dict, self.genset)

    def __call__(self, x, cv=None, genset=None):
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        logger.debug(f"{type(x)=}")
        # if isinstance(x, Mul) or isinstance(x, Add):
        #     raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            if cv is None:
                cv = "y"
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = DoubleSchubertAlgebraElement({(p_x, cv): 1}, genset)
        elif isinstance(x, Permutation):
            if cv is None:
                cv = "y"
            if max([0, *list(x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = DoubleSchubertAlgebraElement({(x, cv): 1}, genset)

        elif isinstance(x, DoubleSchubertAlgebraElement):
            logger.debug("Line record")
            if x.is_Add or x.is_Mul:
                return x.doit()
            if x.genset == genset:
                return x
            raise ValueError("Different generating set")
        # poly
        # elif isinstance(x, sympy.poly):

        else:
            logger.debug(f"{x=}")
            x = sympify(x)
            if cv is None or cv == utils.NoneVar:
                cv = utils.NoneVar
                result = py.mult_poly_py({Permutation([]): 1}, x, genset)
            else:
                result = yz.mult_poly_double({Permutation([]): 1}, x, genset, utils.poly_ring(cv))
            elem = DoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()}, genset)
            logger.debug(f"Returning {elem=}")
        return elem


# def _do_schub_mul(a, b):
#     A = DSx(a)
#     B = DSx(b)
#     return self._from_dict(_mul_schub_dicts(A.coeff_dict, B.coeff_dict))


# def _do_schub_add(a, b):
#     A = DSx(a)
#     B = DSx(b)
#     return self._from_dict(add_perm_dict(A.coeff_dict, B.coeff_dict))


def get_postprocessor(cls):
    if cls is Mul:
        return lambda expr: SchubMul(*expr.args)  # .doit()
    if cls is Add:
        return lambda expr: SchubAdd(*expr.args)  # .doit()
    return None


Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
    "Mul": [get_postprocessor(Mul)],
    "Add": [get_postprocessor(Add)],
}

# add.register_handlerclass((Expr, SchubAdd), SchubAdd)
# mul.register_handlerclass((Expr, SchubMul), SchubMul)


DoubleSchubertPolynomial = DoubleSchubertAlgebraElement


class SchubAdd(Add):
    is_Add = True

    def __new__(cls, *args, evaluate=False, _sympify=True, **_):
        obj = Add.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)
        obj._args = args
        if evaluate:
            return obj.doit()
        return obj

    @property
    def args(self):
        return self._args

    def doit(self):
        ret = self.args[0]
        logger.debug(f"ADD {self.args=}")
        for arg in self.args[1:]:
            logger.debug(f"{arg=} {type(arg)=}")
            logger.debug(f"{ret=} {type(ret)=}")
            ret += sympy.expand(arg)
        return ret

    def _sympystr(self, printer):
        return printer._print_Add(self)

    def expand(self, *_, **__):
        logger.debug(f"Pringles {self.args=}")
        return SchubAdd(*[sympy.expand(arg) for arg in self.args]).doit()


class SchubMul(Mul):
    is_Mul = True

    def __new__(cls, *args, evaluate=False, _sympify=True, **_):
        # args, a, b = Mul.flatten(list(args))
        # if len(args) == 0:
        #     return 1
        obj = Mul.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)
        obj._args = args
        if evaluate:
            return obj.doit()
        return obj

    @property
    def args(self):
        return self._args

    def doit(self):
        ret = self.args[0]
        logger.debug(f"MUL {self.args=}")
        for arg in self.args[1:]:
            logger.debug(f"{arg=} {type(arg)=}")
            ret *= sympy.expand(arg)
        return ret

    def _sympystr(self, printer):
        return printer._print_Mul(self)

    def __neg__(self):
        return SchubMul(sympy.Integer(-1), self)

    def _eval_expand_mul(self, *_, **__):
        logger.debug(f"Pringles {self.args=}")
        return SchubMul(*[sympy.expand(arg) for arg in self.args]).doit()


# Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
#     "Mul": [get_postprocessor(Mul)],
#     "Add": [get_postprocessor(Add)],
# }

DSx = DoubleSchubertAlgebraElement_basis(GeneratingSet("x"))
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


def Sx(x):
    return DSx(x, utils.NoneVar)
