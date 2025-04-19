import multiprocessing
from functools import cache, cached_property

import sympy
from symengine import Add, S, Symbol, expand, sympify
from sympy import Basic
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.printing.str import StrPrinter

import schubmult.rings._quantum_schubert_polynomial_ring as qsr
import schubmult.rings._tensor_schub_ring as tsr
import schubmult.rings._utils as utils
import schubmult.schub_lib.double as yz

# from schubmult.poly_lib.schub_poly import pull_out_var
import schubmult.schub_lib.schub_lib as schub_lib
import schubmult.schub_lib.single as py
from schubmult.perm_lib import Permutation, inv
from schubmult.poly_lib.poly_lib import elem_sym_poly, xreplace_genvars
from schubmult.poly_lib.schub_poly import schubpoly_classical_from_elems, schubpoly_from_elems
from schubmult.poly_lib.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

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


def domul(t1, dict2):
    _vstr, kd, vd, basis, best_effort_positive = t1
    this_dict = {}
    for k, v in dict2:
        did_positive = False
        to_mul = v * vd
        if best_effort_positive:
            try:
                # logger.critical(f"{to_mul=} {kd=} {k=}")
                this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis.cached_positive_product(kd, k[0], _vstr, k[1]).items()})
                did_positive = True
            except Exception:
                # logger.debug("Failed to compute")
                did_positive = False
        if not did_positive:
            this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis.cached_product(kd, k[0], _vstr, k[1]).items()})
    return this_dict


def _mul_schub_dicts(dict1, dict2, basis, best_effort_positive=True):
    by_var = {}

    none_dict = {}
    for k, v in dict1.items():
        if k[1] == utils.NoneVar:  # or k[1] == utils.ZeroVar:
            none_dict[k[0]] = v
        else:
            if k[1] not in by_var:
                by_var[k[1]] = {}
            by_var[k[1]][k[0]] = v

    results = {}
    # import sys

    for _vstr, _dict in by_var.items():
        if BasisSchubertAlgebraElement.do_parallel:
            result_list = []
            mul_funcs = [(_vstr, kd, vd, basis, best_effort_positive) for kd, vd in _dict.items()]
            itemlist = list(dict2.items())

            def add_result(result):
                result_list.append(result)

            with multiprocessing.Pool(processes=6) as pool:
                for mulf in mul_funcs:
                    pool.apply_async(domul, args=(mulf, itemlist), callback=add_result)
                pool.close()
                pool.join()
            for res in result_list:
                results = add_perm_dict(results, res)
        else:
            this_dict = {}
            for k, v in dict2.items():
                for kd, vd in _dict.items():
                    did_positive = False
                    to_mul = v * vd
                    if best_effort_positive:
                        try:
                            # logger.critical(f"{to_mul=} {kd=} {k=}")
                            this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis.cached_positive_product(kd, k[0], _vstr, k[1]).items()})
                            did_positive = True
                        except Exception:
                            # logger.debug("Failed to compute")
                            did_positive = False
                    if not did_positive:
                        this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis.cached_product(kd, k[0], _vstr, k[1]).items()})
            results = add_perm_dict(results, this_dict)

    by_var2 = {}
    none_dict2 = {}
    for k, v in dict2.items():
        if k[1] == utils.NoneVar:  # or k[1] == utils.ZeroVar:
            none_dict2[k[0]] = v
        else:
            if k[1] not in by_var2:
                by_var2[k[1]] = {}
            by_var2[k[1]][k[0]] = v

    for _vstr, _dict in by_var2.items():
        this_dict = {}
        for k, v in none_dict.items():
            if not best_effort_positive:
                this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in basis.double_mul(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
            else:
                this_dict = add_perm_dict(this_dict, {(k1, _vstr): expand(v1) * v for k1, v1 in basis.double_mul(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
        results = add_perm_dict(results, this_dict)

    none_dict, none_dict2 = sorted([none_dict, none_dict2], key=lambda x: -len(x.keys()))
    for k, v in none_dict2.items():
        results = add_perm_dict(results, {(k1, utils.NoneVar): v1 * v for k1, v1 in basis.single_mul(none_dict, k).items()})

    return results


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
        # obj.prune()
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
            if k[1] == utils.ZeroVar or k[1] == utils.NoneVar:
                dict2 = self.basis.mult_poly_single({k[0]: v}, poly, self.genset)
            else:
                dict2 = self.basis.mult_poly_double({k[0]: v}, poly, self.genset, utils.poly_ring(k[1]))
            res_dict2 = add_perm_dict(res_dict2, {(k2, k[1]): v for k2, v in dict2.items()})
        logger.debug(f"{res_dict2=}")
        return self.basis._from_dict(res_dict2)

    # def _cached_sympystr(self, printer):
    #     return printer.doprint(
    #         sympy.Add(
    #             *[
    #                 (self.coeff_dict[k] if k[0] == Permutation([]) else sympy.Mul(self.coeff_dict[k], self.basis.single_element_class(k, self.basis)))
    #                 for k in sorted(self.coeff_dict.keys(), key=lambda bob: (inv(bob[0]), str(bob[1]), *bob[0]))
    #             ],
    #         ),
    #     )

    def _sympystr(self, printer):
        return printer._print_Add(self)
        # return self._cached_sympystr(printer)

    def as_terms(self):
        if len(self.coeff_dict.keys()) == 0:
            return [sympy.sympify(S.Zero)]
        return [
            (sympy.sympify(self.coeff_dict[k]) if k[0] == Permutation([]) else sympy.Mul(self.coeff_dict[k], self.basis.single_element_class(k, self.basis)))
            for k in sorted(self.coeff_dict.keys(), key=lambda bob: (inv(bob[0]), str(bob[1]), *bob[0]))
        ]

    def as_ordered_terms(self, *_, **__):
        return self.as_terms()

    # def _eval_simplify(self, *args, measure, **kwargs):
    #     return self.basis._from_dict({k: sympify(sympy.simplify(v, *args, measure=measure, **kwargs)) for k, v in self.coeff_dict.items()})

    # def __iadd__(self, other):
    #     return self.__add__(other)

    def __add__(self, other):
        # if isinstance(self)
        # # logger.debug(f"{type(other)=} {self.genset=}")
        try:
            other = self.basis(other)
        except Exception:
            return sympify(other) + self.as_polynomial()
        return self.basis._from_dict(add_perm_dict(self.coeff_dict, other.coeff_dict))

    def __radd__(self, other):
        # logger.debug(f"{type(other)=}")
        try:
            other = self.basis(other)
        except Exception:
            # logger.debug(f"{other=} {list(self.genset)=}")
            return sympify(other) + self.as_polynomial()
        return self.basis._from_dict(add_perm_dict(other.coeff_dict, self.coeff_dict))

    def __sub__(self, other):
        # logger.debug(f"{type(other)=}")
        try:
            other = self.basis(other)
        except Exception:
            # logger.debug(f"{other=} {list(self.genset)=}")
            return self.as_polynomial() - sympify(other)
        double_dict = add_perm_dict(self.coeff_dict, {k: -v for k, v in other.coeff_dict.items()})
        return self.basis._from_dict(double_dict)

    def __rsub__(self, other):
        # logger.debug(f"{type(other)=}")
        try:
            other = self.basis(other)
        except Exception:
            # logger.debug(f"{other=} {list(self.genset)=}")
            return sympify(other) - self.as_polynomial()
        double_dict = add_perm_dict(other.coeff_dict, {k: -v for k, v in self.coeff_dict.items()})
        return self.basis._from_dict(double_dict)

    def __neg__(self):
        elem = self
        if self.is_Add or self.is_Mul:
            elem = self.doit()
        double_dict = {k: -sympify(v) for k, v in elem.coeff_dict.items()}
        return self.basis._from_dict(double_dict)

    def __mul__(self, other):
        try:
            o = sympify(other)
            return self.mult_poly(o)
        except Exception:
            try:
                other = self.basis(other)
                return self.basis._from_dict(_mul_schub_dicts(self.coeff_dict, other.coeff_dict, self.basis))
            except Exception:
                return self.as_polynomial() * sympify(other)

    def __rmul__(self, other):
        # logger.debug(f"{type(other)=}")
        try:
            o = sympify(other)
            return self.mult_poly(o)
        except Exception:
            try:
                other = self.basis(other)
                return self.basis._from_dict(_mul_schub_dicts(other.coeff_dict, self.coeff_dict, self.basis))
            except Exception:
                return self.as_polynomial() * sympify(other)

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
        # result = {}
        # fix
        # for k, v in self.coeff_dict.items():
        #     result = add_perm_dict(result, {k1: v1 * v for k1, v1 in self.basis.cached_positive_product(Permutation([]), k[0], cv, k[1]).items()})
        # # result = {(k, cv): v for k, v in self.basis.mul_double(Permutation([]),utils.poly_ring)}
        # return self.basis._from_dict(result)
        return self.basis([], cv) * self

    def as_coefficients_dict(self):
        return sympy.Dict({self.basis.single_element_class(k, self.basis): sympy.sympify(v) for k, v in self.coeff_dict.items()})

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
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

    # default_coeff_var = "y"

    def __new__(cls, _dict, basis):
        return BasisSchubertAlgebraElement.__new__(cls, _dict, basis)

    def divdiff(self, i):
        return self.basis._from_dict({(k[0].swap(i - 1, i), k[1]): v for k, v in self.coeff_dict.items() if i - 1 in k[0].descents()})

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
            # logger.debug(f"{mindex=}")
            # logger.debug(f"{perm=}")
            transf = self.act(perm)
            # # logger.debug(f"{transf=}")
            # # logger.debug(f"{self.expand()=}")
            # # logger.debug(f"{transf.expand().expand()=}")
            # transf2 = transf.coproduct([i for i in range(1,self.max_index()+1)],coeff_var=utils.NoneVar)
            # # logger.debug(f"{transf2=}")
            # for (k1, k2), v in transf2.coeff_dict.items():
            #     result += self.basis._from_dict({k1: v}) * (new**k2[0].inv)
            # don't want to go nuts
            # res_dict = {}
            for k, v in transf.coeff_dict.items():
                perm = k[0]
                coeff_var = k[1]
                coeff_gens = utils.poly_ring(coeff_var)
                # cached mul_poly
                L = schub_lib.pull_out_var(mindex + 1, perm)
                # # logger.debug(f"{perm=} {L=}")
                for index_list, new_perm in L:
                    result += self.basis._from_dict({(new_perm, k[1]): v}).mult_poly(sympy.prod([(new - coeff_gens[index2]) for index2 in index_list]))
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
                    result += self.basis._from_dict({k5: 1}).mult_poly(v5)
                else:
                    result += self.basis._from_dict({k5: v5})
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
    #             poley = sympify(self.basis._from_dict({k: 1}).change_vars(0).expand() * v)
    #             if b_old in poley.free_symbols:
    #                 poley = poley.subs(b_old, b_new)
    #                 new_dict = yz.mult_poly_double({(1, 2): 1}, poley, utils.poly_ring(self._base_var), utils.poly_ring(k[1]))
    #                 new_p = {(koifle, k[1]): voifle for koifle, voifle in new_dict.items()}
    #                 result = add_perm_dict(result, new_p)
    #         elif stuff_to_do:
    #             this_p = self.basis._from_dict({k: v}).change_vars(0)
    #             for kkk, vvv in this_p.coeff_dict.items():
    #                 vvvv = sympify(vvv).subs(b_old, b_new)
    #                 if b_new in sympify(vvvv).free_symbols:
    #                     s_dict = {kkk[0]: 1}
    #                     r_dict = py.mult_poly_py(s_dict, vvvv, utils.poly_ring(self._base_var))
    #                 else:
    #                     r_dict = {kkk[0]: vvvv}
    #                 r_dict = {(kk, 0): voif for kk, voif in r_dict.items()}
    #                 new_p = self.basis._from_dict(r_dict).change_vars(k[1])
    #                 result = add_perm_dict(result, new_p.coeff_dict)
    #         else:
    #             result[k] = result.get(k, 0) + sympify(v).subs(b_old, b_new)
    #     return self.basis._from_dict(result)

    def coproduct(self, indices, coeff_var="y", gname1=None, gname2=None):
        result_dict = {}
        if gname1 is None:
            gname1 = f"{self.genset.label}_A"
        if gname2 is None:
            gname2 = f"{self.genset.label}_B"
        gens2 = MaskedGeneratingSet(self.genset, indices)
        # logger.debug(f"{indices=}")
        gens1 = gens2.complement()
        # logger.debug(f"{gens1.index_mask=}")
        # logger.debug(f"{list(gens1)=}")
        # logger.debug(f"{gens2.index_mask=}")
        # logger.debug(f"{list(gens2)=}")
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

    # def normalize_coefficients(self, coeff_var):
    #     return DSx([1, 2], coeff_var) * self

    # def expand(self, *_a, **_):
    #     if isinstance(self, SchubAdd):
    #         return self.doit().expand()
    #     if isinstance(self, SchubMul):
    #         return self.doit().expand()
    #     return expand(Add(*[v * schubpoly(k[0], self.genset, utils.poly_ring(k[1])) for k, v in self.coeff_dict.items()]))

    def in_SEM_basis(self):
        result = S.Zero
        for k, v in self.coeff_dict.items():
            result += v * schubpoly_from_elems(k[0], self.genset, utils.poly_ring(k[1]), elem_func=self.basis.symbol_elem_func)
        # print(f"{result=}")
        # gens = []
        # for k in range(1, 10):
        #     gens += [sympy.Symbol(f"e_{p}_{k}") for p in range(1,k+1)]
        # #print(f"{gens=}")
        # ply = sympy.poly(sympy.sympify(expand(result)), *gens)
        # floss = 0
        # for m, c in ply.as_dict().items():
        #     floss += c * sympy.prod([gens[i]**m[i] for i in range(len(m))])
        # return floss
        return result

    @cached_property
    def max_gens(self):
        return max([max(k[0].descents()) for k in self.coeff_dict.keys()])


# Atomic Schubert polynomial
class DSchubPoly(DoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, basis):
        return DSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        _coeff_dict = sympy.Dict({(Permutation(k[0]), k[1]): 1})
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
        return self._key[0]

    @property
    def args(self):
        return (sympy.Tuple(*self._key), self._basis)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, basis):
        return DSchubPoly.__xnew__(_class, k, basis)

    def _sympystr(self, printer):
        if self._key[0] == Permutation([]):
            return printer.doprint(1)
        if self._key[1] == 0 or self._key[1] == utils.NoneVar:
            return printer.doprint(f"S{self.genset.label}({printer.doprint(self._key[0])})")
        return printer.doprint(f"DS{self.genset.label}({printer.doprint(self._key[0])}, {_varstr(self._key[1])})")


# def elem_func(p, k, vx, vy):
#     return DSx(elem_func_q(p, k, vx, vy), "y")

# A = schubpoly_from_elems([4,1,3,2], DSx.genset, poly_ring("y"),elem_func)


# None is faster to store
class DoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset):
        return Basic.__new__(cls, genset)

    @property
    def symbol_elem_func(self):
        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return 1
            if p < 0 or p > k:
                return 0
            return sympy.Add(*[(Symbol(f"e_{p - i}_{k}") if p - i > 0 else 1) * elem_sym_poly(i, k + 1 - p, [-v for v in varl2], [0 for a in varl1]) for i in range(p + 1)])

        return elem_func

    # def in_SEM_basis(self, elem):
    #     return

    @property
    def genset(self):
        return self.args[0]

    def _from_dict(self, _dict):
        return DoubleSchubertAlgebraElement(_dict, self)

    @property
    def single_element_class(self):
        return DSchubPoly

    def in_quantum_basis(self, elem):
        result = S.Zero
        for k, v in elem.coeff_dict.items():
            result += v * self.quantum_schubpoly(k[0], k[1])
        return result

    def in_classical_basis(self, elem):
        return elem

    @cache
    def quantum_schubpoly(self, perm, coeff_var="y"):
        return schubpoly_classical_from_elems(perm, self.genset, utils.poly_ring(coeff_var), self.quantum_elem_func(coeff_var))

    @cache
    def cached_product(self, u, v, va, vb):
        return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, va, vb):
        return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_generic_partial_posify(u, v).items()}

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

    def quantum_elem_func(self, coeff_var):
        basis = qsr.QuantumDoubleSchubertAlgebraElement_basis(self.genset)

        def elem_func(p, k, varl1, varl2, xstart=0, ystart=0):
            if p > k:
                return basis(0, coeff_var)
            if p == 0:
                return basis([], coeff_var)
            if p == 1:
                res = basis(varl1[xstart] - varl2[ystart], coeff_var)
                for i in range(1, k):
                    res += basis(varl1[xstart + i] - varl2[ystart + i], coeff_var)
                return res
            if p == k:
                res = basis((varl1[xstart] - varl2[ystart]) * (varl1[xstart + 1] - varl2[ystart]), coeff_var)
                for i in range(2, k):
                    res *= basis(varl1[i + xstart] - varl2[ystart], coeff_var)
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
            logger.debug(f"{res=}")
            return res

        return elem_func

    # @cache
    # def cached_schubpoly_oink(self, u):
    #     return yz.schubpoly(u)

    @cache
    def cached_schubpoly(self, k):
        # return yz.schubpoly(u)
        return schubpoly_classical_from_elems(k[0], self.genset, utils.poly_ring(k[1]), elem_func=elem_sym_poly)

    def __call__(self, x, cv=None):
        # print(f"frivol {x=} {cv=}")
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        # logger.debug(f"{type(x)=}")
        # if isinstance(x, Mul) or isinstance(x, Add):
        #     raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            if cv is None:
                cv = "y"
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self._from_dict({(p_x, cv): 1})
        elif isinstance(x, Permutation):
            if cv is None:
                cv = "y"
            if max([0, *list(x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self._from_dict({(x, cv): 1})

        elif isinstance(x, DoubleSchubertAlgebraElement):
            # logger.debug("Line record")
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
        #         # logger.debug(f"Didn't find any gens in {x=}")
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
            # logger.debug(f"{x=}")
            x = sympify(x)
            if cv is None or cv == utils.NoneVar:
                cv = utils.NoneVar
                result = py.mult_poly_py({Permutation([]): 1}, x, genset)
            else:
                # print("splinterfish")
                result = yz.mult_poly_double({Permutation([]): 1}, x, genset, utils.poly_ring(cv))
            elem = self._from_dict({(k, cv): v for k, v in result.items()})
            logger.debug(f"Returning {elem=}")
        return elem


# def _do_schub_mul(a, b):
#     A = DSx(a)
#     B = DSx(b)
#     return self.basis._from_dict(_mul_schub_dicts(A.coeff_dict, B.coeff_dict))


# def _do_schub_add(a, b):
#     A = DSx(a)
#     B = DSx(b)
#     return self.basis._from_dict(add_perm_dict(A.coeff_dict, B.coeff_dict))


# def get_postprocessor(cls):
#     if cls is Mul:
#         return lambda expr: SchubMul(*expr.args)  # .doit()
#     if cls is Add:
#         return lambda expr: SchubAdd(*expr.args)  # .doit()
#     return None


# Basic._constructor_postprocessor_mapping[BasisSchubertAlgebraElement] = {
#     "Mul": [get_postprocessor(Mul)],
#     "Add": [get_postprocessor(Add)],
# }

# add.register_handlerclass((Expr, SchubAdd), SchubAdd)
# mul.register_handlerclass((Expr, SchubMul), SchubMul)


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
#         # logger.debug(f"ADD {self.args=}")
#         for arg in self.args[1:]:
#             # logger.debug(f"{arg=} {type(arg)=}")
#             # logger.debug(f"{ret=} {type(ret)=}")
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
#         # logger.debug(f"MUL {self.args=}")
#         for arg in self.args[1:]:
#             # logger.debug(f"{arg=} {type(arg)=}")
#             ret *= sympy.expand(arg)
#         return ret

#     def _sympystr(self, printer):
#         return printer._print_Mul(self)

#     def __neg__(self):
#         return SchubMul(sympy.Integer(-1), self)

#     def _eval_expand_mul(self, *_, **__):
#         # logger.debug(f"Pringles {self.args=}")
#         return SchubMul(*[sympy.expand(arg) for arg in self.args]).doit()


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


# def Sx(x):
#     return DSx(x, utils.NoneVar)


class SchubertAlgebraElement_basis(DoubleSchubertAlgebraElement_basis):
    def __new__(cls, genset):
        return DoubleSchubertAlgebraElement_basis.__new__(cls, genset)

    def _from_single_dict(self, _dict):
        return DoubleSchubertAlgebraElement({(k, utils.NoneVar): v for k, v in _dict.items()}, self)

    def __call__(self, x):
        genset = self.genset
        # logger.debug(f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self._from_single_dict({Permutation(x): 1})
        elif isinstance(x, Permutation):
            elem = self._from_single_dict({x: 1})
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
            elem = self._from_single_dict(result)
        return elem


Sx = SchubertAlgebraElement_basis(GeneratingSet("x"))

ybas = SchubertAlgebraElement_basis(GeneratingSet("y"))
