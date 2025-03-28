from functools import cache

import sympy
from symengine import sympify
from sympy import Add, Basic, Mul
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.printing.str import StrPrinter

import schubmult.rings._utils as utils

# import schubmult.rings._schubert_polynomial_ring as spr
import schubmult.schubmult_double as yz
import schubmult.schubmult_py as py

# from sympy.core.operations import AssocOp
from schubmult.perm_lib import add_perm_dict, inv, permtrim

# import utils.NoneVar, utils.ZeroVar, utils.poly_ring

_def_printer = StrPrinter({"order": "none"})
# numpy arrays
# sympy parsing
# quantum
# none is faster
# coeff_var = None
# singleton
# class utils.NoneVarType(sympy.Basic):
#     def __new__(cls, *args):
#         return utils.NoneVarType.__xnew_cached__(cls, args)

#     def __xnew__(cls, *Args):
#         obj = sympy.Basic.__new__(*args)
#         return obj


# def is_nonevar(v):
#     return isinstance(v, utils.NoneVarType)

# _utils.NoneVar = utils.NoneVarType()

# sympy.init_printing(order='none')

# class DSchubPoly(sympy.Symbol):
# as Add


class ExceptionThrower:
    def __bool__(self):
        raise Exception


def _varstr(v):
    if v == utils.NoneVar:
        return "NoneVar"
    if v == utils.ZeroVar:
        return "0"
    return f"'{v}'"


def _mul_schub_dicts(dict1, dict2):
    by_var = {}

    none_dict = {}
    for k, v in dict1.items():
        if k[1] == utils.NoneVar:
            none_dict[k[0]] = v
        else:
            if k[1] not in by_var:
                by_var[k[1]] = {}
            by_var[k[1]][k[0]] = v

    results = {}

    for _vstr, _dict in by_var.items():
        this_dict = {}
        for k, v in dict2.items():
            # print(f"{k=}")
            this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in yz.schubmult(_dict, k[0], utils.poly_ring(_vstr), utils.poly_ring(k[1])).items()})
        results.update(this_dict)

    by_var2 = {}
    none_dict2 = {}
    for k, v in dict2.items():
        if k[1] == utils.NoneVar:
            none_dict2[k[0]] = v
        else:
            if k[1] not in by_var2:
                by_var2[k[1]] = {}
            by_var2[k[1]][k[0]] = v

    for _vstr, _dict in by_var2.items():
        this_dict = {}
        for k, v in none_dict.items():
            this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in yz.schubmult(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
        results = add_perm_dict(results, this_dict)

    this_dict = {}
    for k, v in none_dict2.items():
        this_dict = add_perm_dict(this_dict, {(k1, utils.NoneVar): v1 * v for k1, v1 in py.schubmult(none_dict, k).items()})
    results.update(this_dict)

    return results


class DSchubSymbol(sympy.Symbol):
    def __new__(cls, base_var, k, *args, **kwargs):
        return DSchubSymbol.__xnew_cached__(base_var, k, *args, **kwargs)

    @classmethod
    def __xnew__(cls, base_var, k):
        if k[1] == 0 or k[1] == utils.NoneVar:
            obj = sympy.Symbol.__new__(cls, f"S{base_var}({list(k[0])})", commutative=False)
        else:
            obj = sympy.Symbol.__new__(cls, f"DS{base_var}({list(k[0])}, {_varstr(k[1])})", commutative=False)
        return obj

    @classmethod
    @cache
    def __xnew_cached__(cls, base_var, k):
        return DSchubSymbol.__xnew__(base_var, k)


class DoubleSchubertAlgebraElement(Expr):
    """Algebra with sympy coefficients
    and a dict basis
    """

    _op_priority = 1e200
    # __slots__ = ("_dict", "_parent")
    _kind = NumberKind
    is_Atom = True
    is_commutative = True
    # default_coeff_var = "y"

    def __new__(cls, _dict, parent, *args, **kwargs):
        # print(f"{_dict} {parent=} {args=} {kwargs=}")
        return DoubleSchubertAlgebraElement.__xnew_cached__(cls, sympy.Dict(_dict), parent, *args, **kwargs)

    def __hash__(self):
        return hash(self._doubledict)

    def _latex(self, printer):
        return self._sympystr(printer)

    def _sympyrepr(self, printer):
        return self._sympystr(printer)

    @property
    def _add_handler(self):
        return SchubAdd

    @property
    def _mul_handler(self):
        # print("profilating")
        return SchubMul

    def _sympystr(self, *args):  # noqa: ARG002
        printer = _def_printer
        if not self:
            return printer._print(0)
        # ring = self.ring
        # symbols = ring.symbols
        # ngens = ring.ngens
        # zm = ring.zero_monom
        _dict = self._doubledict
        keys = list(_dict.keys())
        # print("Freftoolnagababarmpy")
        pieces = []
        for k in sorted(keys, key=lambda b: (inv(b[0]), str(b[1]) if b[1] != utils.NoneVar else ".", *b[0])):
            v = _dict[k]
            if sympy.expand(v) != 0:
                pieces += [sympy.Mul(v, DSchubSymbol(self._parent._base_var, k))]

        return printer._print(sympy.Add(*pieces, evaluate=False))

    # def __repr__(self):

    #    return sympy.Add(*self.args, evaluate=False).__repr__()
    @classmethod
    def __xnew__(cls, _class, _dict, parent, *args, **kwargs):
        # pieces = []
        # keys = list(_dict.keys())
        # #print("Freftoolnagababarmpy")
        # for k in sorted(keys, key=lambda b: (inv(b[0]), str(b[1]) if b[1] != utils.NoneVar else str("."), *b[0])):
        #     v = _dict[k]
        #     dvar = "D"
        #     if sympy.expand(v) != 0:
        #         pieces += [
        #             sympy.Mul(v,DSchubSymbol(parent._base_var, k))
        #         ]
        # print(f"{pieces=}")
        # print(f"{args=} {kwargs=} {_dict=} {parent=}")
        # print(f"{args=} {kwargs=}")
        obj = Expr.__new__(_class, _dict, parent)
        # obj.make_args(pieces)
        # obj.args = pieces
        obj._doubledict = _dict
        obj._parent = parent
        return obj

    @classmethod
    @cache
    def __xnew_cached__(cls, _class, _dict, parent, *args, **kwargs):
        return DoubleSchubertAlgebraElement.__xnew__(_class, _dict, parent, *args, **kwargs)

    def _symengine_(self):
        return NotImplemented

    def _eval_simplify_(self):
        # print("Hey pretty baby")
        return self

    def _from_double_dict(self, _doubledict):
        return self.__class__.__new__(self.__class__, _doubledict, self._parent, double_dict=True)

    def __add__(self, other):
        # print("ASFJASJ")
        # print(f"addwinky {self.__class__=} {self=} {other=}")
        # print(f"flarfknockle {self._doubledict=} {self._parent(other)._doubledict=} {other=}")
        ret = add_perm_dict(self._doubledict, self._parent(other)._doubledict)
        return self._from_double_dict(ret)

    def __radd__(self, other):
        # print("is wiz doing radd")
        return self._from_double_dict(add_perm_dict(self._parent(other)._doubledict, self._doubledict))

    def __sub__(self, other):
        # print("ASFJAdsajdSJ")

        return self._from_double_dict(add_perm_dict(self._doubledict, {k: -v for k, v in self._parent(other)._doubledict.items()}))

    def __rsub__(self, other):
        # print("ASFJAdsajdSJ")
        return self._from_double_dict(add_perm_dict(self._parent(other)._doubledict, {k: -v for k, v in self._doubledict.items()}))

    def __neg__(self):
        # print("ASFJAdsajdSJ")

        return self._from_double_dict({k: -v for k, v in self._doubledict.items()})

    def __mul__(self, other):
        # print("ASFJAdsajdSJ")
        # print(f"mulwinky {self.__class__=} {self=} {other=}")
        return self._from_double_dict(_mul_schub_dicts(self._doubledict, self._parent(other)._doubledict))

    def __rmul__(self, other):
        # print("ASFJAdsajdSJ")
        # print(f"rmulwinky {self.__class__=} {self=} {other=}")
        return self._from_double_dict(_mul_schub_dicts(self._parent(other)._doubledict, self._doubledict))

    def __eq__(self, other):
        # elem1 = self.change_vars("y")  # count vars?
        # elem2 = self._parent(other).change_vars("y")
        # done = set()
        # for k, v in elem1._doubledict.items():
        #     done.add(k)
        #     if expand(v - elem2._doubledict.get(k, 0)) != 0:
        #         return False
        # for k, v in elem2._doubledict.items():
        #     if k in done:
        #         continue
        #     if expand(v - elem1._doubledict.get(k, 0)) != 0:
        #         return False
        return sympy.expand(self - other) == 0

    # def __str__(self):
    #     pieces = []
    #     keys = list(self._doubledict.keys())
    #     for k in sorted(keys, key=lambda b: (inv(b[0]), b[1], *b[0])):
    #         v = self._doubledict[k]
    #         dvar = "D"
    #         if sympy.expand(v) != 0:
    #             pieces += [
    #                 sympy.Mul(
    #                     v,
    #                     sympy.Symbol(
    #                         f"{dvar}S{self._parent._base_var}({list(k[0])}, {_varstr(k[1])})",
    #                         commutative=False,
    #                     )
    #                     if k[0] != (1, 2)
    #                     else 1,
    #                 ),
    #             ]
    #     return sympy.sstr(sympy.Add(*pieces, evaluate=False), order="none")

    # def __repr__(self):
    #     return str(self)

    # def _sympystr(self, *args):
    #     return str(self)
    def change_vars(self, cv):
        return DoubleSchubertAlgebraElement(_mul_schub_dicts({((1, 2), cv): 1}, self._doubledict), self)

    def _print_Add(self, *_):
        return "flagelnagel"

    def _sympyrepr(self, *_):
        return str(self)

    def as_coefficients_dict(self):
        # will not allow zeros
        # print("ASFJAdsajdSJ")

        return {k: v for k, v in self._doubledict.items() if sympy.expand(v) != 0}

    def normalize_coefficients(self, coeff_var):
        return self._parent([1, 2], coeff_var) * self

    def expand(self, *_a, **_):
        return sympy.expand(Add(*[yz.schubmult({(1, 2): v}, k[0], utils.poly_ring(self._parent._base_var), utils.poly_ring(k[1])).get((1, 2), 0) for k, v in self._doubledict.items()]))


# None is faster to store
class DoubleSchubertAlgebraElement_basis(Basic):
    coeff_varname = "y"

    def __init__(self, base_var="x"):
        # self._doubledict = _dict
        self._base_var = base_var
        # self._coeff_var = coeff_var if coeff_var else "y"

    def __call__(self, x, cv=None):
        if isinstance(x, list) or isinstance(x, tuple):
            if cv is None:
                cv = "y"
            elem = DoubleSchubertAlgebraElement({(tuple(permtrim(list(x))), cv): 1}, self)
        # elif isinstance(x, spr.SchubertPolynomial):
        #     if x._parent._base_var == self._base_var:
        #         elem_dict = {(x, utils.NoneVar): v for k, v in x._doubledict.items()}
        #         elem = DoubleSchubertAlgebraElement(elem_dict, self)
        #         if cv is not None:
        #             elem = self([1, 2], cv) * elem
        #     else:
        #         return self(x.expand(), cv)
        elif isinstance(x, DoubleSchubertAlgebraElement):
            if x._parent._base_var == self._base_var:
                elem = DoubleSchubertAlgebraElement(x._doubledict, self)
            else:
                return self(x.expand(), cv)
        else:
            if cv is None:
                cv = utils.NoneVar
                result = py.mult_poly({(1, 2): 1}, sympify(x), utils.poly_ring(self._base_var))
            else:
                result = yz.mult_poly({(1, 2): 1}, sympify(x), utils.poly_ring(self._base_var), utils.poly_ring(cv))
            elem = DoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()}, self)
        return elem

    # def _coerce_map_from(self, S):
    #     if isinstance(S, type(DSchub)):
    #         return True
    #     if isinstance(S, type(Schub)):
    #         return True
    #     if isinstance(S, Expr):
    #         return True
    #     return False


def _domul(*args):
    # print(f"domul {args=}")
    ret = 1
    for arg in args:
        if arg.is_Mul:
            for arg0 in arg.args:
                ret *= arg0
        else:
            ret *= arg
    return ret


def _doadd(*args):
    ret = 0
    # print(f"doadd {args=}")
    for arg in args:
        if arg.is_Add:
            for arg0 in arg.args:
                ret += arg0
        else:
            ret += arg
    # print(f"{ret=}")
    return ret


def get_postprocessor(cls):
    # print(f"{cls=} hey baby")
    if cls is Mul:
        return _domul
    if cls is Add:
        return _doadd
    return None


class SchubAdd(Add):
    def __new__(cls, *args, **_):
        if len(args) == 0:
            return 0
        res = args[0]
        for arg in args[1:]:
            res = res + arg
        return res


class SchubMul(Mul):
    def __new__(cls, *args, **_):
        if len(args) == 0:
            return 1
        res = args[0]
        for arg in args[1:]:
            res = res * arg
        return res


Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
    "Mul": [get_postprocessor(Mul)],
    "Add": [get_postprocessor(Add)],
}

# add.register_handlerclass((Expr, SchubAdd), SchubAdd)
# mul.register_handlerclass((Expr, SchubMul), SchubMul)


DSx = DoubleSchubertAlgebraElement_basis()
DoubleSchubertPolynomial = DoubleSchubertAlgebraElement

# def test(*args,**kwargs):
#     # print(f"test {args=} {kwargs=}")
# from sympy import Basic
# Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {'Mul': test}

# print(f"{Basic._constructor_postprocessor_mapping=}")

# is_Add = True
# is_Mul = True
# is_Add
# is_AlgebraicNumber
# is_Atom
# is_Boolean
# is_Derivative
# is_Dummy
# is_Equality
# is_Float
# is_Function
# is_Indexed
# is_Integer
# is_MatAdd
# is_MatMul
# is_Matrix
# is_Mul
# is_Not
# is_Number
# is_NumberSymbol
# is_Order
# is_Piecewise
# is_Point
# is_Poly
# is_Pow
# is_Rational
# is_Relational
# is_Symbol
# is_Vector
# is_Wild
# is_algebraic
# is_algebraic_expr
# is_antihermitian
# is_commutative
# is_comparable
# is_complex
# is_composite
# is_constant
# is_even
# is_extended_negative
# is_extended_nonnegative
# is_extended_nonpositive
# is_extended_nonzero
# is_extended_positive
# is_extended_real
# is_finite
# is_hermitian
# is_hypergeometric
# is_imaginary
# is_infinite
# is_integer
# is_irrational
# is_meromorphic
# is_negative
# is_noninteger
# is_nonnegative
# is_nonpositive
# is_nonzero
# is_number
# is_odd
# is_polar
# is_polynomial
# is_positive
# is_prime
# is_rational
# is_rational_function
# is_real
# is_scalar
# is_symbol
# is_transcendental
# is_zero
# is_polynomial = True
# is_Symbol = True
