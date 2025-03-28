from functools import cache

import symengine
import sympy
from symengine import expand, sympify
from sympy import Add, Basic, Mul
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.printing.str import StrPrinter

import schubmult.rings._utils as utils
import schubmult.schubmult_double as yz
import schubmult.schubmult_py as py
from schubmult.perm_lib import add_perm_dict, inv, permtrim

# import utils.NoneVar, utils.ZeroVar, utils.poly_ring

_def_printer = StrPrinter({"order": "none"})
# numpy arrays
# sympy parsing
# quantum


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

    # print(f"{by_var2=}")

    for _vstr, _dict in by_var2.items():
        this_dict = {}
        for k, v in none_dict.items():
            this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in yz.schubmult(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
        results = add_perm_dict(results, this_dict)

    none_dict, none_dict2 = sorted([none_dict, none_dict2], key=lambda x: -len(x.keys()))
    for k, v in none_dict2.items():
        results = add_perm_dict(results, {(k1, utils.NoneVar): v1 * v for k1, v1 in py.schubmult(none_dict, k).items()})

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

    _base_var = "x"

    _op_priority = 1e200
    # __slots__ = ("_dict", "_parent")
    _kind = NumberKind
    is_Atom = True
    is_commutative = True
    # default_coeff_var = "y"

    def __new__(cls, _dict, *args, **kwargs):
        # print(f"{cls=} {_dict} {parent=} {args=} {kwargs=}")
        return DoubleSchubertAlgebraElement.__xnew_cached__(cls, sympy.Dict(_dict), *args, **kwargs)

    def __hash__(self):
        return hash(tuple(self.args))

    @property
    def _add_handler(self):
        return SchubAdd

    @property
    def _mul_handler(self):
        # print("profilating")
        return SchubMul

    @staticmethod
    def __xnew__(_class, _dict, *args, **kwargs):
        # print("Prong")
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
        obj = Expr.__new__(_class, _dict)
        # obj.make_args(pieces)
        # obj.args = pieces
        obj._doubledict = _dict
        obj._print_sum = sympy.Add(*[sympy.Mul(_dict[k], DSchubSymbol(DSx._base_var, k), evaluate=False) for k in sorted(_dict.keys(), key=lambda bob: (inv(bob[0]), str(bob[1]), *bob[0]))])
        return obj

    def _sympystr(self, printer):
        return printer._print_Add(self._print_sum)

    @staticmethod
    @cache
    def __xnew_cached__(_class, _dict, *args, **kwargs):
        return DoubleSchubertAlgebraElement.__xnew__(_class, _dict, *args, **kwargs)

    def _symengine_(self):
        return NotImplemented

    def _eval_simplify_(self):
        # print("Hey pretty baby")
        return self

    def _from_double_dict(self, _doubledict):
        return DoubleSchubertAlgebraElement(_doubledict, DSx)

    # def __add__(self, other):
    #     # print("ASFJASJ")
    #     # print(f"addwinky {self.__class__=} {self=} {other=}")
    #     # print(f"flarfknockle {self._doubledict=} {DSx(other)._doubledict=} {other=}")
    #     return SchubAdd(self, DSx(other))
    #     #return self._from_double_dict(ret)

    # def __radd__(self, other):
    #     # print("is wiz doing radd")
    #     #return self._from_double_dict(add_perm_dict(DSx(other)._doubledict, self._doubledict))
    #     return SchubAdd(DSx(other),self)

    # def __sub__(self, other):
    #     # print("ASFJAdsajdSJ")
    #     return SchubAdd(self,-DSx(other))
    #     #return self._from_double_dict(add_perm_dict(self._doubledict, {k: -v for k, v in DSx(other)._doubledict.items()}))

    # def __rsub__(self, other):
    #     # print("ASFJAdsajdSJ")
    #     # return self._from_double_dict(add_perm_dict(DSx(other)._doubledict, {k: -v for k, v in self._doubledict.items()}))
    #     return SchubAdd(DSx(other), -self)

    # def __neg__(self):
    #     if self.is_Add:
    #         return SchubAdd( *self.args)
    #     if self.is_Mul:
    #         return SchubMul( -self.args[0], *self.args[1:])
    #     return self._from_double_dict({k: -v for k, v in self._doubledict.items()})

    # def __mul__(self, other):
    #     # print("ASFJAdsajdSJ")
    #     # print(f"mulwinky {self.__class__=} {self=} {other=}")
    #     return SchubMul(self,DSx(other))

    # def __rmul__(self, other):
    #     # print("ASFJAdsajdSJ")
    #     # print(f"rmulwinky {self.__class__=} {self=} {other=}")
    #     #return self._from_double_dict(_mul_schub_dicts(DSx(other)._doubledict, self._doubledict))
    #     return SchubMul(DSx(other),self)

    def __eq__(self, other):
        # elem1 = self.change_vars("y")  # count vars?
        # elem2 = DSx(other).change_vars("y")
        if self.is_Add:
            if other.is_Add:
                if len(self.args) != len(other.args):
                    return False
                for i in range(len(self.args)):
                    if self.args[i] != other.args[i]:
                        return False
                return True
            return False
        if self.is_Mul:
            if other.is_Mul:
                if len(self.args) != len(other.args):
                    return False
                for i in range(len(self.args)):
                    if self.args[i] != other.args[i]:
                        return False
                return True
            return False
        elem1 = self.change_vars(0)
        elem2 = DSx(other).change_vars(0)

        done = set()
        for k, v in elem1._doubledict.items():
            done.add(k)
            if expand(v - elem2._doubledict.get(k, 0)) != 0:
                return False
        for k, v in elem2._doubledict.items():
            if k in done:
                continue
            if expand(v - elem1._doubledict.get(k, 0)) != 0:
                return False
        return True

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
    #                         f"{dvar}S{DSx._base_var}({list(k[0])}, {_varstr(k[1])})",
    #                         commutative=False,
    #                     )
    #                     if k[0] != (1, 2)
    #                     else 1,
    #                 ),
    #             ]
    #     return sympy.sstr(sympy.Add(*pieces, evaluate=False), order="none")

    # def __repr__(self):
    #     return str(self)
    # def _as_ordered_terms(self, *args, **kwargs):

    # def _sympystr(self, *args):
    #     return str(self)
    def _eval_simplify(self, *args, **kwargs):  # noqa: ARG002
        return self._from_double_dict({k: sympify(sympy.simplify(v)) for k, v in self._doubledict.items()})

    def change_vars(self, cv):
        return self._from_double_dict(_mul_schub_dicts({((1, 2), cv): 1}, self._doubledict))

    def _print_Add(self, *_):
        return "flagelnagel"

    def _sympyrepr(self, *_):
        return str(self)

    def as_coefficients_dict(self):
        # will not allow zeros
        # print("ASFJAdsajdSJ")

        return {k: v for k, v in self._doubledict.items() if sympy.expand(v) != 0}

    def normalize_coefficients(self, coeff_var):
        return DSx([1, 2], coeff_var) * self

    def expand(self, *_a, **_):
        if isinstance(self, SchubAdd):
            return self.doit().expand()
        if isinstance(self, SchubMul):
            return self.doit().expand()
        return sympy.sympify(expand(symengine.Add(*[yz.schubmult({(1, 2): v}, k[0], utils.poly_ring(DSx._base_var), utils.poly_ring(k[1])).get((1, 2), 0) for k, v in self._doubledict.items()])))


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
            elem = DoubleSchubertAlgebraElement({(tuple(permtrim(list(x))), cv): 1})
        # elif isinstance(x, spr.SchubertPolynomial):
        #     if x._parent._base_var == self._base_var:
        #         elem_dict = {(x, utils.NoneVar): v for k, v in x._doubledict.items()}
        #         elem = DoubleSchubertAlgebraElement(elem_dict, self)
        #         if cv is not None:
        #             elem = self([1, 2], cv) * elem
        #     else:
        #         return self(x.expand(), cv)
        elif isinstance(x, DoubleSchubertAlgebraElement):
            if x.is_Add or x.is_Mul:
                return x
            if x._base_var == self._base_var:
                elem = DoubleSchubertAlgebraElement(x._doubledict)  # , self)
            else:
                return self(x.expand(), cv)
        else:
            if cv is None:
                cv = utils.NoneVar
                result = py.mult_poly({(1, 2): 1}, sympify(x), utils.poly_ring(self._base_var))
            else:
                result = yz.mult_poly({(1, 2): 1}, sympify(x), utils.poly_ring(self._base_var), utils.poly_ring(cv))
            elem = DoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()})
        return elem

    # def _coerce_map_from(self, S):
    #     if isinstance(S, type(DSchub)):
    #         return True
    #     if isinstance(S, type(Schub)):
    #         return True
    #     if isinstance(S, Expr):
    #         return True
    #     return False


def _do_schub_mul(a, b):
    A = DSx(a)
    B = DSx(b)
    return A._from_double_dict(_mul_schub_dicts(A._doubledict, B._doubledict))


def _do_schub_add(a, b):
    A = DSx(a)
    B = DSx(b)
    return A._from_double_dict(add_perm_dict(A._doubledict, B._doubledict))


def _domul(*args):
    # print(f"domul {args=}")
    return SchubMul(*args)


def _doadd(*args):
    return SchubAdd(*args)


def get_postprocessor(cls):
    # print(f"{cls=} hey baby")
    if cls is Mul:
        return lambda expr: SchubMul(*expr.args, evaluate=True).doit()
    if cls is Add:
        return lambda expr: SchubAdd(*expr.args, evaluate=True).doit()
    return None

    # def _sympystr(self, printer):
    #     return _def_printer._print(f"SchubMul({self.args}")

    # @classmethod
    # def doit(cls, expr):
    #     # print("Mul pip monger")
    # def __new__(cls, *args, **_):
    #     if len(args) == 0:
    #         return 1
    #     res = args[0]
    #     for arg in args[1:]:
    #         res = res.__mul__(arg)
    #     return res
    # def _sympystr(self, printer):
    #     printer._print_Mul(*self.args, order='none')

    # def __str__(self):
    #     return self._sympystr(printer=_def_printer)

    # def __repr__(self):
    #     return str(self)


# Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
#     "Mul": [get_postprocessor(Mul)],
#     "Add": [get_postprocessor(Add)],
# }

# add.register_handlerclass((Expr, SchubAdd), SchubAdd)
# mul.register_handlerclass((Expr, SchubMul), SchubMul)


DSx = DoubleSchubertAlgebraElement_basis()


def Sx(x):
    return DSx(x, utils.NoneVar)


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


class SchubAdd(DoubleSchubertAlgebraElement, Add):
    is_Add = True

    def __new__(cls, *args, evaluate=True, _sympify=True):
        obj = Add.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)
        if evaluate:
            return obj.doit()
        return obj

    def doit(self):
        ret = self.args[0]
        for arg in self.args[1:]:
            if arg.is_Add or arg.is_Mul:
                arg = arg.doit()
            ret = _do_schub_add(ret, arg)
        return ret

    # def _sympystr(self, printer):
    #     return _def_printer._print(f"SchubAdd({self.args}")


class SchubMul(DoubleSchubertAlgebraElement, Mul):
    is_Mul = True

    def __new__(cls, *args, evaluate=True, _sympify=True):
        args, a, b = Mul.flatten(list(args))
        obj = Mul.__new__(cls, *args, evaluate=evaluate, _sympify=True)

        if evaluate:
            return obj.doit()
        return obj

    def doit(self):
        ret = self.args[0]
        for arg in self.args[1:]:
            if arg.is_Add or arg.is_Mul:
                arg = arg.doit()
            ret = _do_schub_mul(ret, arg)
        return ret


Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
    "Mul": [get_postprocessor(Mul)],
    "Add": [get_postprocessor(Add)],
}

# add.register_handlerclass((Add, SchubAdd), SchubAdd)
# mul.register_handlerclass((Mul, SchubMul), SchubMul)
