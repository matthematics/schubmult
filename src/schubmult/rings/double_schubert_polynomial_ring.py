from symengine import Add, Symbol, expand, sympify, Basic

import schubmult.schubmult_double as yz
from schubmult.perm_lib import add_perm_dict, permtrim

from ._utils import poly_ring
from .schubert_polynomial_ring import Schub, SchubertPolynomial


class DoubleDictAlgebraElement:
    """Algebra with sympy coefficients
    and a dict basis
    """

    def __init__(self, _dict, base_var="x", coeff_var=None):
        self._dict = _dict
        self._base_var = base_var
        self._coeff_var = coeff_var if coeff_var else "y"

    def __add__(self, other):
        if isinstance(other, DoubleDictAlgebraElement):
            return DoubleDictAlgebraElement(add_perm_dict(self._dict, other._dict), self._base_var, self._coeff_var)
        else:
            return DoubleDictAlgebraElement(add_perm_dict(self._dict, self(other)._dict), self._base_var, self._coeff_var)

    def __sub__(self, other):
        if isinstance(other, DoubleDictAlgebraElement):
            return self.__add__(DoubleDictAlgebraElement({k: -v for k, v in other._dict.items()}, self._base_var, self._coeff_var))
        else:
            return self.__add__(DoubleDictAlgebraElement({k: -v for k, v in self(other)._dict.items()}, self._base_var, self._coeff_var))

    def __neg__(self):
        return DoubleDictAlgebraElement({k: -v for k, v in self._dict.items()}, self._base_var, self._coeff_var)

    def _sympify_(self):
        raise TypeError("You cannot sympify this")

    def __mul__(self, other):
        elem = other
        print(f"seref perfnool {elem=}")
        if isinstance(other, Basic) or self._parent._coerce_map_from(other._parent): #not isinstance(other, DictAlgebraElement):
            elem = self._parent(other)
        elif other._parent._coerce_map_from(self._parent):
            print("boinle")
            return other.__mul__(self)
        else:
            raise TypeError("Noip {type(self)}")
        print(f"{self=} {elem=}")
        ret = {}
        for k, v in elem._dict.items():
            for k0, v0 in self._dict.items():
                ret = {(k1, k0[1]): v1 for k1, v1 in add_perm_dict(ret, yz.schubmult({k0[0]: v0 * v}, k[0], poly_ring(k0[1]), poly_ring(k[1]))).items()}
        return DoubleDictAlgebraElement(ret, self._base_var, self._coeff_var)

    def __str__(self):
        pieces = []
        for k, v in self._dict.items():
            if expand(v) != 0:
                pieces += [v * Symbol(f"DS{self._base_var}({list(k[0])}, {k[1]})")]
        return str(Add(*pieces))

    def __repr__(self):
        return self.__str__()


class DoubleDictAlgebraElement_basis:
    def __init__(self, base_var="x", coeff_var=None):
        # self._dict = _dict
        self._base_var = base_var
        self._coeff_var = coeff_var if coeff_var else "y"

    def __call__(self, *x):
        cv = None
        print(f"bagel {x=}")
        if len(x) == 1:
            x = x[0]
        elif len(x) == 2 and (isinstance(x[0], list) or isinstance(x[0], tuple)):
            x, cv = x[0], x[1]
        else:
            print(f"bongdunket {x=}")
        print(f"{x=} {cv=}")
        if isinstance(x, list) or isinstance(x, tuple):
            print(f"{x=} {cv=}")
            # checking the input to avoid symmetrica crashing Sage, see trac 12924
            cv = self._coeff_var if cv is None else cv
            elem = DoubleDictAlgebraElement({(tuple(permtrim(list(x))), cv): 1}, self._base_var, cv)
        elif isinstance(x, DoubleDictAlgebraElement):
            if x._base_var == self._base_var:
                elem = DoubleDictAlgebraElement(x._dict, self._base_var, x._coeff_var)
            else:
                return self(x.expand())
        elif isinstance(x, SchubertPolynomial):
            if x._base_var == self._base_var:
                print(f"{x=} knoff boogle")
                elem_dict = {}
                for k, v in x._dict.items():
                    res = yz.schubmult(
                        {(1, 2): v},
                        k,
                        poly_ring(self._coeff_var),
                        [
                            0
                            for i in range(100)
                        ],
                    )
                print(f"proif bunkgle {res=}")
                for k0, c0 in res.items():
                    elem_dict[(k0, self._coeff_var)] = elem_dict.get(
                        (k0, self._coeff_var),
                        0)+ c0
                print(f"blor {elem_dict=}")
                elem = DoubleDictAlgebraElement(elem_dict,self._base_var,self._coeff_var)
            else:
                return self(x.expand())
        else:
            result = yz.mult_poly({(1, 2): 1}, sympify(x), poly_ring(self._base_var), poly_ring(self._coeff_var))
            elem = DoubleDictAlgebraElement({(k, self._coeff_var): v for k, v in result.items()}, self._base_var, self._coeff_var)
        elem._parent = self
        print(f"{elem=} {x=}")
        return elem

    def _coerce_map_from(self, S):
        if isinstance(S, type(DSchub)):
            return True
        if isinstance(S, type(Schub)):
            return True
        if isinstance(S, Basic):
            return True
        return False




DSchub = DoubleDictAlgebraElement_basis()
DoubleSchubertPolynomial = DoubleDictAlgebraElement
