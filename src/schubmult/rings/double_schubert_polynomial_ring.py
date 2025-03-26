from symengine import Add, Symbol, expand, sympify

import schubmult.schubmult_double as yz
from schubmult.perm_lib import add_perm_dict, permtrim

from ._utils import poly_ring
from .schubert_polynomial_ring import SchubertPolynomial


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

    def __rmul__(self, other):
        elem = other
        if not isinstance(other, DoubleDictAlgebraElement):
            elem = DoubleDictAlgebraElement_basis(other)

        ret = {}
        for k, v in elem._dict.items():
            ret = add_perm_dict(ret, yz.schubmult({k0: v0 * v for k0, v0 in self._dict.items()}, k, poly_ring(self._coeff_var), poly_ring(other._coeff_var)))
        return DoubleDictAlgebraElement(ret, self._gens)

    def __str__(self):
        pieces = []
        for k, v in self._dict.items():
            if expand(v) != 0:
                pieces += [v * Symbol(f"S({list(k)})")]
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
        if len(x) == 1:
            x = x[0]
        elif len(x) == 2 and (isinstance(x[0], list) or isinstance(x[0], tuple)):
            x, cv = x[0], x[1]

        if isinstance(x, list) or isinstance(x, tuple):
            # checking the input to avoid symmetrica crashing Sage, see trac 12924
            elem = DoubleDictAlgebraElement({(tuple(permtrim(list(x))), self._coeff_var): 1}, self._base_var, cv if cv else self._coeff_var)
        elif isinstance(x, DoubleDictAlgebraElement):
            if x._base_var == self._base_var:
                elem = DoubleDictAlgebraElement(x._dict, self._base_var, x._coeff_var)
            else:
                return self(x.expand())
        elif isinstance(x, SchubertPolynomial):
            if x.base_var == self._base_var:
                elem_dict = {}
                for k, v in x._dict.items():
                    res = yz.schubmult(
                        {(1, 2): v},
                        k,
                        poly_ring(self._coeff_var),
                        [0 for i in range(100)],
                    )
                    for k0, c0 in res.items():
                        elem_dict[(k0, self._coeff_var)] = (
                            elem_dict.get(
                                (k0, self._coeff_var),
                                0,
                            )
                            + c0
                        )
                elem = DoubleDictAlgebraElement(elem_dict, self._base_var, self._coeff_var)
            else:
                elem = self(x.expand())
        else:
            result = yz.mult_poly({(1, 2): 1}, sympify(x), poly_ring(self._base_var), poly_ring(self._coeff_var))
            elem = DoubleDictAlgebraElement({(k, self._coeff_var): v for k, v in result.items()}, self._base_var, self._coeff_var)
        return elem


DSchub = DoubleDictAlgebraElement_basis()
DoubleSchubertPolynomial = DoubleDictAlgebraElement
