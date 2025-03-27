import sympy
from symengine import expand, sympify

import schubmult.schubmult_double as yz
from schubmult.perm_lib import add_perm_dict, permtrim

from ._utils import poly_ring
from .schubert_polynomial_ring import SchubertPolynomial

# numpy arrays
# sympy parsing
# quantum


def _mul_schub_dicts(dict1, dict2):
    by_var = {}

    for k, v in dict1.items():
        if k[1] not in by_var:
            by_var[k[1]] = {}
        by_var[k[1]][k[0]] = v

    results = {}

    for _vstr, _dict in by_var.items():
        this_dict = {}
        for k, v in dict2.items():
            this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in yz.schubmult(_dict, k[0], poly_ring(_vstr), poly_ring(k[1])).items()})
        results.update(this_dict)
    return results


class DoubleDictAlgebraElement:
    """Algebra with sympy coefficients
    and a dict basis
    """

    def __init__(self, _dict, parent):
        self._dict = _dict
        self._parent = parent

    def __add__(self, other):
        return DoubleDictAlgebraElement(add_perm_dict(self._dict, self._parent(other)._dict), self._parent)

    def __radd__(self, other):
        return DoubleDictAlgebraElement(add_perm_dict(self._parent(other)._dict, self._dict), self._parent)

    def __sub__(self, other):
        return DoubleDictAlgebraElement(add_perm_dict(self._dict, {k: -v for k, v in self._parent(other)._dict.items()}), self._parent)

    def __rsub__(self, other):
        return DoubleDictAlgebraElement(add_perm_dict(self._parent(other)._dict, {k: -v for k, v in self._dict.items()}), self._parent)

    def __neg__(self):
        return DoubleDictAlgebraElement({k: -v for k, v in self._dict.items()}, self._parent)

    def __mul__(self, other):
        return DoubleDictAlgebraElement(_mul_schub_dicts(self._dict, self._parent(other)._dict), self._parent)

    def __rmul__(self, other):
        return DoubleDictAlgebraElement(_mul_schub_dicts(self._parent(other)._dict, self._dict), self._parent)

    def __eq__(self, other):
        one = self._parent([1, 2])
        elem1 = one * self
        elem2 = one * self._parent(other)
        print(f"{elem1=}")
        print(f"{elem2=}")
        done = set()
        for k, v in elem1._dict.items():
            done.add(k)
            if expand(v - elem2._dict.get(k, 0)) != 0:
                return False
        for k, v in elem2._dict.items():
            if k in done:
                continue
            if expand(v - elem1._dict.get(k, 0)) != 0:
                return False
        return True

    def __str__(self):
        pieces = []
        for k, v in self._dict.items():
            if expand(v) != 0:
                pieces += [sympy.Mul(v, sympy.Symbol(f"DS{self._parent._base_var}({list(k[0])}, '{k[1]}')", commutative=False))]
        return sympy.sstr(sympy.Add(*pieces), order=lambda order: pieces)  # use sstr

    def __repr__(self):
        return self.__str__()

    def expand(self):
        ret = 0
        keys = list(self._dicts.keys())
        for k in sorted(keys):
            v = self._dict[k]
            ret += yz.schubmult({(1, 2): v}, k[0], poly_ring(self._parent._base_var), poly_ring(k[1])).get((1, 2), 0)
        return ret


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
        if cv is None:
            cv = self._coeff_var
        if isinstance(x, list) or isinstance(x, tuple):
            elem = DoubleDictAlgebraElement({(tuple(permtrim(list(x))), cv): 1}, self)
        elif isinstance(x, DoubleDictAlgebraElement):
            if x._parent._base_var == self._base_var:
                elem = DoubleDictAlgebraElement(x._dict, self)
            else:
                return self(x.expand())
        elif isinstance(x, SchubertPolynomial):
            if x._parent._base_var == self._base_var:
                elem_dict = {}
                for k, v in x._dict.items():
                    res = yz.schubmult(
                        {(1, 2): v},
                        k,
                        poly_ring(self._coeff_var),
                        [0 for i in range(100)],
                    )
                for k0, c0 in res.items():
                    elem_dict[(k0, self._coeff_var)] = elem_dict.get((k0, self._coeff_var), 0) + c0
                elem = DoubleDictAlgebraElement(elem_dict, self)
            else:
                return self(x.expand())
        else:
            result = yz.mult_poly({(1, 2): 1}, sympify(x), poly_ring(self._base_var), poly_ring(self._coeff_var))
            elem = DoubleDictAlgebraElement({(k, self._coeff_var): v for k, v in result.items()}, self)
        elem._parent = self
        return elem

    # def _coerce_map_from(self, S):
    #     if isinstance(S, type(DSchub)):
    #         return True
    #     if isinstance(S, type(Schub)):
    #         return True
    #     if isinstance(S, Basic):
    #         return True
    #     return False


DSchub = DoubleDictAlgebraElement_basis()
DoubleSchubertPolynomial = DoubleDictAlgebraElement
