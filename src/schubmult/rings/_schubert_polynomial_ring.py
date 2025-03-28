import sympy
from symengine import sympify
from sympy import Add, Basic, Mul

import schubmult.rings._double_schubert_polynomial_ring as dsb
import schubmult.rings._utils as utils
import schubmult.schubmult_py as py
from schubmult.perm_lib import permtrim


class SchubertAlgebraElement(dsb.DoubleSchubertAlgebraElement):
    """Algebra with sympy coefficients
    and a dict basis

    """

    def __hash__(self):
        return hash(self._singledict)

    def __new__(cls, _dict, parent, *args, **kwargs):
        # print(f"{cls=} {_dict=} {args=} {kwargs=}")
        # print(f"pifflanfadonky {_dict=}")
        if kwargs.get("double_dict", False):
            # print(f"{_dict=}")
            _doubledict = _dict
            _dict = {k[0]: v for k, v in _doubledict.items()}
        else:
            _doubledict = {(k, utils.NoneVar): v for k, v in _dict.items()}
            # print(f"frafunka {_dict=} {kwargs=}")
            # print(f"{_doubledict=}")
        obj = dsb.DoubleSchubertAlgebraElement.__new__(cls, _doubledict, parent, *args, **kwargs)
        obj._singledict = sympy.Dict(_dict)
        # print(f"{obj._doubledict=} {obj._singledict=}")
        # print(f"{obj._doubledict=} {obj._parent=}")
        # obj._doubledict = sympy.Dict(_doubledict)
        return obj

    def __eq__(self, other):
        elem1 = self.change_vars(0)  # count vars?
        elem2 = self._parent(other).change_vars(0)
        done = set()
        for k, v in elem1._doubledict.items():
            done.add(k)
            if sympy.expand(v - elem2._doubledict.get(k, 0)) != 0:
                return False
        for k, v in elem2._doubledict.items():
            if k in done:
                continue
            if sympy.expand(v - elem1._doubledict.get(k, 0)) != 0:
                return False
        return True


class SchubertAlgebraElement_basis(dsb.DoubleSchubertAlgebraElement_basis):
    def __init__(self, base_var="x"):
        self._base_var = base_var
        super().__init__(base_var)

    def __call__(self, x):
        # print(f"prevnool {x=} {x.__class__=} {isinstance(x, SchubertAlgebraElement)=}")
        if isinstance(x, list) or isinstance(x, tuple):
            elem = SchubertAlgebraElement({tuple(permtrim(list(x))): 1}, self)
        elif isinstance(x, SchubertAlgebraElement):
            if x._parent._base_var == self._base_var:
                # print(f"piffdonk {x=} {x._singledict=}")
                elem = SchubertAlgebraElement(x._singledict, self)
            else:
                return self(x.expand())
        elif isinstance(x, dsb.DoubleSchubertAlgebraElement):
            elem = dsb.DoubleSchubertAlgebraElement(x._doubledict, self)
        else:
            try:
                result = py.mult_poly({(1, 2): 1}, sympify(x), utils.poly_ring(self._base_var))  # this will raise an error if no good
                elem = SchubertAlgebraElement(result, self)
            except Exception:
                return NotImplemented
        return elem


Sx = SchubertAlgebraElement_basis()
SchubertPolynomial = SchubertAlgebraElement

Basic._constructor_postprocessor_mapping[SchubertAlgebraElement] = {"Add": [dsb.get_postprocessor(Add)], "Mul": [dsb.get_postprocessor(Mul)]}
