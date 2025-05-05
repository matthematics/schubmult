import sys

import symengine.lib.symengine_wrapper as sw
from sympy import sstr


class SymPolyWrap(sw.PyFunction, sw.Symbol):

    def __new__(cls, *args):
        return sw.PyFunction.__new__(cls, *args)

    def __init__(self, *args):
        super(sw.Symbol, self).__init__(sstr(self.pyobject()))

    def __repr__(self):
        return sstr(self.pyobject())

    def __str__(self):
        return sstr(self.pyobject())

    @property
    def genvars(self):
        return self.pyobject().genvars

    @property
    def coeffvars(self):
        return self.pyobject().coeffvars
    
    def _symengine_(self):
        return self
    
    def _sympy_(self):
        return self.pyobject()

    @property
    def func(self):
        return self.pyobject().func