"""`LazyPrintable`: drop-in for SymPy's ``Printable`` mixin that defers importing SymPy until an
object is actually printed.

SymPy's printers discover these classes through the ``_sympystr``/``_latex``/``_pretty`` hooks, so
subclassing ``Printable`` itself is unnecessary. The IPython ``_repr_*_`` methods forward to
``Printable``'s at call time, so ``init_printing``'s class-level patches still take effect.
"""


class LazyPrintable:
    __slots__ = ()

    def __str__(self):
        from sympy.printing.str import sstr

        return sstr(self, order=None)

    __repr__ = __str__

    def _repr_latex_(self):
        from sympy.printing.defaults import Printable

        return Printable._repr_latex_(self)

    def _repr_png_(self):
        from sympy.printing.defaults import Printable

        return Printable._repr_png_(self)

    def _repr_svg_(self):
        from sympy.printing.defaults import Printable

        return Printable._repr_svg_(self)
