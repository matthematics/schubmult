<a id="schubmult.utils._printable"></a>

# schubmult.utils.\_printable

`LazyPrintable`: drop-in for SymPy's ``Printable`` mixin that defers importing SymPy until an
object is actually printed.

SymPy's printers discover these classes through the ``_sympystr``/``_latex``/``_pretty`` hooks, so
subclassing ``Printable`` itself is unnecessary. The IPython ``_repr_*_`` methods forward to
``Printable``'s at call time, so ``init_printing``'s class-level patches still take effect.

