"""Display atoms for ring basis elements.

Every ring's ``printing_term(key)`` returns a `PrintingTerm` subclass instance: an inert atom
that knows how to render itself for ``str``, pretty printing, and LaTeX. Instances are interned via
cached ``__xnew_cached__`` constructors so equal keys give identical objects. The subclasses cover
single/double Schubert (``S``/``DS``), quantum (``QS``/``QDS``, ``QPS``/``QPDS``), Grothendieck
(``G``/``DG``), separated descents (``Xi``), and a `GenericPrintingTerm` ``name(key)`` fallback.

A `PrintingTerm` is a plain Python object, so defining and creating them needs no SymPy. SymPy sees
one through ``_sympy_``, which returns a copy whose class also derives from SymPy's ``Expr`` (an
``Expr`` atom with ``args == ()``, so SymPy never traverses into it).
"""

import types
from functools import cache

import schubmult.symbolic as ssymb
from schubmult.combinatorics.permutation import Permutation
from schubmult.utils._printable import LazyPrintable


class _Hidden:
    """Descriptor that makes an inherited attribute look absent (``hasattr`` is False)."""

    def __get__(self, obj, objtype=None):
        raise AttributeError


@cache
def _sympy_class(cls):
    """``cls`` combined with SymPy's ``Expr``; ``cls``'s methods take precedence except arithmetic."""
    (expr,) = types.resolve_bases((ssymb.Expr,))

    def body(ns):
        ns["__module__"] = cls.__module__
        ns["__qualname__"] = cls.__qualname__
        # SymEngine recurses on _sympy_ before checking for sympy.Basic, so it must not be visible here
        ns["_sympy_"] = _Hidden()
        for name in _ARITHMETIC:
            ns[name] = getattr(expr, name)

    return types.new_class(cls.__name__, (cls, expr), exec_body=body)


_ARITHMETIC = ("__add__", "__radd__", "__sub__", "__rsub__", "__mul__", "__rmul__", "__truediv__", "__rtruediv__", "__pow__", "__neg__")


# Atomic Schubert polynomial
class PrintingTerm(LazyPrintable):
    """Base display atom carrying a key, generating set, coefficient generating set, and prefix."""

    is_Atom = True
    is_number = False

    def __new__(cls, k, genset, coeff_genset, prefix=""):
        obj = object.__new__(cls)
        obj._prefix = prefix
        obj._key = k
        obj._genset = genset
        obj._coeff_genset = coeff_genset
        obj._perm = k
        return obj

    def _sympy_(self):
        obj = self.__dict__.get("_sympy_obj")
        if obj is None:
            sympy_cls = _sympy_class(type(self))
            # skip the plain __new__ chain: allocate through Expr (Basic.__new__), the next base after LazyPrintable
            obj = super(LazyPrintable, sympy_cls).__new__(sympy_cls)
            obj.__dict__.update(self.__dict__)
            self._sympy_obj = obj
        return obj

    def __add__(self, other):
        return self._sympy_() + other

    def __radd__(self, other):
        return other + self._sympy_()

    def __sub__(self, other):
        return self._sympy_() - other

    def __rsub__(self, other):
        return other - self._sympy_()

    def __mul__(self, other):
        return self._sympy_() * other

    def __rmul__(self, other):
        return other * self._sympy_()

    def __truediv__(self, other):
        return self._sympy_() / other

    def __rtruediv__(self, other):
        return other / self._sympy_()

    def __pow__(self, other):
        return self._sympy_() ** other

    def __neg__(self):
        return -self._sympy_()

    def __hash__(self):
        return hash((self._key, self._genset, self._coeff_genset, "AbS", self._prefix))

    # def __init__(self, k, genset, coeff_genset):
    #     super().__init__(k, genset, coeff_genset)

    # def __reduce__(self):
    #     return (self.__class__, self.args)

    # def _eval_expand_basic(self, *args, **kwargs):
    #     return self._ring_elem.as_polynomial()

    # def __reduce__(self):
    #     print("Fleff")
    #     return (self.__class__, self.args)

    # def expand(self, deep=True, *args, **kwargs):
    #     if not deep:
    #         return self._from_dict({k: expand(v) for k, v in self._ring_elem.items()})
    #     return sympify(expand(sympify(self.as_polynomial())))

    # def as_polynomial(self):
    #     return self._cached_schubpoly(self._key)

    @property
    def genset(self):
        return self._genset

    @property
    def coeff_genset(self):
        return self._coeff_genset

    @property
    def args(self):
        # Return empty tuple since all subclasses are Atoms - SymPy shouldn't traverse into them
        return ()


class GenericPrintingTerm(PrintingTerm):
    """Displays a key as ``name(key)`` (e.g. ``AGx(perm, n)``, ``N(2, 1)``); the identity key prints as ``1``."""

    is_Atom = True

    _pretty_schub_char = "𝔖"

    def __hash__(self):
        return hash((self._key, self._genset, self._coeff_genset, "dasiub", self._name))

    def __new__(cls, k, name):
        return GenericPrintingTerm.__xnew_cached__(cls, k, name)

    @staticmethod
    def __xnew__(_class, k, name):
        obj = PrintingTerm.__new__(_class, k, None, None, None)
        obj._name = name
        return obj

    def _sympystr(self, printer):
        key = self._key
        if self._key == Permutation([]):
            return printer.doprint(ssymb.S.One)
        return printer.doprint(f"{self._name}{printer.doprint(key)}")

    def __reduce__(self):
        return (self.__class__, self.args)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, name):
        return GenericPrintingTerm.__xnew__(_class, k, name)


class TypedPrintingTerm(PrintingTerm):
    """Displays a key by delegating to the key's own printer (used for keys that are themselves
    printable objects such as RC graphs).
    """

    is_Atom = True

    _pretty_schub_char = "𝔖"

    def __hash__(self):
        return hash((self._key, "dasiub"))

    def __new__(cls, k):
        return TypedPrintingTerm.__xnew_cached__(cls, k)

    @property
    def args(self):
        # Return the key wrapped in a tuple to prevent sympy from trying to traverse it
        return (self._key,)

    @staticmethod
    def __xnew__(_class, k):
        return PrintingTerm.__new__(_class, k, None, None, None)

    def _sympystr(self, printer):
        key = self._key
        return printer._print(key)

    def _pretty(self, printer):
        key = self._key
        return printer._print(key)

    def _latex(self, printer):
        key = self._key
        return printer._print(key)

    def __reduce__(self):
        return (self.__class__, self.args)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k):
        return TypedPrintingTerm.__xnew__(_class, k)


class DSchubPoly(PrintingTerm):
    """Schubert polynomial term: ``S<genset>(perm)`` or ``DS<genset>(perm, <coeff_genset>)``."""

    is_Atom = True

    _pretty_schub_char = "𝔖"

    def __hash__(self):
        return hash((self._key, self._genset, self._coeff_genset, "dasiub"))

    def __new__(cls, k, genset, coeff_genset, prefix=""):
        return DSchubPoly.__xnew_cached__(cls, k, genset, coeff_genset, prefix)

    @staticmethod
    def __xnew__(_class, k, genset, coeff_genset, prefix):
        return PrintingTerm.__new__(_class, k, genset, coeff_genset, prefix)

    def _sympystr(self, printer):
        key = self._key
        if self._key == Permutation([]):
            return printer.doprint(ssymb.S.One)
        if self._coeff_genset is None:
            return printer.doprint(f"{self._prefix}S{self._genset}({printer.doprint(key)})")
        return printer.doprint(f"{self._prefix}DS{self._genset}({printer.doprint(key)}, {self._coeff_genset})")

    def _pretty(self, printer):
        key = self._key
        gl = self._genset
        if key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(key)
        if self._coeff_genset is None:
            return printer._print_Function(ssymb.Function(f"{self._prefix}{self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(gl)))
        return printer._print_Function(ssymb.Function(f"{self._prefix}{self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(f"{self._genset}; {self._coeff_genset}")))

    def _latex(self, printer):
        key = self._key
        gl = self._genset
        if key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(key)
        if self._coeff_genset is None:
            return printer._print_Function(ssymb.Function(f"{self._prefix}" + "\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(ssymb.Symbol(gl)))
        return printer._print_Function(ssymb.Function(f"{self._prefix}" + "\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(ssymb.Symbol(f"{{{self._genset}}}; {{{self._coeff_genset}}}")))

    def __reduce__(self):
        return (self.__class__, self.args)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset, coeff_genset, prefix):
        return DSchubPoly.__xnew__(_class, k, genset, coeff_genset, prefix)


class SepDescSchubPoly(PrintingTerm):
    """Separated-descents term for the key ``(perm, numvars)``: ``Xi_{perm}^{numvars}``."""

    is_Atom = True

    _pretty_schub_char = "Ξ"

    def __hash__(self):
        return hash((self._key, self._genset, self._coeff_genset, "SepDesc"))

    def __new__(cls, k, genset, coeff_genset):
        return SepDescSchubPoly.__xnew_cached__(cls, k, genset, coeff_genset)

    @staticmethod
    def __xnew__(_class, k, genset, coeff_genset):
        # k is a tuple (perm, length)
        return PrintingTerm.__new__(_class, k, genset, coeff_genset, "")

    @property
    def args(self):
        # Return empty tuple since we're an Atom - SymPy shouldn't traverse into us
        return ()

    def _sympystr(self, printer):
        perm, length = self._key
        if perm == Permutation([]):
            return printer.doprint(ssymb.S.One)
        return printer.doprint(f"Xi_{{{printer.doprint(perm)}}}^{{{length}}}")

    def _pretty(self, printer):
        perm, length = self._key
        if perm == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(perm)
        superscript = printer._print(length)
        return printer._print(ssymb.Symbol(f"{self._pretty_schub_char}_{subscript}^{superscript}"))

    def _latex(self, printer):
        perm, length = self._key
        if perm == Permutation([]):
            return printer._print(ssymb.S.One)
        # Format permutation as tuple for LaTeX
        perm_str = ''.join(str(x) for x in perm)
        return f"\\Xi_{{{perm_str}}}^{{{length}}}"

    def __reduce__(self):
        return (self.__class__, self.args)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset, coeff_genset):
        return SepDescSchubPoly.__xnew__(_class, k, genset, coeff_genset)


class QDSchubPoly(PrintingTerm):
    """Quantum Schubert term: ``QS<genset>(perm)`` or ``QDS<genset>(perm, <coeff_genset>)``."""

    is_Atom = True

    def __hash__(self):
        return hash((self._key, self._genset, self._coeff_genset, "AbSsf", "q"))

    _pretty_schub_char = "𝕼𝔖"

    def __new__(cls, k, genset, coeff_genset):
        return QDSchubPoly.__xnew_cached__(cls, k, genset, coeff_genset)

    @staticmethod
    def __xnew__(_class, k, genset, coeff_genset):
        return PrintingTerm.__new__(_class, k, genset, coeff_genset)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset, coeff_genset):
        return QDSchubPoly.__xnew__(_class, k, genset, coeff_genset)

    def _sympystr(self, printer):
        if self.coeff_genset is None:
            return printer.doprint(f"QS{self.genset}({printer.doprint(self._perm)})")
        return printer.doprint(f"QDS{self.genset}({printer.doprint(self._perm)}, {self.coeff_genset})")

    def _pretty(self, printer):
        if self._key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer.doprint(int("".join([str(i) for i in self._key])))
        if self.coeff_genset is None:
            return printer._print_Function(ssymb.Function(f"{self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(self.genset)))
        return printer._print_Function(ssymb.Function(f"{self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(f"{self.genset}; {self.coeff_genset}")))

    def _latex(self, printer):
        if self._key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = ssymb.sstr(self._key)
        if self.coeff_genset is None:
            return printer._print_Function(ssymb.Function("\\widetilde{\\mathfrak{S}}" + f"_{'{' + subscript + '}'}")(ssymb.Symbol(self.genset)))
        return printer._print_Function(ssymb.Function("\\widetilde{\\mathfrak{S}}" + f"_{'{' + subscript + '}'}")(ssymb.Symbol(f"{self.genset}; {self.coeff_genset}")))

    def __reduce__(self):
        return (self.__class__, self.args)


class PQDSchubPoly(PrintingTerm):
    """Parabolic quantum Schubert term, tagged with the index composition:
    ``QPS<genset>(comp)(perm)`` or ``QPDS<genset>(comp)(perm, <coeff_genset>)``.
    """

    is_Atom = True

    def __hash__(self):
        return hash((self._key, self._genset, self._coeff_genset, self.index_comp, "asfafAbS"))

    _pretty_schub_char = "𝕼𝔖"

    def __new__(cls, k, genset, coeff_genset, index_comp):
        return PQDSchubPoly.__xnew_cached__(cls, k, genset, coeff_genset, index_comp)

    @staticmethod
    def __xnew__(_class, k, genset, coeff_genset, index_comp):
        obj = PrintingTerm.__new__(_class, k, genset, coeff_genset)
        obj._perm = k
        obj._key = k
        obj._index_comp = index_comp
        return obj

    @property
    def index_comp(self):
        return self._index_comp

    @property
    def args(self):
        return (tuple(self._key), self._basis, tuple(self._index_comp))

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset, coeff_genset, index_comp):
        return PQDSchubPoly.__xnew__(_class, k, genset, coeff_genset, index_comp)

    def _sympystr(self, printer):
        if self.coeff_genset is None:
            return printer.doprint(f"QPS{self.genset}{(tuple(self.index_comp))}({printer.doprint(self._perm)})")
        return printer.doprint(f"QPDS{self.genset}{tuple(self.index_comp)}({printer.doprint(self._perm)}, {self.coeff_genset})")

    def _pretty(self, printer):
        if self._key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(int("".join([str(i) for i in self._key])))
        if self.coeff_genset is None:
            return printer._print_Function(ssymb.Function(f" {self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(f"{self.genset} | {self.index_comp}")))
        return printer._print_Function(ssymb.Function(f" {self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(f"{self.genset}; {self.coeff_genset} | {self.index_comp}")))

    def _latex(self, printer):
        if self._key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(self._key)
        supscript = printer._print(tuple(self.index_comp))
        if self.coeff_genset is None:
            return printer._print_Function(ssymb.Function("\\widetilde{\\mathfrak{S}}" + f"^{'{'}{supscript}{'}'}_{'{' + subscript + '}'}")(ssymb.Symbol(self.genset)))
        return printer._print_Function(
            ssymb.Function("\\widetilde{\\mathfrak{S}}" + f"^{'{'}{supscript}{'}'}_{'{' + subscript + '}'}")(ssymb.Symbol(f"{self.genset}; {self.coeff_genset}")),
        )


class GrothendieckPoly(PrintingTerm):
    """Grothendieck term ``G<genset>(perm)``, or ``G<genset>(perm, numvars)`` when the key carries a variable count."""

    is_Atom = True

    _pretty_schub_char = "𝔊"

    def __hash__(self):
        return hash((self._key, self._genset, "GrothendieckPoly"))

    def __new__(cls, k, genset, prefix=""):
        return GrothendieckPoly.__xnew_cached__(cls, k, genset, prefix)

    @staticmethod
    def __xnew__(_class, k, genset, prefix):
        return PrintingTerm.__new__(_class, k, genset, None, prefix)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset, prefix):
        return GrothendieckPoly.__xnew__(_class, k, genset, prefix)

    def _split_key(self):
        if isinstance(self._key, tuple) and len(self._key) == 2 and isinstance(self._key[0], Permutation) and isinstance(self._key[1], int):
            return self._key[0], self._key[1]
        return self._key, None

    def _sympystr(self, printer):
        perm, numvars = self._split_key()
        if perm == Permutation([]):
            return printer.doprint(ssymb.S.One)
        if numvars is None:
            return printer.doprint(f"{self._prefix}G{self._genset}({printer.doprint(perm)})")
        return printer.doprint(f"{self._prefix}G{self._genset}({printer.doprint(perm)}, {numvars})")

    def _pretty(self, printer):
        perm, numvars = self._split_key()
        if perm == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(perm)
        if numvars is None:
            return printer._print_Function(ssymb.Function(f"{self._prefix}{self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(self._genset)))
        return printer._print_Function(ssymb.Function(f"{self._prefix}{self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(f"{self._genset}; {numvars}")))

    def _latex(self, printer):
        perm, numvars = self._split_key()
        if perm == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(perm)
        if numvars is None:
            return printer._print_Function(ssymb.Function(f"{self._prefix}\\mathfrak{{G}}_{'{' + subscript + '}'}")(ssymb.Symbol(self._genset)))
        return printer._print_Function(ssymb.Function(f"{self._prefix}\\mathfrak{{G}}_{'{' + subscript + '}'}")(ssymb.Symbol(f"{self._genset}; {numvars}")))

    def __reduce__(self):
        return (self.__class__, (self._key, self._genset, self._prefix))


class DoubleGrothendieckPoly(PrintingTerm):
    """Double Grothendieck term ``DG<genset>(perm, <coeff_genset>)``."""

    is_Atom = True

    _pretty_schub_char = "𝔊"

    def __hash__(self):
        return hash((self._key, self._genset, self._coeff_genset, "DoubleGrothendieckPoly"))

    def __new__(cls, k, genset, coeff_genset, prefix=""):
        return DoubleGrothendieckPoly.__xnew_cached__(cls, k, genset, coeff_genset, prefix)

    @staticmethod
    def __xnew__(_class, k, genset, coeff_genset, prefix):
        return PrintingTerm.__new__(_class, k, genset, coeff_genset, prefix)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset, coeff_genset, prefix):
        return DoubleGrothendieckPoly.__xnew__(_class, k, genset, coeff_genset, prefix)

    def _sympystr(self, printer):
        key = self._key
        if key == Permutation([]):
            return printer.doprint(ssymb.S.One)
        return printer.doprint(f"{self._prefix}DG{self._genset}({printer.doprint(key)}, {self._coeff_genset})")

    def _pretty(self, printer):
        key = self._key
        if key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(key)
        return printer._print_Function(ssymb.Function(f"{self._prefix}{self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(f"{self._genset}; {self._coeff_genset}")))

    def _latex(self, printer):
        key = self._key
        if key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(key)
        return printer._print_Function(ssymb.Function(f"{self._prefix}\\mathfrak{{G}}_{'{' + subscript + '}'}")(ssymb.Symbol(f"{{{self._genset}}}; {{{self._coeff_genset}}}")))

    def __reduce__(self):
        return (self.__class__, (self._key, self._genset, self._coeff_genset, self._prefix))
