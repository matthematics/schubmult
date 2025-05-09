from functools import cache


import schubmult.symbolic as ssymb
from schubmult.perm_lib import Permutation

# from schubmult.rings.backend import expand, sympify
# from schubmult.symbolic import S

_pretty_schub_char = "𝔖"  # noqa: RUF001
# Atomic Schubert polynomial
class AbstractSchubPoly(ssymb.SymengineExpr):
    is_Atom = True

    def __new__(cls, k, genset, coeff_genset):
        obj = ssymb.SymengineExpr.__new__(cls)
        obj._key = k
        obj._genset = genset
        obj._coeff_genset = coeff_genset
        obj._perm = k
        return obj

    def __init__(self, k, genset, coeff_genset):
        super().__init__(k, genset, coeff_genset)
    
    # def __reduce__(self):
    #     return (self.__class__, self.args)

    # def _eval_expand_basic(self, *args, **kwargs):  # noqa: ARG002
    #     return self._ring_elem.as_polynomial()

    # def __reduce__(self):
    #     return (self.__class__, self.args)

    # def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
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
        return ((*self._key,), self._genset, self._coeff_genset)


class DSchubPoly(AbstractSchubPoly):
    is_Atom = True

    def __new__(cls, k, genset, coeff_genset):
        return DSchubPoly.__xnew_cached__(cls, k, genset, coeff_genset)

    @staticmethod
    def __xnew__(_class, k, genset, coeff_genset):
        return AbstractSchubPoly.__new__(_class, k, genset, coeff_genset)

    def _sympystr(self, printer):
        key = self._key
        if self._key == Permutation([]):
            return printer.doprint(1)
        if self._coeff_genset is None:
            return printer.doprint(f"S{self._genset}({printer.doprint(key)})")
        return printer.doprint(f"DS{self._genset}({printer.doprint(key)}, {self._coeff_genset})")

    def _pretty(self, printer):
        key = self._key
        gl = self._genset
        if key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(int("".join([str(i) for i in key])))
        if self._coeff_genset is None:
            return printer._print_Function(Function(f"{_pretty_schub_char}_{subscript}")(Symbol(gl)))
        return printer._print_Function(Function(f"{_pretty_schub_char}_{subscript}")(Symbol(f"{self._genset}; {self._coeff_genset}")))

    def _latex(self, printer):
        key = self._key
        gl = self._genset
        if key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(key)
        if self._coeff_genset is None:
            return printer._print_Function(Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(Symbol(gl)))
        return printer._print_Function(
            Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(Symbol(f"{{{self._genset}}}; {{{self._coeff_genset}}}"))
        )

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset, coeff_genset):
        return DSchubPoly.__xnew__(_class, k, genset, coeff_genset)
