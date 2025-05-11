from functools import cache

import schubmult.symbolic as ssymb
from schubmult.perm_lib import Permutation

# from schubmult.rings.backend import expand, sympify
# from schubmult.symbolic import S


# Atomic Schubert polynomial
class AbstractSchubPoly(ssymb.Expr):
    is_Atom = True
    is_number = False

    def __new__(cls, k, genset, coeff_genset):
        obj = ssymb.Expr.__new__(cls)
        obj._key = k
        obj._genset = genset
        obj._coeff_genset = coeff_genset
        obj._perm = k
        return obj

    def __hash__(self):
        return hash((self._key,self._genset,self._coeff_genset,"AbS"))
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
        return ((*self._key,), self._genset, self._coeff_genset)


class DSchubPoly(AbstractSchubPoly):
    is_Atom = True

    _pretty_schub_char = "𝔖"  # noqa: RUF001

    def __hash__(self):
        return hash((self._key,self._genset,self._coeff_genset,"dasiub"))

    def __new__(cls, k, genset, coeff_genset):
        return DSchubPoly.__xnew_cached__(cls, k, genset, coeff_genset)

    @staticmethod
    def __xnew__(_class, k, genset, coeff_genset):
        return AbstractSchubPoly.__new__(_class, k, genset, coeff_genset)

    def _sympystr(self, printer):
        key = self._key
        if self._key == Permutation([]):
            return printer.doprint(ssymb.S.One)
        if self._coeff_genset is None:
            return printer.doprint(f"S{self._genset}({printer.doprint(key)})")
        return printer.doprint(f"DS{self._genset}({printer.doprint(key)}, {self._coeff_genset})")

    def _pretty(self, printer):
        key = self._key
        gl = self._genset
        if key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(int("".join([str(i) for i in key])))
        if self._coeff_genset is None:
            return printer._print_Function(ssymb.Function(f"{self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(gl)))
        return printer._print_Function(ssymb.Function(f"{self.__class__._pretty_schub_char}_{subscript}")(ssymb.Symbol(f"{self._genset}; {self._coeff_genset}")))

    def _latex(self, printer):
        key = self._key
        gl = self._genset
        if key == Permutation([]):
            return printer._print(ssymb.S.One)
        subscript = printer._print(key)
        if self._coeff_genset is None:
            return printer._print_Function(ssymb.Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(ssymb.Symbol(gl)))
        return printer._print_Function(ssymb.Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(ssymb.Symbol(f"{{{self._genset}}}; {{{self._coeff_genset}}}")))

    def __reduce__(self):
        return (self.__class__, self.args)

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset, coeff_genset):
        return DSchubPoly.__xnew__(_class, k, genset, coeff_genset)


class QDSchubPoly(AbstractSchubPoly):
    is_Atom = True

    def __hash__(self):
        return hash((self._key,self._genset,self._coeff_genset,"AbSsf","q"))

    _pretty_schub_char = "𝕼𝔖"  # noqa: RUF001

    def __new__(cls, k, genset, coeff_genset):
        return QDSchubPoly.__xnew_cached__(cls, k, genset, coeff_genset)

    @staticmethod
    def __xnew__(_class, k, genset, coeff_genset):
        return AbstractSchubPoly.__new__(_class, k, genset, coeff_genset)

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


class PQDSchubPoly(AbstractSchubPoly):
    is_Atom = True

    def __hash__(self):
        return hash((self._key,self._genset,self._coeff_genset,self.index_comp,"asfafAbS"))

    _pretty_schub_char = "𝕼𝔖"  # noqa: RUF001

    def __new__(cls, k, basis, index_comp):
        return PQDSchubPoly.__xnew_cached__(cls, k, basis, index_comp)

    @staticmethod
    def __xnew__(_class, k, basis, index_comp):
        obj = AbstractSchubPoly.__new__(_class, k, basis)
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
    def __xnew_cached__(_class, k, basis, index_comp):
        return PQDSchubPoly.__xnew__(_class, k, basis, index_comp)

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
