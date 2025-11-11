from functools import cached_property

from schubmult.rings.abstract_schub_poly import GenericPrintingTerm

# If you use symbolic coefficients like S.One, import S from your symbolic module
from schubmult.schub_lib.perm_lib import Permutation
from schubmult.symbolic import Add, Mul, S

# Utility for combining dictionaries with permutation keys
from schubmult.utils.perm_utils import add_perm_dict

# If you use _tensor_product_of_dicts_first in coproduct
# If you use self.monomial_basis and transition_schubert, import them from their respective modules
# If transition_schubert is a method of PolynomialBasis or another class, import accordingly
# For example, if it's a method of PolynomialBasis, you don't need to import it separately
# If you use TensorRing in change_tensor_basis
from .schubert_ring import Sx


class PolynomialBasis:

    @property
    def genset(self): ...

    def is_key(self, x): ...

    def as_key(self, x): ...

    def attach_key(self, dct):
        return {self.as_key(k): v for k, v in dct.items()}


    @property
    def zero_monom(self): ...

    @property
    def monomial_basis(self): ...

    # def product(self, key1, key2, coeff=S.One):
    #     sb = SchubertPolyBasis(numvars=self.numvars, ring=self.ring)
    #     dct = self.transition(sb)({key1: coeff})
    #     dct2 = self.transition(sb)({key2: S.One})
    #     res_dict = {}
    #     for k, v in dct.items():
    #         for k2, v2 in dct2.items():
    #             res_dict = add_perm_dict(res_dict, sb.product(k, k2, v * v2))
    #     return res_dict

    # @classmethod
    # def coproduct(self, key, coeff=S.One): ...

    def transition(self, other_basis): ...

    def printing_term(self, k): ...

    # @property
    # def zero_monom(self): ...

    @staticmethod
    def compose_transition(tkeyfunc, output):
        ret = {}
        for key, v in output.items():
            ret = add_perm_dict(ret, {k: v * v0 for k, v0 in tkeyfunc(key).items()})
        return ret

    def change_tensor_basis(self, tensor_elem, basis1, basis2):
        from .tensor_ring import TensorRing

        ring1 = tensor_elem.ring.rings[0]
        ring2 = tensor_elem.ring.rings[1]
        Tring2 = TensorRing(ring1.__class__(basis=basis1), ring2.__class__(basis=basis2))
        res = Tring2.zero
        for (key1, key2), v in tensor_elem.items():
            new_elem1 = ring1(*key1).change_basis(basis1)
            new_elem2 = ring2(*key2).change_basis(basis2)
            res += v * Tring2.ext_multiply(new_elem1, new_elem2)
        return res

    # @classmethod
    # @cache
    # def coproduct(self, key):
    #     from ._mul_utils import _tensor_product_of_dicts_first

    #     return PolynomialBasis.compose_transition(
    #         lambda x: _tensor_product_of_dicts_first(self.monomial_basis.transition(self)(x[0]), self.monomial_basis.transition(self)(x[1])),
    #         PolynomialBasis.compose_transition(lambda y: self.monomial_basis.coproduct(y), self.transition(self.monomial_basis)(key)),
    #     )

    # @classmethod
    # @cache

    def expand(self, dct):
        return self.monomial_basis.expand(self.transition(self.monomial_basis)(dct))

    def coproduct(self, key):
        from ._mul_utils import _tensor_product_of_dicts_first

        monom_version = self.transition(self.monomial_basis)({key: S.One})
        ret_dict = {}
        for k, v in monom_version.items():
            this_coproduct = {k0: v0 * v for k0, v0 in self.monomial_basis.coproduct(k).items()}
            for (mk1, mk2), v2 in this_coproduct.items():
                ret_dict = add_perm_dict(
                    ret_dict,
                    _tensor_product_of_dicts_first(self.monomial_basis.transition({mk1: v2}), self.monomial_basis.transition({mk2: S.One})),
                )
        return ret_dict

    def product(self, key1, key2, coeff=S.One):
        mnb = MonomialBasis(numvars=self.numvars, genset=self.genset)
        left = self.transition(mnb)({key1: coeff})
        right = self.transition(mnb)({key2: S.One})

        ret = {}

        for key_schub_right, v in right.items():
            # print(f"{key_schub_right=} {v=}")
            for key_schub_left, v2 in left.items():
                # print(f"{key_schub_left=} {v2=}")
                ret = add_perm_dict(ret, mnb.transition(self)(mnb.product(key_schub_left, key_schub_right, v * v2)))
        return ret


class MonomialBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        # if len(x) > self.numvars:
        #     raise ValueError(f"Length of key {x} exceeds the number of variables {self.numvars}.")
        # return (*x, *([0] * (self.numvars - len(x))))
        return tuple(x)

    # def attach_key(self, dct):
    #     return dct

    def printing_term(self, k):
        return GenericPrintingTerm(str(self.expand_monom(k)), "")

    def coproduct(self, key):
        result_dict = {}
        key = self.as_key(key)
        for i in range(self.numvars + 1):
            result_dict[(key[:i], key[i:])] = S.One
        return result_dict

    @property
    def monomial_basis(self):
        return self

    # @property
    # def numvars(self):
    #     return self._numvars

    @property
    def genset(self):
        return self._genset

    def __init__(self, genset):
        self._genset = genset

    def product(self, key1, key2, coeff=S.One):
        # print(f"{key1=} {key2=} {coeff=}")
        # if len(key1) != len(key2):
        #     return {}
        # if len(key1) != self.numvars:
        #     _altbasis = self.with_numvars(key1[1])
        # else:
        #     _altbasis = self
        if len(key1) != len(key2):
            return {}
        return {tuple(a + b for a, b in zip(key1, key2)): coeff}

    def expand_monom(self, monom):
        return Mul(*[self.genset[i + 1] ** monom[i] for i in range(len(monom))])

    def expand(self, dct):
        return Add(*[v * self.expand_monom(k) for k, v in dct.items()])

    def transition(self, other_basis):
        if isinstance(other_basis, MonomialBasis):
            return lambda x: other_basis.attach_key(x)
        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: other_basis.attach_key(other_basis.ring.from_expr(Add(*[v * self.expand_monom(k) for k, v in x.items()])))
        if isinstance(other_basis, SepDescPolyBasis):
            bonky_basis = SchubertPolyBasis(numvars=other_basis.numvars, ring=other_basis.ring)
            return lambda x: other_basis.attach_key(bonky_basis.transition(other_basis)(bonky_basis.attach_key(bonky_basis.ring.from_expr(Add(*[v * self.expand_monom(k) for k, v in x.items()])))))
        if isinstance(other_basis, ElemSymPolyBasis):
            spb = SchubertPolyBasis(numvars=other_basis.numvars, ring=other_basis.ring)
            return lambda x: spb.transition(other_basis)(self.transition(spb)(x))
        return None

    def from_expr(self, expr):
        from .variables import genset_dict_from_expr

        return self.attach_key(genset_dict_from_expr(expr, self.genset))

    @property
    def zero_monom(self):
        return self.as_key([])


class SchubertPolyBasis(PolynomialBasis):
    def __hash__(self):
        return hash(("schooboo", self.ring))

    # def with_numvars(self, numvars):
    #     return self.__class__(ring=self.ring)

    # @property
    # def numvars(self):
    #     return self._numvars

    def coproduct(self, key):
        length = key[1]
        res = {}
        for len0 in range(length + 1):
            cprd = dict(self.ring(key[0]).coproduct(*list(range(1, len0 + 1))))
            res = add_perm_dict(res, {((k1, len0), (k2, length - len0)): v for (k1, k2), v in cprd.items()})
        return res

    @property
    def genset(self):
        return self.ring.genset

    def printing_term(self, k):
        return self.ring.printing_term(k)

    def is_key(self, x):
        return isinstance(x, list | tuple | Permutation)

    def as_key(self, x):
        return (Permutation(x), len(Permutation(x).trimcode))

    def __init__(self, ring=None):
        # self._numvars = numvars
        self.ring = ring
        if self.ring is None:
            self.ring = Sx([]).ring

    @cached_property
    def monomial_basis(self):
        return MonomialBasis(genset=self.ring.genset)

    def product(self, key1, key2, coeff=S.One):
        if key1[1] != key2[1]:
            return {}
        def pair_length(dct, length):
            return {(k, length): v for k, v in dct.items()}

        return pair_length(self.ring.mul(self.ring.from_dict({key1[0]: coeff}), self.ring(key2[0])), key1[1])

    def transition_sepdesc(self, dct, other_basis):
        from schubmult.abc import e

        k = other_basis.k
        if k > other_basis.numvars:
            k = other_basis.numvars
        res_dict = {}

        dct2 = self.transition_elementary(dct, ElemSymPolyBasis(numvars=other_basis.numvars, ring=other_basis.ring))
        for (part, n), coeff0 in dct2.items():
            part1 = self.ring.one
            part2 = self.ring.one
            for i in range(k - 1):
                try:
                    part1 *= e(part[i], i + 1, self.ring.genset[1:])
                except IndexError:
                    pass
            for i in range(k - 1, len(part)):
                try:
                    part2 *= e(part[i], n if i + 1 > n else i + 1, self.ring.genset[1:])
                except IndexError:
                    pass
            for v, coeff1 in part1.items():
                for u, coeff2 in part2.items():
                    res_dict[other_basis.as_key([u, v])] = res_dict.get(other_basis.as_key([u, v]), S.Zero) + coeff1 * coeff2 * coeff0
        # for (w, n), coeff in dct.items():
        #     if len(w) <= k:
        #         res_dict[other_basis.as_key(Permutation([]), w)] = res_dict.get(other_basis.as_key(Permutation([]), w), S.Zero) + coeff
        #     else:
        #         cd = list(range(other_basis.numvars, 0, -1))
        #         if len(w) > other_basis.numvars:
        #             cd = ([cd[0]] * (len(w) - other_basis.numvars)) + cd
        #         mu = ~uncode(cd)

        #         mu0 = ~uncode(trimcode(~mu)[:-k+1])
        #         mu1 = ~uncode(trimcode(~mu)[-k+1:])
        #         # Lebniz
        return res_dict

    def transition_elementary(self, dct, other_basis):
        from schubmult.rings import FA
        from schubmult.schub_lib.perm_lib import uncode

        def elem_func(p, k, *_):
            return FA(p, k)

        elem = self.ring.from_dict(dct)
        res = {}
        for k, v in elem.items():
            if k[0].inv != 0:
                # d = other_basis.numvars - 1
                cd = list(range(other_basis.numvars, 0, -1))
                if len(k[0]) > other_basis.numvars:
                    cd = ([cd[0]] * (len(k[0]) - other_basis.numvars)) + cd
                w0 = uncode(cd)
                # print(f"{w0.code=} {k[0].code=}")
                rp = self.ring(k[0]).cem_rep(elem_func=elem_func, mumu=w0)
                for part, v2 in rp.items():
                    # part0 = part[2::2]
                    funny_bacon = [0]
                    for i in range(0, len(part), 2):
                        degree = part[i]
                        numvars = part[i + 1]
                        if numvars == 0:
                            continue
                        if numvars == other_basis.numvars:
                            if len(funny_bacon) < numvars:
                                funny_bacon += [0] * (numvars - len(funny_bacon))
                                funny_bacon[numvars - 1] = degree
                            else:
                                funny_bacon.append(degree)
                        else:
                            if len(funny_bacon) < numvars:
                                funny_bacon += [0] * (numvars - len(funny_bacon))
                            funny_bacon[numvars - 1] += degree
                    funny_bacon = funny_bacon[: other_basis.numvars - 1] + sorted(funny_bacon[other_basis.numvars - 1 :])
                    key = other_basis.as_key(funny_bacon)
                    res[key] = res.get(key, S.Zero) + v * v2
            else:
                key = other_basis.zero_monom
                res[key] = res.get(key, S.Zero) + v
        return res

    @property
    def zero_monom(self):
        return (Permutation([]), 0)

    def from_expr(self, expr):
        return self.attach_key(self.ring.from_expr(expr))

    def to_monoms(self, key):
        from .variables import genset_dict_from_expr

        def pad_tuple(tup, length):
            return (*tup, *(0,) * (length - len(tup)))

        dct = {pad_tuple(k, key[1]): v for k, v in genset_dict_from_expr(self.ring.from_dict({key[0]: S.One}).as_polynomial(), self.genset).items()}
        return dct

    def transition(self, other_basis):
        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: other_basis.attach_key({k[0]: v for k, v in x.items()})
        if isinstance(other_basis, MonomialBasis):
            def sum_dct(*dcts):
                res = {}
                for dct in dcts:
                    res = add_perm_dict(res, dct)
                return res
            def mul_dict(dct, coeff):
                return {k: v * coeff for k, v in dct.items()}

            return lambda x: sum_dct(*[mul_dict(self.to_monoms(k),v) for k, v in x.items()])
        if isinstance(other_basis, ElemSymPolyBasis):
            return lambda x: self.transition_elementary(x, other_basis)
        if isinstance(other_basis, SepDescPolyBasis):
            return lambda x: self.transition_sepdesc(x, other_basis)
            # raise NotImplementedError("The bonky cheese need to implement the elemnify")
        return None


class SepDescPolyBasis(PolynomialBasis):
    def __hash__(self):
        return hash(self.numvars, self.ring, self.k, "fatbacon")

    @property
    def zero_monom(self):
        return self.as_key([[], []])

    def with_numvars(self, numvars):
        return self.__class__(numvars=numvars, ring=self.ring, k=self.k)

    def product(self, key1, key2, coeff=S.One):
        mnb = SchubertPolyBasis(numvars=self.numvars, ring=self.ring)
        left = self.transition(mnb)({key1: coeff})
        right = self.transition(mnb)({key2: S.One})

        ret = {}

        for key_schub_right, v in right.items():
            # print(f"{key_schub_right=} {v=}")
            for key_schub_left, v2 in left.items():
                # print(f"{key_schub_left=} {v2=}")
                ret = add_perm_dict(ret, mnb.transition(self)(mnb.product(key_schub_left, key_schub_right, v * v2)))
        return ret

    @property
    def numvars(self):
        return self._numvars

    @property
    def genset(self):
        return self.ring.genset

    @property
    def k(self):
        return self._k

    def printing_term(self, k):
        return GenericPrintingTerm(k, "K")

    def is_key(self, x):
        return isinstance(x, list | tuple)

    def as_key(self, x):
        return (Permutation(x[0]), Permutation(x[1]), self.k, self.numvars)

    def __init__(self, numvars, k, ring=None):
        self._numvars = numvars
        self._k = k
        self.ring = ring
        if self.ring is None:
            self.ring = Sx([]).ring

    @cached_property
    def monomial_basis(self):
        return MonomialBasis(numvars=self.numvars, genset=self.ring.genset)

    # def product(self, key1, key2, coeff=S.One):
    #     if key1[1] != key2[1]:
    #         return {}
    #     if key1[1] != self.numvars:
    #         _altbasis = self.with_numvars(key1[1])
    #     else:
    #         _altbasis = self
    #     return _altbasis.attach_key(self.ring.mul(self.ring.from_dict({key1[0]: coeff}), self.ring(key2[0])))

    def transition_schubert(self, dct, other_basis):
        # print(f"{dct=}")
        # from schubmult.schub_lib.perm_lib import uncode
        # from schubmult.rings import FA

        # def elem_func(p, k, *_):
        #     return FA(p, k)

        # elem = self.ring.from_dict(dct)
        # res = {}
        # for k, v in elem.items():
        #     if k[0].inv != 0:
        #         # d = other_basis.numvars - 1
        #         cd = list(range(other_basis.numvars, 0, -1))
        #         if len(k[0]) > other_basis.numvars:
        #             cd = ([cd[0]] * (len(k[0]) - other_basis.numvars)) + cd
        #         # spot = len(cd) - d - 1
        #         # cd = [cd[spot]] * (spot) + cd[spot:]
        #         # print(f"{cd=}")
        #         w0 = uncode(cd)
        #         # print(f"{w0.code=} {k[0].code=}")
        #         rp = self.ring(k[0]).cem_rep(elem_func=elem_func, mumu=w0)
        #         for part, v2 in rp.items():
        #             # part0 = part[2::2]
        #             funny_bacon = [0]
        #             for i in range(0, len(part), 2):
        #                 degree = part[i]
        #                 numvars = part[i + 1]
        #                 if numvars == 0:
        #                     continue
        #                 if numvars == other_basis.numvars:
        #                     if len(funny_bacon) < numvars:
        #                         funny_bacon += [0] * (numvars - len(funny_bacon))
        #                         funny_bacon[numvars - 1] = degree
        #                     else:
        #                         funny_bacon.append(degree)
        #                 else:
        #                     if len(funny_bacon) < numvars:
        #                         funny_bacon += [0] * (numvars - len(funny_bacon))
        #                     funny_bacon[numvars - 1] += degree
        #             funny_bacon = funny_bacon[: other_basis.numvars - 1] + sorted(funny_bacon[other_basis.numvars - 1 :])
        #             key = other_basis.as_key(funny_bacon)
        #             res[key] = res.get(key, S.Zero) + v * v2
        #     else:
        #         key = other_basis.zero_monom
        #         res[key] = res.get(key, S.Zero) + v
        # return res
        out_ret = self.ring.zero
        for (k1, k2, _, __), v in dct.items():
            out_ret += self.ring.from_dict({k1: v}) * self.ring.from_dict({k2: S.One})
        return other_basis.attach_key(dict(out_ret))

    def transition_sepdesc(x, other_basis):  # noqa: ARG002
        return {x: S.One}

    def from_expr(self, expr):
        return self.attach_key(self.ring.from_expr(expr))

    def transition(self, other_basis):
        if isinstance(other_basis, SepDescPolyBasis):
            return lambda x: self.transition_sepdesc(x, other_basis)
        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: self.transition_schubert(x, other_basis)
        if isinstance(other_basis, MonomialBasis):
            bonky_basis = SchubertPolyBasis(self.numvars, self.ring)
            return lambda x: bonky_basis.transition(other_basis)(self.transition_schubert(x, bonky_basis))
        if isinstance(other_basis, ElemSymPolyBasis):
            return lambda x: other_basis.transition_elementary(self.transition_schubert(x, other_basis))
            # raise NotImplementedError("The bonky cheese need to implement the elemnify")
        return None


class ElemSymPolyBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return ((*x,), self.numvars)

    def printing_term(self, k):
        return GenericPrintingTerm(k, "E")

    def __hash__(self):
        return hash((self.numvars, self.ring, "dangit_bobjo"))

    def with_numvars(self, numvars):
        return self.__class__(numvars=numvars, ring=self.ring)

    @cached_property
    def monomial_basis(self):
        return MonomialBasis(numvars=self.numvars, genset=self.ring.genset)

    @property
    def numvars(self):
        return self._numvars

    @property
    def genset(self):
        return self.ring.genset

    def __init__(self, numvars, ring=None):
        self._numvars = numvars
        self.ring = ring
        if self.ring is None:
            self.ring = Sx([]).ring

    def transition_schubert(self, dct):
        from schubmult.abc import e

        res = self.ring.zero
        for (k, n), v in dct.items():
            to_add = self.ring.one
            for i, a in enumerate(k[:n]):
                to_add = self.ring.elem_mul(to_add, e(a, i + 1, self.ring.genset[1:]))
            for a in enumerate(k[n:]):
                to_add = self.ring.elem_mul(to_add, e(a, n, self.ring.genset[1:]))
            res += v * to_add
        return res

    def transition_monomial(self, dct):
        from schubmult.abc import e
        from schubmult.symbolic import expand_func

        res = S.Zero
        for (k, n), v in dct.items():
            # print(f"{(k,n)=} {v=}")
            to_add = S.One
            for i, a in enumerate(k[:n]):
                # if a > i +1:
                # print(f"{a=} {i=} a>i+1 waffle")
                to_add *= expand_func(e(a, i + 1, self.ring.genset[1:]))
                # print(f"{to_add=}")
            for a in k[n:]:
                to_add *= expand_func(e(a, n, self.ring.genset[1:]))
                # print(f"{to_add=    }")
            res += v * to_add
        # print(f"{res=}")
        return res

    def transition(self, other_basis):
        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: other_basis.attach_key(self.transition_schubert(x))
        if isinstance(other_basis, MonomialBasis):
            from .variables import genset_dict_from_expr

            return lambda x: other_basis.attach_key(genset_dict_from_expr(self.transition_monomial(x), other_basis.genset))
        if isinstance(other_basis, ElemSymPolyBasis):
            return lambda x: x
        return None

    @property
    def zero_monom(self):
        return self.as_key([])
