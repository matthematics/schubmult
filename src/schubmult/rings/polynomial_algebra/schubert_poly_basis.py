from functools import cached_property

from schubmult.combinatorial_reps.permutation import Permutation
from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict

from ..schubert.schubert_ring import Sx
from .base_polynomial_basis import PolynomialBasis


class SchubertPolyBasis(PolynomialBasis):
    def __hash__(self):
        return hash(("schooboo", self.ring))

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
        self.ring = ring
        if self.ring is None:
            self.ring = Sx([]).ring

    @cached_property
    def monomial_basis(self):
        from .monomial_basis import MonomialBasis

        return MonomialBasis(genset=self.ring.genset)

    def product(self, key1, key2, coeff=S.One):
        if key1[1] != key2[1]:
            return {}

        def pair_length(dct, length):
            return {(k, length): v for k, v in dct.items()}

        return pair_length(self.ring.mul(self.ring.from_dict({key1[0]: coeff}), self.ring(key2[0])), key1[1])

    def transition_sepdesc(self, dct, other_basis):
        from schubmult.abc import e

        from .elem_sym_poly_basis import ElemSymPolyBasis

        k = other_basis.k
        res_dict = {}

        dct2 = self.transition_elementary(dct, ElemSymPolyBasis(ring=other_basis.ring))
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
        return res_dict

    def transition_elementary(self, dct, other_basis):
        from schubmult.combinatorial_reps.permutation import uncode
        from schubmult.rings import FA

        res = {}
        elem = self.ring.from_dict(dct)

        def elem_func(p, k, *_):
            return FA(p, k)

        for k, v in elem.items():
            num_vars = k[1]
            if k[0].inv != 0:
                cd = list(range(num_vars, 0, -1))
                if len(k[0]) > num_vars:
                    cd = ([cd[0]] * (len(k[0]) - num_vars)) + cd
                w0 = uncode(cd)
                rp = self.ring(k[0]).cem_rep(elem_func=elem_func, mumu=w0)
                for part, v2 in rp.items():
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
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        def pad_tuple(tup, length):
            return (*tup, *(0,) * (length - len(tup)))

        dct = {pad_tuple(k, key[1]): v for k, v in genset_dict_from_expr(self.ring.from_dict({key[0]: S.One}).as_polynomial(), self.genset).items()}
        return dct

    def transition(self, other_basis):
        from .elem_sym_poly_basis import ElemSymPolyBasis
        from .monomial_basis import MonomialBasis
        from .sepdesc_poly_basis import SepDescPolyBasis

        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: dict(x)
        if isinstance(other_basis, MonomialBasis):
            def sum_dct(*dcts):
                res = {}
                for dct in dcts:
                    res = add_perm_dict(res, dct)
                return res

            def mul_dict(dct, coeff):
                return {k: v * coeff for k, v in dct.items()}

            return lambda x: sum_dct(*[mul_dict(self.to_monoms(k), v) for k, v in x.items()])
        if isinstance(other_basis, ElemSymPolyBasis):
            return lambda x: self.transition_elementary(x, other_basis)
        if isinstance(other_basis, SepDescPolyBasis):
            return lambda x: self.transition_sepdesc(x, other_basis)
        return None
