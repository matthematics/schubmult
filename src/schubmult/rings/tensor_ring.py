from functools import cache

from sympy import Tuple

from schubmult.symbolic import Mul, S, sympy_Mul
from schubmult.utils.logging import get_logger

from .printing import PrintingTerm
from .schubert.base_schubert_ring import BaseSchubertElement, BaseSchubertRing

logger = get_logger(__name__)


class TensorRing(BaseSchubertRing):
    # tensor ring
    def __eq__(self, other):
        return type(self) is type(other) and self.rings == other.rings

    @property
    def args(self):
        return ()

    def coproduct_on_basis(self, k):
        """Compute coproduct of a basis element in the tensor ring.

        For k = (k_1, k_2, ..., k_n) in ring_1 ⊗ ring_2 ⊗ ... ⊗ ring_n,
        Δ(k_1 ⊗ k_2 ⊗ ... ⊗ k_n) = (⊗ Δ(k_i))

        This properly interlaces the individual coproducts.
        """
        tring = self @ self
        assert len(tring.rings) == len(self.rings) * 2, "Coproduct ring should have twice as many factors as original ring"
        # Get all individual coproducts
        coprods = [self.rings[i].coproduct_on_basis(k[i]) for i in range(len(self.rings))]

        # Build result by taking all combinations
        result_dict = {}

        # For each combination of left/right factors from each coproduct
        def recurse_coproducts(idx, left_parts, right_parts, coeff):
            if idx == len(coprods):
                # Base case: we've processed all factors
                # Create flat key: (k1_L, k2_L, ..., kn_L, k1_R, k2_R, ..., kn_R)
                key = tuple(left_parts + right_parts)
                assert len(key) == len(tring.rings), f"Key length {len(key)} != tring.rings length {len(tring.rings)}"
                if key in result_dict:
                    result_dict[key] += coeff
                else:
                    result_dict[key] = coeff
                return

            # Process coproduct of factor idx
            for (left_i, right_i), coeff_i in coprods[idx].items():
                # Ring is not a TensorRing after flattening, so keys are single elements
                left_to_add = [left_i]
                right_to_add = [right_i]
                recurse_coproducts(idx + 1, left_parts + left_to_add, right_parts + right_to_add, coeff * coeff_i)

        recurse_coproducts(0, [], [], 1)
        return tring.from_dict(result_dict)

    def from_rc_graph_tensor(self, rc_graph_tensor):
        return self.ext_multiply(self.rings[0].from_rc_graph(rc_graph_tensor[0]), self.rings[1].from_rc_graph(rc_graph_tensor[1]))

    # def sub(self, elem, other):
    #     print("Mooger")
    #     return self.from_dict(add_perm_dict(elem, {k: -v for k,v in other.items()}))

    def __init__(self, *rings):
        self._rings = rings
        ring_list = [*self._rings]
        while any(isinstance(r, TensorRing) for r in ring_list):
            new_ring_list = []
            for r in ring_list:
                if isinstance(r, TensorRing):
                    new_ring_list.extend(r.rings)
                else:
                    new_ring_list.append(r)
            ring_list = new_ring_list
        self._rings = tuple(ring_list)
        genset = set()
        for r in self._rings:
            try:
                genset.update(set(r.genset))
            except AttributeError:
                pass
        genset = tuple(genset)
        coeff_genset = set()
        for r in self._rings:
            try:
                if r.coeff_genset.label:
                    coeff_genset.update(set(r.coeff_genset))
            except AttributeError:
                pass
        coeff_genset = tuple(coeff_genset)
        super().__init__(list(genset), list(coeff_genset))
        self.zero_monom = tuple([self.rings[i].zero_monom for i in range(len(self.rings))])
        # self.dtype = type("TensorRingElement", (TensorRingElement,), {"ring": self})

    def dtype(self):
        elem = TensorRingElement()
        elem.ring = self
        return elem

    def __hash__(self):
        return hash(self.rings)

    @property
    def rings(self):
        return self._rings

    def rmul(self, elem1, elem2):
        # print(f"{dict(elem1)=} {elem2=}")
        # print(f"{self.zero_monom=} {type(elem2)=}")
        # return self.from_dict({self.zero_monom: elem2}) * elem1
        return self.from_dict({k: v * elem2 for k, v in elem1.items()})
        # except Exception:
        #     # import traceback
        #     # traceback.print_exc()
        #     raise

    def from_dict(self, element):
        dct = self.dtype()
        dct.update({k: v for k, v in element.items() if v != 0})
        return dct

    def mul(self, elem1, elem2):
        """Multiply two elements in the tensor ring.

        (a1 ⊗ a2 ⊗ ... ⊗ an) * (b1 ⊗ b2 ⊗ ... ⊗ bn) = (a1*b1) ⊗ (a2*b2) ⊗ ... ⊗ (an*bn)
        """
        ret_dict = {}

        for k1, v1 in elem1.items():
            for k2, v2 in elem2.items():
                # Compute products of each factor
                factor_products = []
                for i in range(len(self.rings)):
                    prod_i = self.rings[i].from_dict({k1[i]: S.One}) * self.rings[i].from_dict({k2[i]: S.One})
                    factor_products.append(prod_i)

                # Now take all combinations of terms from each factor product
                def recurse_terms(idx, current_key, current_coeff):
                    if idx == len(factor_products):
                        # We have a complete key
                        key = tuple(current_key)
                        total_coeff = v1 * v2 * current_coeff
                        if key in ret_dict:
                            ret_dict[key] += total_coeff
                        else:
                            ret_dict[key] = total_coeff
                        return

                    # Process all terms from factor idx
                    for key_i, coeff_i in factor_products[idx].items():
                        recurse_terms(idx + 1, [*current_key, key_i], current_coeff * coeff_i)

                recurse_terms(0, [], S.One)

        return self.from_dict(ret_dict)

    def _coerce_add(self, x):  # noqa: ARG002
        return None

    def _coerce_mul(self, x):
        # have to pull out gens
        if not isinstance(x.ring, TensorRing):
            if set(x.ring.genset) == set(self.genset):
                return x.coproduct(*[x.ring.genset.index(v) for v in self.rings[0].genset[1:]])
        return None

    @property
    def one(self):
        return self.from_dict({self.zero_monom: S.One})

    @cache
    def cached_schubpoly(self, k):
        return Mul(*[self.rings[i].cached_schubpoly(k[i]) for i in range(len(self.rings))])

    def printing_term(self, k):
        return TensorBasisElement(k, self)

    def from_comp_ring(self, t):
        dct = {}
        for k, v in t.items():
            new_k = list(self.zero_monom)
            if isinstance(t.ring, TensorRing):
                for i in range(len(t.ring.rings)):
                    new_k[self.rings.index(t.ring.rings[i])] = k[i]
            else:
                new_k[self.rings.index(t.ring)] = k
            dct[tuple(new_k)] = v
        return self.from_dict(dct)

    # def from_sympy(self, x):
    #     dct = {}
    #     elem1 = self.rings[0].from_sympy(x)
    #     dct = {(k,): v for k,v in elem1.items()}
    #     for i in range(1, len(self.rings)):
    #         dct_new = {}
    #         for k, v in dct.items():
    #             dct_new = add_perm_dict(dct_new,{(*k, k1): v1 for k1, v1 in self.rings[i].from_sympy(v).items()})
    #         dct = dct_new
    #     return self.from_dict(dct)

    def ext_multiply(self, elem1, elem2):
        ret = self.zero

        for key, val in elem1.items():
            for key2, val2 in elem2.items():
                if isinstance(elem1.ring, TensorRing):
                    if isinstance(elem2.ring, TensorRing):
                        ret += self.from_dict({(*key, *key2): val * val2})
                    else:
                        ret += self.from_dict({(*key, key2): val * val2})
                else:
                    if isinstance(elem2.ring, TensorRing):
                        ret += self.from_dict({(key, *key2): val * val2})
                    else:
                        ret += self.from_dict({(key, key2): val * val2})
        return ret

    def __call__(self, x):
        if isinstance(x, tuple):
            return self.from_dict({x: self.domain.one})
        return self.from_expr(x)


class TensorBasisElement(PrintingTerm):
    is_commutative = False
    precedence = 50

    def __new__(cls, k, basis):
        return TensorBasisElement.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj = PrintingTerm.__new__(_class, k, None, None)
        obj._key = k
        obj.ring = basis
        assert len(k) == len(basis.rings), f"Key length must match number of factors in tensor ring, got {len(k)} and {len(basis.rings)}"
        return obj

    def __hash__(self):
        return hash((self._key, self.ring))

    @staticmethod
    def __xnew_cached__(_class, k, basis):
        return TensorBasisElement.__xnew__(_class, k, basis)

    def _sympystr(self, printer):
        return " # ".join([printer._print(self.ring.rings[i].printing_term(self._key[i])) for i in range(len(self._key))])

    def _pretty(self, printer):
        return printer._print_TensorProduct(Tuple(*[self.ring.rings[i].printing_term(self._key[i]) for i in range(len(self._key))]))

    def _latex(self, printer):
        return printer._print_TensorProduct(Tuple(*[self.ring.rings[i].printing_term(self._key[i]) for i in range(len(self._key))]))


class TensorRingElement(BaseSchubertElement):
    def __init__(self):
        pass

    def coproduct(self):
        """Override coproduct to use the correct target ring."""
        tring = self.ring @ self.ring
        result = tring.zero
        for k, v in self.items():
            result += v * self.ring.coproduct_on_basis(k)
        return result

    @property
    def free_symbols(self):
        ret = set()
        for k, v in self.items():
            ret.update(v.free_symbols)
            for i in range(len(k)):
                ret.update(self.ring.rings[i](k[i]).free_symbols)
        return ret

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [
            self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k))
            for k in self.keys()  # sorted(self.keys(), key=lambda kkt: [(kk.inv, tuple(kk)) for kk in kkt])
        ]
