from schubmult.symbolic import CoercionFailed, S, sympify, sympy_Mul

from .base_ring import BaseRing, BaseRingElement
from .printing import PrintingTerm


class DirectProductRing(BaseRing):
    """Direct product of an arbitrary (fixed) number of rings.

    An element is stored as a dict mapping keys ``(i, k)`` to coefficients,
    where ``i`` is the component index and ``k`` is a basis key from the
    *i*-th ring.  Addition and multiplication are componentwise: terms from
    different components never interact, and cross-component products are
    zero.

    Parameters
    ----------
    ``*rings`` : BaseRing
        One or more constituent rings.

    Examples
    --------
    >>> D = DirectProductRing(R1, R2, R3)
    >>> a = D.from_component(0, some_R1_element)
    >>> b = D.from_component(2, some_R3_element)
    >>> a + b          # lives in components 0 and 2
    >>> a * b          # zero (different components)
    >>> D[0]           # R1
    >>> D.project(a, 0)  # back to R1
    """

    def __eq__(self, other):
        return type(self) is type(other) and self.rings == other.rings

    def __hash__(self):
        return hash(("DirectProductRing", self.rings))

    def __init__(self, *rings):
        if len(rings) == 0:
            raise ValueError("DirectProductRing requires at least one ring")
        super().__init__()
        self._rings = tuple(rings)
        self.zero_monom = (0, self._rings[0].zero_monom)

    def __getitem__(self, i):
        """Return the *i*-th constituent ring."""
        return self._rings[i]

    def __len__(self):
        return len(self._rings)

    def dtype(self):
        elem = DirectProductRingElement()
        elem.ring = self
        return elem

    @property
    def rings(self):
        return self._rings

    @property
    def one(self):
        dct = {}
        for i, ring_i in enumerate(self._rings):
            dct[(i, ring_i.zero_monom)] = S.One
        return self.from_dict(dct)

    def component_one(self, i):
        """Return the identity element supported only on component *i*."""
        return self.from_dict({(i, self._rings[i].zero_monom): S.One})

    def from_dict(self, element):
        dct = self.dtype()
        dct.update({k: v for k, v in element.items() if v != 0})
        return dct

    def from_component(self, i, elem):
        """Lift an element of ``self[i]`` into the direct product."""
        return self.from_dict({(i, k): v for k, v in elem.items()})

    def project(self, elem, i):
        """Project onto component *i*, returning an element of ``self[i]``."""
        ring_i = self._rings[i]
        return ring_i.from_dict({k: v for (j, k), v in elem.items() if j == i})

    def rmul(self, elem, scalar):
        return self.from_dict({k: v * scalar for k, v in elem.items()})

    def mul(self, elem1, elem2):
        """Componentwise multiplication.

        Only terms in the same component interact; cross-component products
        are zero.
        """
        result = {}
        for (i1, k1), v1 in elem1.items():
            for (i2, k2), v2 in elem2.items():
                if i1 != i2:
                    continue
                ring_i = self._rings[i1]
                prod_i = ring_i.from_dict({k1: S.One}) * ring_i.from_dict({k2: S.One})
                for k, v in prod_i.items():
                    key = (i1, k)
                    result[key] = result.get(key, S.Zero) + v1 * v2 * v
        return self.from_dict(result)

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        if isinstance(element, BaseRingElement):
            raise CoercionFailed("Not a domain element")
        return sympify(element)

    def printing_term(self, k):
        return DirectProductBasisElement(k, self)

    def __call__(self, *elems):
        """Construct an element from one element per component.

        ``D(e0, e1, ..., en)`` lifts each ``ei`` (an element of ``D[i]``)
        into the direct product and sums them.
        """
        if len(elems) != len(self._rings):
            raise ValueError(
                f"Expected {len(self._rings)} elements (one per component), got {len(elems)}",
            )
        result = self.zero
        for i, elem in enumerate(elems):
            if isinstance(elem, BaseRingElement):
                result += self.from_component(i, elem)
            else:
                result += self.from_dict({(i, self._rings[i].zero_monom): sympify(elem)})
        return result


class DirectProductBasisElement(PrintingTerm):
    is_commutative = False
    precedence = 50

    def __new__(cls, k, ring):
        return DirectProductBasisElement.__xnew__(cls, k, ring)

    @staticmethod
    def __xnew__(_class, k, ring):
        obj = PrintingTerm.__new__(_class, k, None, None)
        obj._key = k
        obj.ring = ring
        return obj

    def __hash__(self):
        return hash((self._key, self.ring))

    def _sympystr(self, printer):
        i, inner_key = self._key
        inner_term = self.ring[i].printing_term(inner_key)
        return f"({printer._print(inner_term)})_{i}"

    def _pretty(self, printer):
        i, inner_key = self._key
        inner_term = self.ring[i].printing_term(inner_key)
        return printer._print(inner_term)

    def _latex(self, printer):
        i, inner_key = self._key
        inner_term = self.ring[i].printing_term(inner_key)
        return f"({printer._print(inner_term)})_{{{i}}}"


class DirectProductRingElement(BaseRingElement):
    def __init__(self):
        pass

    def __getitem__(self, key):
        """Index by component integer or by ``(i, k)`` basis key.

        * ``elem[i]`` — project onto component *i* (returns a ``self.ring[i]`` element).
        * ``elem[(i, k)]`` — coefficient lookup (standard dict behaviour).
        """
        if isinstance(key, int):
            return self.ring.project(self, key)
        return dict.__getitem__(self, key)

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [
            self[(i, k)] if (i, k) == self.ring.zero_monom else sympy_Mul(self[(i, k)], self.ring.printing_term((i, k)))
            for i, k in sorted(dict.keys(self))
        ]

    @property
    def free_symbols(self):
        ret = set()
        for (i, k), v in self.items():
            if hasattr(v, 'free_symbols'):
                ret.update(v.free_symbols)
            ring_elem = self.ring[i].from_dict({k: S.One})
            if hasattr(ring_elem, 'free_symbols'):
                ret.update(ring_elem.free_symbols)
        return ret
