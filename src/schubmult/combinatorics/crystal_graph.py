from itertools import zip_longest

from sympy.printing.defaults import Printable


class CrystalGraph(Printable):
    def raising_operator(self, index):
        """The raising operator for the crystal graph."""
        raise NotImplementedError

    def lowering_operator(self, index):
        """The lowering operator for the crystal graph."""
        raise NotImplementedError

    @property
    def crystal_weight(self):
        """The weight of the crystal graph element."""
        raise NotImplementedError

    def phi(self, i):
        if i < 1 or i >= self.crystal_length():
            if i - 1 < self.crystal_length():
                return self.epsilon(i) + self.crystal_weight[i - 1]
            return self.epsilon(i)
        return self.epsilon(i) + self.crystal_weight[i - 1] - self.crystal_weight[i]


    def epsilon(self, i):
        if i > 0 and i < self.crystal_length():
            cnt = 0
            rc = self
            while rc is not None:
                rc = rc.raising_operator(i)
                if rc is not None:
                    cnt += 1
            return cnt
        return 0

    def to_lowest_weight(self, length=None):
        """Return the lowest weight element in the connected component."""
        g = self
        lower_seq = []
        found = True
        if length is None:
            length = self.crystal_length()
        g = self.to_highest_weight(length=length)[0]
        while found:
            found = False
            for row in range(length - 1, 0, -1):
                g0 = g.lowering_operator(row)
                if g0 is not None:
                    found = True
                    g = g0
                    lower_seq.append(row)
                    break
        return (g, tuple(lower_seq))

    def to_highest_weight(self, length=None, start=1):
        """Return the highest weight element in the connected component."""
        g = self
        raise_seq = []
        found = True
        if length is None:
            length = self.crystal_length()
        while found:
            found = False
            for row in range(start, length):
                g0 = g.raising_operator(row)
                if g0 is not None:
                    found = True
                    g = g0
                    raise_seq.append(row)
                    break
        return (g, tuple(raise_seq))

    def crystal_length(self):
        """Return the length of the crystal element."""
        raise NotImplementedError

    def reverse_raise_seq(self, raise_seq):
        rc = self
        for row in reversed(raise_seq):
            rc = rc.lowering_operator(row)
            if rc is None:
                raise ValueError(f"Cannot reverse raising sequence {raise_seq} on {self}")
        return rc

    def reverse_lower_seq(self, lower_seq):
        rc = self
        for row in reversed(lower_seq):
            rc = rc.raising_operator(row)
            if rc is None:
                return None
        return rc

    def crystal_reflection(self, index):
        new_rc = self
        e = self.epsilon(index)
        f = self.phi(index)
        if e > f:
            for _ in range(e - f):
                new_rc = new_rc.raising_operator(index)
                if new_rc is None:
                    return self
        elif f > e:
            for _ in range(f - e):
                new_rc = new_rc.lowering_operator(index)
                if new_rc is None:
                    return self
        return new_rc


    @property
    def full_crystal(self):
        hw, _ = self.to_highest_weight()
        return hw.crystal_beneath

    @property
    def crystal_beneath(self):
        crystal = set()
        stack = [self]
        while len(stack) > 0:
            elem = stack.pop()
            crystal.add(elem)
            for i in range(1, elem.crystal_length()):
                new_elem = elem.lowering_operator(i)
                if new_elem is not None:
                    stack.append(new_elem)
        return crystal

    def crystal_above(self, length=None):
        crystal = set()
        stack = [self]
        if length is None:
            length = self.crystal_length()
        while len(stack) > 0:
            elem = stack.pop()
            crystal.add(elem)
            for i in range(1, length):
                new_elem = elem.raising_operator(i)
                if new_elem is not None:
                    stack.append(new_elem)
        return crystal

    def truncated_crystal(self, length, start=1):
        hw, _ = self.to_highest_weight(length=length)
        crystal = set()
        stack = [hw]
        while len(stack) > 0:
            elem = stack.pop()
            crystal.add(elem)
            for i in range(start, length):
                new_elem = elem.lowering_operator(i)
                if new_elem is not None:
                    stack.append(new_elem)
        return crystal

    @property
    def params(self):
        param_list = []
        rc = self
        for i in range(1, self.crystal_length()):
            ai = 0
            while rc.raising_operator(i) is not None:
                rc = rc.raising_operator(i)
                ai += 1
            param_list.append(ai)
        return tuple(param_list)

    def weight_bump(self):
        ...

    def weight_reflection(self, index):
        res = self
        try:
            res = res.crystal_reflection(index)
        except Exception:
            res = None
        if res is None:
            res = self.weight_bump()
            return res.crystal_reflection(index)
        return res

    @property
    def is_highest_weight(self):
        for row in range(1, self.crystal_length()):
            if self.raising_operator(row) is not None:
                # print(f"IKINPROVENOTHIGHESTWEIGHT {row=}")
                # print(self.raising_operator(row))
                return False
        return True

    @property
    def is_lowest_weight(self):
        for row in range(1, self.crystal_length()):
            if self.lowering_operator(row) is not None:
                return False
        return True

    @property
    def dual(self):
        return CrystalGraphDual(self)

    @property
    def reverse(self):
        return CrystalGraphReverse(self)

class CrystalGraphDual(CrystalGraph):
    def __init__(self, base_crystal):
        self.base_crystal = base_crystal

    @property
    def crystal_weight(self):
        return tuple(-w for w in self.base_crystal.crystal_weight)

    def raising_operator(self, index):
        lowered = self.base_crystal.lowering_operator(index)
        if lowered is None:
            return None
        return CrystalGraphDual(lowered)

    def lowering_operator(self, index):
        raised = self.base_crystal.raising_operator(index)
        if raised is None:
            return None
        return CrystalGraphDual(raised)

    def crystal_length(self):
        return self.base_crystal.crystal_length()

class CrystalGraphReverse(CrystalGraph):
    def __init__(self, base_crystal):
        self.base_crystal = base_crystal

    @property
    def crystal_weight(self):
        return self.base_crystal.crystal_weight[::-1]

    def raising_operator(self, index):
        n = self.crystal_length()
        raised = self.base_crystal.raising_operator(n + 1 - index)
        if raised is None:
            return None
        return CrystalGraphReverse(raised)

    def lowering_operator(self, index):
        n = self.crystal_length()
        lowered = self.base_crystal.lowering_operator(n + 1 - index)
        if lowered is None:
            return None
        return CrystalGraphReverse(lowered)

    def crystal_length(self):
        return self.base_crystal.crystal_length()


# There is a decomposition here into subcrystals
# NOT COMMUTATIVE TENSOR PRODUCT
class CrystalGraphTensor(CrystalGraph):

    @property
    def crystal_weight(self):
        result = self.factors[0].crystal_weight
        for factor in self.factors[1:]:
            result = tuple(a + b for a, b in zip_longest(result, factor.crystal_weight, fillvalue=0))
        return result

    def __hash__(self):
        return hash(self.factors)

    def __getitem__(self, i):
        return self.factors[i]

    def __len__(self):
        return len(self.factors)

    def __iter__(self):
        return iter(self.factors)

    def __add__(self, other):
        if isinstance(other, CrystalGraphTensor):
            return CrystalGraphTensor(*self.factors, *other.factors)
        if isinstance(other, tuple):
            return CrystalGraphTensor(*self.factors, *other)
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, tuple):
            return CrystalGraphTensor(*other, *self.factors)
        return NotImplemented

    def __eq__(self, other):
        return type(self) is type(other) and self.factors == other.factors

    @property
    def args(self):
        return self.factors

    def weight_bump(self):
        return type(self)(*(f.weight_bump() for f in self.factors))

    def all_highest_weights(self):
        import itertools
        full_crystals = [f.full_crystal for f in self.factors]
        highest_weights = set()
        for elements in itertools.product(*full_crystals):
            highest_weights.add(CrystalGraphTensor(*elements).to_highest_weight()[0])
        return highest_weights

    def _sympystr(self, printer):
        return printer.stringify(self.factors, sep=" # ")

    def _pretty(self, printer):
        return printer._print_TensorProduct(self)

    def _get_args_for_traditional_printer(self):
        return None, self.factors

    def __init__(self, *factors):
        self.factors = factors

    def crystal_length(self):
        return max(factor.crystal_length() for factor in self.factors)

    def _left_folded_ep_phi(self, index):
        """Compute (epsilon, phi) for left-folded sub-tensors T_k = (...((b_0 ⊗ b_1) ⊗ b_2) ⊗ ... ⊗ b_k)."""
        eps = self.factors[0].epsilon(index)
        phi = self.factors[0].phi(index)
        result = [(eps, phi)]
        for k in range(1, len(self.factors)):
            eps_k = self.factors[k].epsilon(index)
            phi_k = self.factors[k].phi(index)
            eps, phi = eps + max(0, eps_k - phi), phi_k + max(0, phi - eps_k)
            result.append((eps, phi))
        return result

    def lowering_operator(self, index):
        n = len(self.factors)
        ep = self._left_folded_ep_phi(index)
        for k in range(n - 1, 0, -1):
            if self.factors[k].epsilon(index) < ep[k - 1][1]:
                continue
            result = self.factors[k].lowering_operator(index)
            if result is None:
                return None
            new_factors = list(self.factors)
            new_factors[k] = result
            return CrystalGraphTensor(*new_factors)
        result = self.factors[0].lowering_operator(index)
        if result is None:
            return None
        new_factors = list(self.factors)
        new_factors[0] = result
        return CrystalGraphTensor(*new_factors)

    def raising_operator(self, index):
        n = len(self.factors)
        ep = self._left_folded_ep_phi(index)
        for k in range(n - 1, 0, -1):
            if self.factors[k].epsilon(index) > ep[k - 1][1]:
                result = self.factors[k].raising_operator(index)
                if result is None:
                    return None
                new_factors = list(self.factors)
                new_factors[k] = result
                return CrystalGraphTensor(*new_factors)
        result = self.factors[0].raising_operator(index)
        if result is None:
            return None
        new_factors = list(self.factors)
        new_factors[0] = result
        return CrystalGraphTensor(*new_factors)

    def epsilon(self, i):
        return self._left_folded_ep_phi(i)[-1][0]

    def phi(self, i):
        return self._left_folded_ep_phi(i)[-1][1]
