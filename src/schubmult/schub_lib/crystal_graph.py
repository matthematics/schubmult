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
        while found:
            found = False
            for row in range(1, length):
                g0 = g.lowering_operator(row)
                if g0 is not None:
                    found = True
                    g = g0
                    lower_seq.append(row)
                    break
        return (g, tuple(lower_seq))

    def to_highest_weight(self, length=None):
        """Return the highest weight element in the connected component."""
        g = self
        raise_seq = []
        found = True
        if length is None:
            length = self.crystal_length()
        while found:
            found = False
            for row in range(1, length):
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
        elif f > e:
            for _ in range(f - e):
                new_rc = new_rc.lowering_operator(index)
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

    def truncated_crystal(self, length):
        hw, _ = self.to_highest_weight(length=length)
        crystal = set()
        stack = [hw]
        while len(stack) > 0:
            elem = stack.pop()
            crystal.add(elem)
            for i in range(1, length):
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
        return tuple(a + b for a, b in zip_longest(self.factors[0].crystal_weight, self.factors[1].crystal_weight, fillvalue=0))

    def __hash__(self):
        return hash(self.factors)

    def __eq__(self, other):
        return type(self) is type(other) and self.factors == other.factors



    @property
    def args(self):
        return self.factors

    def weight_bump(self):
        return type(self)(self.factors[0].weight_bump(), self.factors[1].weight_bump())

    def all_highest_weights(self):
        import itertools
        right_full_crystal = self.factors[1].full_crystal
        full_crystal = self.factors[0].full_crystal
        highest_weights = set()
        for left_element, right_element in itertools.product(full_crystal, right_full_crystal):
            highest_weights.add(CrystalGraphTensor(left_element, right_element).to_highest_weight()[0])
        return highest_weights

    def _sympystr(self, printer):
        return printer.stringify(self.factors, sep=" # ")

    def _pretty(self, printer):
        return printer._print_TensorProduct(self)

    def _get_args_for_traditional_printer(self):
        return None, self.factors

    def __init__(self, *factors):
        if len(factors) == 1:
            self.factors = factors[0]
        else:
            self.factors = factors

    def crystal_length(self):
        return max(factor.crystal_length() for factor in self.factors)

    def lowering_operator(self, index):
        if self.factors[1].epsilon(index) < self.factors[0].phi(index):
            tz = CrystalGraphTensor(self.factors[0].lowering_operator(index), self.factors[1])
            if tz.factors[0] is None:
                return None
            return tz
        tz = CrystalGraphTensor(self.factors[0], self.factors[1].lowering_operator(index))
        if tz.factors[1] is None:
            return None
        return tz

    def raising_operator(self, index):
        if self.factors[1].epsilon(index) > self.factors[0].phi(index):
            tz = CrystalGraphTensor(self.factors[0], self.factors[1].raising_operator(index))
            if tz.factors[1] is None:
                return None
            return tz
        tz = CrystalGraphTensor(self.factors[0].raising_operator(index), self.factors[1])
        if tz.factors[0] is None:
            return None
        return tz

    def epsilon(self, i):
        return max(self.factors[0].epsilon(i),self.factors[1].epsilon(i) - self.factors[0].crystal_weight[i - 1] + self.factors[0].crystal_weight[i])

    def phi(self, i):
        return max(self.factors[1].phi(i), self.factors[0].phi(i) + self.factors[1].crystal_weight[i - 1] - self.factors[1].crystal_weight[i])

    # def phi(self, i):
    #     return self.factors[0].phi(i) + max(0, self.factors[1].phi(i) - self.factors[0].epsilon(i))

    # def epsilon(self, i):
    #     return self.factors[1].epsilon(i) + max(0, self.factors[0].epsilon(i) - self.factors[1].phi(i))
