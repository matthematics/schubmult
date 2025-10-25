from functools import cache

from sympy.printing.defaults import Printable


class CrystalGraph(Printable):
    def raising_operator(self, index):
        """The raising operator for the crystal graph."""
        raise NotImplementedError

    def lowering_operator(self, index):
        """The lowering operator for the crystal graph."""
        raise NotImplementedError


    def phi(self, i):
        if i == 0:
            return 0
        if hasattr(self, "_cached_phi") and i in self._cached_phi:
            return self._cached_phi[i]
        rc = self
        cnt = 0
        while rc is not None:
            rc = rc.lowering_operator(i)
            if rc is not None:
                cnt += 1
        if hasattr(self, "_cached_phi"):
            self._cached_phi[i] = cnt
        else:
            self._cached_phi = {i: cnt}
        return cnt


    def epsilon(self, i):
        if i == 0:
            return 0
        if hasattr(self, "_cached_epsilon") and i in self._cached_epsilon:
            return self._cached_epsilon[i]

        rc = self
        cnt = 0
        while rc is not None:
            rc = rc.raising_operator(i)
            if rc is not None:
                cnt += 1
        if hasattr(self, "_cached_epsilon"):
            self._cached_epsilon[i] = cnt
        else:
            self._cached_epsilon = {i: cnt}
        return cnt

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
    def full_crsytal(self):
        hw, _ = self.to_highest_weight()
        crystal = set()
        stack = [hw]
        while len(stack) > 0:
            elem = stack.pop()
            crystal.add(elem)
            for i in range(1, elem.crystal_length()):
                new_elem = elem.lowering_operator(i)
                if new_elem is not None:
                    stack.append(new_elem)
        return crystal


# There is a decomposition here into subcrystals
class CrystalGraphTensor(CrystalGraph):
    @property
    def args(self):
        return self.factors

    def all_highest_weights(self):
        right_element, _ = self.factors[1].to_highest_weight()
        full_crystal = self.factors[0].full_crystal
        highest_weights = set()
        for left_element in full_crystal:
            good = True
            for i in range(1, self.crystal_length()):
                if left_element.epsilon(i) > right_element.phi(i):
                    good = False
                    break
            if good:
                highest_weights.add(CrystalGraphTensor(left_element, right_element))
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

    def __eq__(self, other):
        return type(self) is type(other) and self.factors == other.factors

    def crystal_length(self):
        return max(factor.crystal_length() for factor in self.factors)

    def phi(self, index):
        return self.factors[0].phi(index) + max(0, self.factors[1].phi(index) - self.factors[0].epsilon(index))

    def epsilon(self, index):
        return self.factors[1].epsilon(index) + max(0, self.factors[0].epsilon(index) - self.factors[1].phi(index))

    def lowering_operator(self, index):
        if self.factors[0].epsilon(index) < self.factors[1].phi(index):
            tz = CrystalGraphTensor(self.factors[0], self.factors[1].lowering_operator(index))
            if tz.factors[1] is None:
                return None
            return tz
        tz = CrystalGraphTensor(self.factors[0].lowering_operator(index), self.factors[1])
        if tz.factors[0] is None:
            return None
        return tz

    def raising_operator(self, index):
        if self.factors[0].epsilon(index) > self.factors[1].phi(index):
            tz = CrystalGraphTensor(self.factors[0].raising_operator(index), self.factors[1])
            if tz.factors[0] is None:
                return None
            return tz
        tz = CrystalGraphTensor(self.factors[0], self.factors[1].raising_operator(index))
        if tz.factors[1] is None:
            return None
        return tz
