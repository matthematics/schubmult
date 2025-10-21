

class CrystalGraph:
    def raising_operator(self, index):
        """The raising operator for the crystal graph."""
        raise NotImplementedError

    def lowering_operator(self, index):
        """The lowering operator for the crystal graph."""
        raise NotImplementedError

    def phi(self, index):
        """The phi function for the crystal graph."""
        raise NotImplementedError

    def epsilon(self, index):
        """The epsilon function for the crystal graph."""
        raise NotImplementedError

    def to_lowest_weight(self):
        """Return the lowest weight element in the connected component."""
        g = self
        lower_seq = []
        found = True
        while found:
            found = False
            for row in range(1, g.crystal_length()):
                g0 = g.lowering_operator(row)
                if g0 is not None:
                    found = True
                    g = g0
                    lower_seq.append(row)
                    break
        return (g, tuple(lower_seq))

    def to_highest_weight(self):
        """Return the highest weight element in the connected component."""
        g = self
        raise_seq = []
        found = True
        while found:
            found = False
            for row in range(1, g.crystal_length()):
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
                return None
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



class CrystalGraphTensor(CrystalGraph):
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
