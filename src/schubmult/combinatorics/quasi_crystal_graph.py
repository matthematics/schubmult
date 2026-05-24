from itertools import zip_longest

from .crystal_graph import CrystalGraph, CrystalGraphTensor


class QuasiCrystalGraph(CrystalGraph):
    def quasi_raising_operator(self, index):
        """The raising operator for the crystal graph."""
        raise NotImplementedError

    def quasi_lowering_operator(self, index):
        """The lowering operator for the crystal graph."""
        raise NotImplementedError

    def to_quasi_lowest_weight(self):
        """Return the lowest weight element in the same quasi-crystal component."""
        current = self
        found = True
        lower_seq = []
        while found:
            found = False
            for i in range(1, self.crystal_length()):
                next_element = current.quasi_lowering_operator(i)
                if next_element is None:
                    continue
                current = next_element
                lower_seq.append(i)
                found = True
                break
        return current, tuple(lower_seq)

    def reverse_quasi_lower_seq(self, seq):
        """Apply the reverse of the given lowering sequence."""
        current = self
        for i in reversed(seq):
            next_element = current.quasi_raising_operator(i)
            if next_element is None:
                raise ValueError(f"Invalid raising sequence {seq} for element {self}")
            current = next_element
        return current

    @classmethod
    def _wrap(cls, element):
        if isinstance(element, cls):
            return element
        wrapper = cls()
        wrapper.quasi_lowering_operator = element.quasi_lowering_operator
        wrapper.quasi_raising_operator = element.quasi_raising_operator
        wrapper.epsilon = element.epsilon
        wrapper.phi = element.phi
        wrapper.crystal_weight = element.crystal_weight
        wrapper.crystal_length = element.crystal_length
        return wrapper


# There is a decomposition here into subcrystals
# NOT COMMUTATIVE TENSOR PRODUCT
class QuasiCrystalGraphTensor(QuasiCrystalGraph):
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
        if isinstance(other, QuasiCrystalGraphTensor):
            return QuasiCrystalGraphTensor(*self.factors, *other.factors)
        if isinstance(other, tuple):
            return QuasiCrystalGraphTensor(*self.factors, *other)
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, tuple):
            return QuasiCrystalGraphTensor(*other, *self.factors)
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

        full_quasi_crystals = [f.full_quasi_crystal for f in self.factors]
        highest_weights = set()
        for elements in itertools.product(*full_quasi_crystals):
            highest_weights.add(QuasiCrystalGraphTensor(*elements).to_highest_weight()[0])
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
        return max([factor.crystal_length() for factor in self.factors], default=0)

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

    def quasi_lowering_operator(self, index):
        n = len(self.factors)
        ep = self._left_folded_ep_phi(index)
        if len(ep) > 1 and ep[-2][0] > 0 and ep[-1][1] > 0:
            return None
        for k in range(n - 1, 0, -1):
            if self.factors[k].epsilon(index) < ep[k - 1][1]:
                continue
            result = self.factors[k].quasi_lowering_operator(index)
            if result is None:
                return None
            new_factors = list(self.factors)
            new_factors[k] = result
            return QuasiCrystalGraphTensor(*new_factors)
        result = self.factors[0].quasi_lowering_operator(index)
        if result is None:
            return None
        new_factors = list(self.factors)
        new_factors[0] = result
        return QuasiCrystalGraphTensor(*new_factors)

    def quasi_raising_operator(self, index):
        n = len(self.factors)
        ep = self._left_folded_ep_phi(index)
        if len(ep) > 1 and ep[-2][0] > 0 and ep[-1][1] > 0:
            return None
        for k in range(n - 1, 0, -1):
            if self.factors[k].epsilon(index) > ep[k - 1][1]:
                result = self.factors[k].quasi_raising_operator(index)
                if result is None:
                    return None
                new_factors = list(self.factors)
                new_factors[k] = result
                return QuasiCrystalGraphTensor(*new_factors)
        result = self.factors[0].quasi_raising_operator(index)
        if result is None:
            return None
        new_factors = list(self.factors)
        new_factors[0] = result
        return QuasiCrystalGraphTensor(*new_factors)

    def epsilon(self, i):
        return self._left_folded_ep_phi(i)[-1][0]

    def phi(self, i):
        return self._left_folded_ep_phi(i)[-1][1]
