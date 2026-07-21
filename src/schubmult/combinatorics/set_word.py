from .crystal_graph import CrystalGraph, CrystalGraphTensor


class SetLetter(CrystalGraph, frozenset):
    """
    A set of ints with a sqrt(gl_n) crystal structure.
    """

    def __new__(cls, iterable, length=None):
        obj = frozenset.__new__(cls, iterable)
        length = 100
        if length is not None:
            obj._length = length
        else:
            obj._length = max(iterable) if iterable else 0
        if obj._length > 0:
            weight = [0] * obj._length
            for i in obj:
                weight[i - 1] += 1
            obj._weight = tuple(weight)
        else:
            obj._weight = ()
        return obj

    def __init__(self, iterable, length=None):
        pass

    def __repr__(self):
        return f"SetLetter({super().__repr__()})"

    def crystal_length(self):
        return self._length

    @property
    def crystal_weight(self):
        return self._weight

    def raising_operator(self, i):
        if i >= self._length or i < 1:
            return None
        if i + 1 in self and i not in self:
            return SetLetter(set(self) | {i}, length=self._length)
        if i in self and i + 1 in self:
            return SetLetter(set(self) - {i + 1}, length=self._length)
        return None

    def lowering_operator(self, i):
        if i >= self._length or i < 1:
            return None
        if i in self and i + 1 not in self:
            return SetLetter(set(self) | {i + 1}, length=self._length)
        if i in self and i + 1 in self:
            return SetLetter(set(self) - {i}, length=self._length)
        return None

class SetWord(CrystalGraphTensor):
    """
    A tuple of SetLetters with a sqrt(gl_n) crystal structure.
    """

    def __init__(self, *args):
        super().__init__(*args)

    def __repr__(self):
        return f"SetWord({self.factors.__repr__()})"

    def to_wc_graph(self, rows):
        from .wc_graph import WCGraph
        wc = WCGraph([()]).resize(rows)
        for j in range(1, len(self.factors) + 1):
            for i in self.factors[j - 1]:
                wc = wc.toggle_ref_at(i, j)
        return wc

    @classmethod
    def from_wc_graph(cls, wc):
        columns = wc.cols
        rows = wc.perm.max_descent
        factors = []
        for j in range(1, columns + 1):
            letter = set()
            for i in range(1, rows + 1):
                if wc.has_element(i, j):
                    letter.add(i)
            factors.append(SetLetter(letter, length=rows))
        return cls(*factors)
