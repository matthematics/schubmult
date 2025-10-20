from functools import cache

from .rc_graph import RCGraph


class RCGraphCrystalChute(RCGraph):

    def __hash__(self):
        return hash((RCGraphCrystalChute, tuple(self)))
    
    def __eq__(self, other):
        return type(self) is type(other) and tuple(self) == tuple(other)


    def __new__(cls, rows):
        return RCGraphCrystalChute.__xnew_cached__(cls, tuple(rows))

    @staticmethod
    @cache
    def __xnew_cached__(_class, rows):
        return RCGraph.__xnew__(_class, rows)

    @staticmethod
    def __xnew__(_class, rows):
        return RCGraph.__new__(_class, rows)

    def pairing(self, row):
        pairs = {}
        for index in range(len(self[row - 1])):
            for index2 in range(len(self[row]) - 1, -1, -1):
                if index2 in pairs.values():
                    continue
                if self[row - 1][index] - row <= self[row][index2] - row - 1:
                    pairs[index] = index2
                    break
        return pairs

    def lowering_operator(self, row):
        if row < 1 or row >= len(self):
            return None
        pairs = self.pairing(row)
        if len(pairs) == len(self[row - 1]):
            return None
        leftmost_unpaired = max(set(range(len(self[row - 1]))) - set(pairs.keys()))
        col_unpaired = self[row - 1][leftmost_unpaired] - row + 1
        m = 1
        while m < col_unpaired and self.has_element(row, col_unpaired - m) and self.has_element(row + 1, col_unpaired - m):
            m += 1
        new_row_i = [a for a in self[row - 1] if a != self[row - 1][leftmost_unpaired]]
        new_row_ip1 = sorted([*self[row], col_unpaired - m + row], reverse=True)
        rc = type(self)([*self[:row], tuple(new_row_i), tuple(new_row_ip1), *self[row + 2:]])
        return rc

    def raising_operator(self, row):
        if row < 1 or row >= len(self):
            return None

        pairs0 = self.pairing(row)
        pairs = {}
        for k, v in pairs0.items():
            pairs[v] = k
        if len(pairs) == len(self[row]):
            return None
        rightmost_unpaired = min(set(range(len(self[row]))) - set(pairs.keys()))
        col_unpaired = self[row][rightmost_unpaired] - row
        n = col_unpaired + 1
        while self.has_element(row + 1, n):
            n += 1
        new_row_i = sorted([*self[row - 1], n + row - 1], reverse=True)
        new_row_ip1 = [a for a in self[row] if a != self[row][rightmost_unpaired]]
        rc = type(self)([*self[:row], tuple(new_row_i), tuple(new_row_ip1), *self[row + 2:]])
        return rc
