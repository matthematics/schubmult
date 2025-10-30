# filepath: /home/matthematics/schubmult/src/schubmult/rings/plactic.py
from functools import cache

from schubmult.schub_lib.perm_lib import Permutation

from ..utils._grid_print import GridPrint
from .crystal_graph import CrystalGraph


class Plactic(GridPrint, CrystalGraph):


    def __iter__(self):
        for row in self._word:
            yield from row

    # in order of row row
    @property
    def iter_boxes(self):
        for i in range(len(self._word) - 1, -1, -1):
            for j in range(len(self._word[i])):
                yield (i, j)

    def up_jdt_slide(self, row, col):
        """
        Perform an upward jeu de taquin slide starting from the given (row, col)
        position (0-indexed). Returns a new Plactic tableau.
        """
        new_word = [list(r) for r in self._word]
        i, j = row, col
        if self[i, j] != 0:
            raise ValueError(f"up_jdt_slide starting position must be empty, got {self[i, j]=}")
        if self[i - 1, j] == 0 and self[i, j - 1] == 0:
            raise ValueError("up_jdt_slide starting position has no valid moves")
        if i == len(new_word):
            new_word.append([0] * (j + 1))
        if j == len(new_word[i]):
            new_word[i].append(0)
        while True:
            up = new_word[i - 1][j] if i - 1 >= 0 and j < len(new_word[i - 1]) else None
            left = new_word[i][j - 1] if j - 1 >= 0 else None

            if up is None and left is None:
                # No more moves possible
                break
            if up is None:
                # Can only move left
                new_word[i][j] = left
                j -= 1
            elif left is None:
                # Can only move up
                new_word[i][j] = up
                i -= 1
            else:
                # Both moves possible; choose the larger entry
                if up >= left:
                    new_word[i][j] = up
                    i -= 1
                else:
                    new_word[i][j] = left
                    j -= 1

        new_word[i][j] = 0

        return Plactic(tuple(tuple(r) for r in new_word))

    def down_jdt_slide(self, row, col):
        """
        Perform a jeu de taquin slide starting from the given (row, col)
        position (0-indexed). Returns a new Plactic tableau.
        """
        new_word = [list(r) for r in self._word]
        i, j = row, col
        if self[i, j] != 0:
            raise ValueError(f"down_jdt_slide starting position must be empty, got {self[i, j]=}")
        if self[i + 1, j] == 0 and self[i, j + 1] == 0:
            raise ValueError("down_jdt_slide starting position has no valid moves")
        while True:
            down = new_word[i + 1][j] if i + 1 < len(new_word) and j < len(new_word[i + 1]) else None
            right = new_word[i][j + 1] if j + 1 < len(new_word[i]) else None

            if down is None and right is None:
                # No more moves possible
                break
            if down is None:
                # Can only move right
                new_word[i][j] = right
                j += 1
            elif right is None:
                # Can only move down
                new_word[i][j] = down
                i += 1
            else:
                # Both moves possible; choose the smaller entry
                if down < right:
                    new_word[i][j] = down
                    i += 1
                else:
                    new_word[i][j] = right
                    j += 1

        # Remove the last empty cell
        new_word[i].pop()
        # Remove any empty rows at the bottom
        while new_word and len(new_word[-1]) == 0:
            new_word.pop()

        return Plactic(tuple(tuple(r) for r in new_word))

    def __init__(self, word=()):
        # ensure we always store a flat tuple of rows (each row a tuple of ints)
        self._word = tuple(tuple(r) for r in word)

    def shiftup(self, k):
        """Return a new Plactic with all entries increased by k."""
        new_word = tuple(tuple(int(a) + k for a in row) for row in self._word)
        return Plactic(new_word)

    @classmethod
    def all_ss_tableaux(cls, shape, max_entry):
        """Generate all semistandard tableaux of given shape with entries <= max_entry."""
        tableaux = set()
        current_tableau = [[0 for _ in range(row_len)] for row_len in shape]

        def _recurse(row, col):
            if row == len(shape):
                # Base case: A complete tableau is found
                tableaux.add(cls([tuple(r) for r in current_tableau]))
                return

            if col == shape[row]:
                # Move to the next row if current row is filled
                _recurse(row + 1, 0)
                return

            # Determine the minimum possible value for the current cell
            min_val = 1
            if col > 0:
                min_val = max(min_val, current_tableau[row][col - 1])
            if row > 0:
                min_val = max(min_val, current_tableau[row - 1][col] + 1)

            for val in range(min_val, max_entry + 1):
                current_tableau[row][col] = val
                _recurse(row, col + 1)
                # Backtrack (no explicit removal needed as it's overwritten in next iteration)

        _recurse(0, 0)
        return tableaux

    @property
    def row_word(self):
        """Return the row-reading word as a flat tuple."""
        word = tuple(a for row in reversed(self._word) for a in row if a != 0)
        return word

    @property
    def column_word(self):
        wrd = []
        for i in range(self.cols):
            for j in range(self.rows - 1, -1, -1):
                if i < len(self._word[j]):
                    wrd.append(self._word[j][i])
        return tuple(wrd)

    def transpose(self):
        """Return the transpose of this Plactic tableau."""
        new_word = []
        for i in range(self.cols):
            new_row = []
            for j in range(self.rows):
                if i < len(self._word[j]):
                    new_row.append(self._word[j][i])
            new_word.append(tuple(new_row))
        return Plactic(tuple(new_word))

    @property
    def rows(self):
        return len(self._word)

    @property
    def cols(self):
        if len(self._word) == 0:
            return 0
        return max(len(r) for r in self._word)

    def __getitem__(self, key):
        # FLIPPED FOR PRINTING
        if isinstance(key, int):
            return self._word[key]
        if isinstance(key, tuple):
            i, j = key
            if len(self._word) == 0 or i >= len(self._word) or j >= len(self._word[i]):
                return None
            return self._word[i][j]
        is_slice = isinstance(key, slice)

        if is_slice:
            return tuple(self._word[n] for n in range(len(self))[key])

        raise ValueError(f"Bad indexing {key=}")

    @staticmethod
    def _remap_word(word, fn):
        """Apply fn to every entry of a tableau-like tuple-of-rows."""
        return tuple(tuple(fn(x) for x in row) for row in tuple(word))

    def invert(self):
        """
        Return a Plactic whose entries are remapped so that standard
        (increasing) insertion order applies. If reverse_semistandard is True
        we negate entries (so larger original becomes smaller).
        """
        return Plactic(Plactic._remap_word(self._word, lambda x: -int(x)))

    def __mul__(self, other):
        """
        Plactic product: insert entries of `other` in row-reading order
        (top-to-bottom, left-to-right) into a copy of self.
        """
        if not isinstance(other, Plactic):
            return NotImplemented
        word = [*self.row_word, *other.row_word]
        pl = Plactic()
        return pl.rs_insert(*word)

    @property
    def shape(self):
        return tuple(len(row) for row in self._word)

    @classmethod
    def from_word(cls, word):
        # accept any iterable-of-rows and normalize to tuple-of-tuples
        return cls(tuple(tuple(r) for r in word))

    def __hash__(self):
        return hash(self._word)

    @staticmethod
    def _rs_insert(word, letter, i=0):
        """
        Row-insertion: insert `letter` into `word` (tuple-of-tuples) starting at row i.
        Returns a new tuple-of-tuples representing the tableau after insertion.

        Algorithm:
        - If row i doesn't exist, append a new row with the letter.
        - Otherwise find the first entry in row i that is > letter, replace it (bump)
          and recursively insert the bumped value into row i+1.
        - If no entry > letter, append letter to the end of row i.
        """
        # normalize input
        word = tuple(tuple(r) for r in word)

        if i >= len(word):
            # append new row containing the letter
            return (*word, (int(letter),))

        row = list(word[i])
        x = int(letter)

        # find first entry strictly greater than x to bump
        for k, val in enumerate(row):
            if val > x:
                bumped = val
                row[k] = x
                new_word = (*word[:i], tuple(row), *word[i + 1 :])
                # recursively insert bumped into next row
                return Plactic._rs_insert(new_word, bumped, i=i + 1)

        # nothing to bump -> append to end of this row
        row.append(x)
        return (*word[:i], tuple(row), *word[i + 1 :])

    def rs_insert(self, *letters):
        """Insert one or more letters in sequence (row-insertion) and return a new Plactic."""
        # accept either rs_insert(a,b,...) or rs_insert([a,b,...])
        if len(letters) == 1 and isinstance(letters[0], (list, tuple)):
            seq = list(letters[0])
        else:
            seq = list(letters)

        new_word = self._word
        for letter in seq:
            new_word = Plactic._rs_insert(new_word, int(letter))
        return Plactic(new_word)

    def __eq__(self, other):
        if isinstance(other, Plactic):
            return self._word == other._word
        return False

    # -- CrystalGraph API wrappers (delegate to the NilPlactic <-> RCGraph machinery) --
    # These adapt the RCGraph crystal operators/queries to Plactic tableaux by
    # converting the tableau to its highest-weight RC graph (via NilPlactic.hw_rc),
    # delegating to the RCGraph implementation, then mapping results back to Plactic.

    def raising_operator(self, i):
        """Crystal raising operator e_i on the Plactic tableau (delegates to RCGraph)."""
        word = [*self.row_word]
        opening_stack = []
        closing_stack = []
        for index in range(len(word)):
            if word[index] == i + 1:
                opening_stack.append(index)
            elif word[index] == i:
                if len(opening_stack) > 0:
                    opening_stack.pop()
                else:
                    closing_stack.append(index)
        if len(opening_stack) == 0:
            return None
        index_to_change = opening_stack[0]
        word[index_to_change] = i
        result = Plactic().rs_insert(*word)
        assert result.shape == self.shape, f"{result.shape=} {self.shape=}"
        return result

    def lowering_operator(self, i):
        """Crystal lowering operator f_i on the Plactic tableau (delegates to RCGraph)."""
        word = [*self.row_word]
        opening_stack = []
        closing_stack = []
        for index in range(len(word)):
            if word[index] == i + 1:
                opening_stack.append(index)
            elif word[index] == i:
                if len(opening_stack) > 0:
                    opening_stack.pop()
                else:
                    closing_stack.append(index)
        if len(closing_stack) == 0:
            return None
        index_to_change = closing_stack[-1]
        word[index_to_change] = i + 1
        result = Plactic().rs_insert(*word)
        assert result.shape == self.shape, f"{result.shape=} {self.shape=} {result=} {self=}"
        return result
    # def _bracket_positions(self, i: int):
    #     """
    #     Bracketing scan for indices i (as '+') and i+1 (as '-') on the tableau's
    #     row-reading word. Returns (plus_stack, minus_unpaired, word_list) where:
    #       - plus_stack is a list of indices of unmatched i (these are the '+' positions)
    #       - minus_unpaired is a list of indices of unmatched i+1 (these are the '-' positions)
    #       - word_list is a mutable list copy of the row-reading word used for modification.
    #     The scan follows the standard left->right cancellation rule.
    #     """
    #     word = list(self.row_word)
    #     plus_stack = []
    #     minus_unpaired = []
    #     for idx, a in enumerate(word):
    #         if a == i:
    #             plus_stack.append(idx)
    #         elif a == i + 1:
    #             if plus_stack:
    #                 # cancel with the most recent unmatched '+'
    #                 plus_stack.pop()
    #             else:
    #                 # unmatched '-'
    #                 minus_unpaired.append(idx)
    #     return plus_stack, minus_unpaired, word

    # def raising_operator(self, i):
    #     """
    #     Crystal raising operator e_i: change the leftmost unmatched (i+1) -> i,
    #     if any; rebuild the Plactic from the modified row-reading word and return it.
    #     Returns None if epsilon_i == 0.
    #     """
    #     plus, minus, word = self._bracket_positions(i)
    #     if len(minus) == 0:
    #         return None
    #     index_to_change = minus[0]  # leftmost unmatched i+1
    #     word[index_to_change] = i
    #     result = Plactic().rs_insert(*word)
    #     # shape-preserving check (crystal operators preserve shape)
    #     assert result.shape == self.shape, f"shape changed by e_{i}: {self.shape} -> {result.shape}"
    #     return result

    # def lowering_operator(self, i):
    #     """
    #     Crystal lowering operator f_i: change the rightmost unmatched i -> i+1,
    #     if any; rebuild the Plactic from the modified row-reading word and return it.
    #     Returns None if phi_i == 0.
    #     """
    #     plus, minus, word = self._bracket_positions(i)
    #     if len(plus) == 0:
    #         return None
    #     index_to_change = plus[-1]  # rightmost unmatched i
    #     word[index_to_change] = i + 1
    #     result = Plactic().rs_insert(*word)
    #     assert result.shape == self.shape, f"shape changed by f_{i}: {self.shape} -> {result.shape}"
    #     return result

    @property
    def crystal_weight(self):
        """Return the crystal weight of this tableau (delegated to RCGraph)."""

    def crystal_length(self):
        """Return the length/number of rows used for the crystal"""

    @classmethod
    def yamanouchi(cls, shape):
        """
        Return the Yamanouchi (highest-weight) tableau of the given shape.
        """
        new_word = []
        for i in range(len(shape)):
            new_word.append([0] * shape[i])
            for j in range(shape[i]):
                new_word[i][j] = i + 1
        return cls(tuple(tuple(row) for row in new_word))

    def rectify(self):
        if len(self._word) == 0:
            return self
        if len(self._word[0]) == 0:
            return self.__class__([self._word[0], *self.__class__(self.word[1:]).rectify()._word])
        if self._word[0][0] != 0:
            return self
        index = 0
        while index < len(self._word[0]) and self._word[0][index] == 0:
            index += 1
        if index == len(self._word[0]):
            if len(self._word) == 1:
                return self
            if self._word[1][0] == 0:
                return self.__class__([self._word[0], *self.__class__(self._word[1:]).rectify()._word])
            return self.down_jdt_slide(0, 0).rectify()
        index -= 1
        return self.down_jdt_slide(0, index).rectify()

    @classmethod
    def superstandard(cls, shape):
        if shape is None:
            return None
        new_word = []
        index = 1
        for i in range(len(shape)):
            new_word.append([0] * shape[i])
            for j in range(shape[i]):
                new_word[i][j] = index
                index += 1
        return cls(tuple(tuple(row) for row in new_word))

    @property
    def is_semistandard(self):
        for i in range(len(self._word)):
            for j in range(len(self._word[i])):
                if j > 0 and self._word[i][j] < self._word[i][j - 1]:
                    return False
                if i > 0 and j < len(self._word[i - 1]) and self._word[i][j] <= self._word[i - 1][j]:
                    return False
        return True

    def reverse_rsk(self, recording_tableau):
        """
        Inverse RSK (row-insertion) for the pair (P,Q) where `self` is P and
        `recording_tableau` is the standard recording tableau Q of the same shape.

        Returns the original word as a list of integers (in insertion order).
        """
        # mutable copies of P and Q rows
        P = [list(r) for r in self._word]
        Q = [list(r) for r in recording_tableau._word]

        total = sum(len(r) for r in Q)
        word_rev = []

        for t in range(total, 0, -1):
            # find position (i,j) of t in Q
            pos_i = pos_j = None
            for i, row in enumerate(Q):
                for j, val in enumerate(row):
                    if val == t:
                        pos_i, pos_j = i, j
                        break
                if pos_i is not None:
                    break
            if pos_i is None:
                raise ValueError(f"recording tableau missing entry {t}")

            i, j = pos_i, pos_j

            # remove the box from Q
            Q[i].pop(j)
            if len(Q[i]) == 0:
                Q.pop(i)

            # remove the corresponding entry from P and start reverse bumping
            x = P[i].pop(j)
            if len(P[i]) == 0:
                P.pop(i)

            # move upwards: in row r find rightmost entry < x, replace it and continue,
            # if none found in a row, stop and output x
            for r in range(i - 1, -1, -1):
                replaced = False
                for k in range(len(P[r]) - 1, -1, -1):
                    if P[r][k] < x:
                        y = P[r][k]
                        P[r][k] = x
                        x = y
                        replaced = True
                        break
                if not replaced:
                    break

            word_rev.append(int(x))

        # reverse collected letters to obtain original insertion order
        return list(reversed(word_rev))

    @classmethod
    def rsk_insert(cls, *letters):
        """
        Perform ordinary RSK (row insertion) on the given sequence of letters,
        starting from this Plactic as the initial P-tableau. Returns a pair
        (P_tableau, Q_tableau) where both are Plactic instances and Q is the
        standard recording tableau with entries 1..m (in insertion order).

        Usage:
          P, Q = Plactic().rsk_insert(3,1,2,1)
          or
          P, Q = Plactic().rsk_insert([3,1,2,1])
        """
        # normalize letters accept either rsk_insert(a,b,...) or rsk_insert(iterable)
        if len(letters) == 1 and isinstance(letters[0], (list, tuple)):
            seq = list(letters[0])
        else:
            seq = list(letters)

        # mutable rows for P and Q
        P_rows = []
        Q_rows = []

        def _insert_letter(rows, x):
            """Insert x into mutable rows (list of list) and return final (i,j)."""
            i = 0
            while True:
                if i >= len(rows):
                    rows.append([])
                row = rows[i]
                if len(row) == 0 or x >= max(row):
                    row.append(x)
                    return i, len(row) - 1
                # bump smallest entry > x
                for idx, val in enumerate(row):
                    if val > x:
                        bumped = val
                        row[idx] = x
                        x = bumped
                        break
                i += 1

        # perform insertions, updating Q at appended cell each time
        for k, letter in enumerate(seq, start=1):
            i, j = _insert_letter(P_rows, int(letter))
            # ensure Q has row i
            while i >= len(Q_rows):
                Q_rows.append([])
            # append placeholder zeros for earlier existing length mismatches
            # (Q_rows[i] should be same length as P_rows[i]-1 before append)
            # append the insertion index to the new/extended cell
            if j == len(Q_rows[i]):
                Q_rows[i].append(k)
            else:
                # In typical RSK we always append a new box at the final row.
                # However guard against unexpected states by inserting at position j.
                Q_rows[i].insert(j, k)
                # ensure consistency of row lengths in Q with P:
                # if this caused earlier rows to mismatch, leave as-is (shouldn't happen)
        # convert to Plactic objects
        P = cls(tuple(tuple(r) for r in P_rows))
        Q = cls(tuple(tuple(r) for r in Q_rows))
        return P, Q


