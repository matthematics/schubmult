# filepath: /home/matthematics/schubmult/src/schubmult/rings/plactic.py
from schubmult.perm_lib import Permutation

from ._grid_print import GridPrint
from .crystal_graph import CrystalGraph


class Plactic(GridPrint, CrystalGraph):
    def __init__(self, word=()):
        # ensure we always store a flat tuple of rows (each row a tuple of ints)
        self._word = tuple(tuple(r) for r in word)

    def shiftup(self, k):
        """Return a new Plactic with all entries increased by k."""
        new_word = tuple(tuple(int(a) + k for a in row) for row in self._word)
        return Plactic(new_word)


    @property
    def row_word(self):
        """Return the row-reading word as a flat tuple."""
        word =  tuple(a for row in reversed(self._word) for a in row)
        assert len(word) == sum(len(r) for r in self._word), f"Length mismatch in row_word {len(word)=} vs {self._word=}"
        return word

    @property
    def column_word(self):
        wrd = []
        for i in range(self.cols):
            for j in range(self.rows - 1, -1, -1):
                if i < len(self._word[j]):
                    wrd.append(self._word[j][i])
        return tuple(wrd)

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
                return " "
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

    @staticmethod
    def from_word(word):
        # accept any iterable-of-rows and normalize to tuple-of-tuples
        return Plactic(tuple(tuple(r) for r in word))

    def __hash__(self):
        return hash(self._word)

    @staticmethod
    def _rs_insert(word, letter, i=0):
        """
        Row-standard insertion variant that returns a new word as a tuple of rows
        (each row itself a tuple of ints). `word` must be a tuple-of-tuples.
        """
        # Ensure word is tuple-of-tuples (idempotent if already correct)
        word = tuple(tuple(r) for r in word)

        if i == len(word):
            row_i = ()
        else:
            row_i = word[i]
        x0 = letter

        # If row empty or letter is >= max(row) then append to this row
        if len(row_i) == 0 or x0 >= max(row_i):
            new_row = (*row_i, x0)
            return (*word[:i], new_row, *word[i + 1 :])

        # otherwise bump the smallest entry > x0
        x1 = min(a for a in row_i if a > x0)
        new_first_row = list(row_i)
        idx = new_first_row.index(x1)
        new_first_row[idx] = x0
        new_word = (*word[:i], tuple(new_first_row), *word[i + 1 :])
        # recursively insert the bumped value into the next row
        return Plactic._rs_insert(new_word, x1, i=i + 1)

    def rs_insert(self, *letters):
        """Insert a letter/entry into this Plactic tableau and return a new Plactic."""
        new_word = self._word
        for letter in letters:
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
        assert result.shape == self.shape, f"{result.shape=} {self.shape=}"
        return result

    @property
    def crystal_weight(self):
        """Return the crystal weight of this tableau (delegated to RCGraph)."""


    def crystal_length(self):
        """Return the length/number of rows used for the crystal"""

    def yamanouchi(self):
        """
        Return the Yamanouchi (highest-weight) tableau of the same shape
        as this Plactic tableau.
        """
        new_word = []
        for i in range(len(self._word)):
            new_word.append([0] * len(self._word[i]))
            for j in range(len(self._word[i])):
                new_word[i][j] = i + 1
        return Plactic(tuple(tuple(row) for row in new_word))


class NilPlactic(Plactic):

    def hw_rc(self, length=None):
        """
        Return the highest-weight RCGraph corresponding to this NilPlactic
        tableau via Edelman–Greene insertion.
        """
        from schubmult.rings.rc_graph import RCGraph

        perm_word = list(reversed(self.row_word))
        graph = []
        last_letter = -1
        row = []
        for letter in perm_word:
            if last_letter != -1 and letter > last_letter:
                # need to add empty rows in between
                graph.append(tuple(row))
                row = []
            row.append(letter)
            last_letter = letter
        graph.append(tuple(row))
        graph = RCGraph(graph).normalize()
        if length < len(graph):
            raise ValueError(f"Requested length {length} too small for RCGraph of size {len(graph)}")
        if length > len(graph):
            graph = graph.resize(length)
        assert graph.perm == ~self.perm, f"{graph.perm=} {self.perm=}"
        assert graph.p_tableau == self, f"{graph.p_tableau=} {self=} {graph=}"
        return graph

    # 

    # def ring_product(self, other, length):
    #     """
    #     NilPlactic product: insert entries of `other` in Edelman–Greene
    #     row-reading order (top-to-bottom, left-to-right) into a copy of self.
    #     """
    #     if not isinstance(other, NilPlactic):
    #         return NotImplemented
    #     st = (self.hw_rc()).prod_with_rc(other.hw_rc())
    #     ret = {}
    #     for rc, coeff in st.items():
    #         ret[rc.p_tableau] = ret.get(rc.p_tableau, 0) + coeff
    #     return ret

    @property
    def perm(self):
        return Permutation.ref_product(*self.row_word)

    @staticmethod
    def _ed_insert(word, letter):
        """
        Row insertion for NilPlactic. `word` is a tuple-of-rows (each a tuple).
        Returns a tuple-of-rows.
        """
        word = tuple(tuple(r) for r in word)
        if len(word) == 0:
            return ((letter,),)

        row = tuple(word[0])
        index = 0
        while index < len(row) and row[index] < letter:
            index += 1

        # append to first row
        if index == len(row):
            return ((*row, letter), *tuple(word[1:]))

        # equal-case handling
        if row[index] == letter:
            if index < len(row) - 1 and row[index + 1] == letter + 1:
                # bump into next row
                return (row, *NilPlactic._ed_insert(tuple(word[1:]), letter + 1))
            # skip equal-run
            while index < len(row) and row[index] == letter:
                index += 1
            if index == len(row):
                return ((*row, letter), *tuple(word[1:]))

        # bump case
        new_row = list(row)
        bump = new_row[index]
        new_row[index] = letter
        return (tuple(new_row), *NilPlactic._ed_insert(tuple(word[1:]), bump))

    @staticmethod
    def ed_insert_rsk(word, word2, letter, letter2, i=0):
        """
        Edelman–Greene style two-row insertion.
        word and word2 are tuples-of-rows; always return pair of tuple-of-rows.
        """
        word = tuple(tuple(r) for r in word)
        word2 = tuple(tuple(r) for r in word2)

        # determine current rows (safe when word2 shorter)
        if i >= len(word):
            row_i = ()
        else:
            row_i = word[i]
        if i >= len(word2):
            row2_i = ()
        else:
            row2_i = word2[i]

        x0 = letter

        # append case (no bump in this row)
        if len(row_i) == 0 or x0 >= max(row_i):
            new_word = (*word[:i], (*row_i, x0), *word[i + 1 :])
            new_word2 = (*word2[:i], (*row2_i, letter2), *word2[i + 1 :])
            return new_word, new_word2

        # find bump
        x1 = min(a for a in row_i if a > x0)


        # normal replace + recurse into next row
        if x1 != x0 + 1 or x0 not in row_i:
            new_first_row = list(row_i)
            new_first_row[new_first_row.index(x1)] = x0
            new_word = (*word[:i], tuple(new_first_row), *word[i + 1 :])
            return NilPlactic.ed_insert_rsk(new_word, word2, x1, letter2, i=i + 1)

        # special case: continue bumping without changing current row
        return NilPlactic.ed_insert_rsk(word, word2, x1, letter2, i=i + 1)

    def __hash__(self):
        return hash(self._word)

    def __eq__(self, other):
        if isinstance(other, NilPlactic):
            return self._word == other._word
        return False




