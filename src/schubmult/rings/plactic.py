# filepath: /home/matthematics/schubmult/src/schubmult/rings/plactic.py
from schubmult.perm_lib import Permutation

from ._grid_print import GridPrint


class Plactic(GridPrint):
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
        return tuple(a for row in reversed(self._word) for a in row)

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
        for letter in word:
            pl = pl.rs_insert(letter)
        return pl

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
            new_row = tuple([*row_i, x0])
            return tuple([*word[:i], new_row, *word[i + 1 :]])

        # otherwise bump the smallest entry > x0
        x1 = min(a for a in row_i if a > x0)
        new_first_row = list(row_i)
        idx = new_first_row.index(x1)
        new_first_row[idx] = x0
        new_word = tuple([*word[:i], tuple(new_first_row), *word[i + 1 :]])
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


class NilPlactic(Plactic):

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
            return (tuple([*row, letter]),) + tuple(word[1:])

        # equal-case handling
        if row[index] == letter:
            if index < len(row) - 1 and row[index + 1] == letter + 1:
                # bump into next row
                return (row,) + NilPlactic._ed_insert(tuple(word[1:]), letter + 1)
            # skip equal-run
            while index < len(row) and row[index] == letter:
                index += 1
            if index == len(row):
                return (tuple([*row, letter]),) + tuple(word[1:])

        # bump case
        new_row = list(row)
        bump = new_row[index]
        new_row[index] = letter
        return (tuple(new_row),) + NilPlactic._ed_insert(tuple(word[1:]), bump)

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
            new_word = tuple([*word[:i], tuple([*row_i, x0]), *word[i + 1 :]])
            new_word2 = tuple([*word2[:i], tuple([*row2_i, letter2]), *word2[i + 1 :]])
            return new_word, new_word2

        # find bump
        x1 = min(a for a in row_i if a > x0)

        # normal replace + recurse into next row
        if x1 != x0 + 1 or x0 not in row_i:
            new_first_row = list(row_i)
            new_first_row[new_first_row.index(x1)] = x0
            new_word = tuple([*word[:i], tuple(new_first_row), *word[i + 1 :]])
            return NilPlactic.ed_insert_rsk(new_word, word2, x1, letter2, i=i + 1)

        # special case: continue bumping without changing current row
        return NilPlactic.ed_insert_rsk(word, word2, x1, letter2, i=i + 1)

    def __eq__(self, other):
        if isinstance(other, NilPlactic):
            return self._word == other._word
        return False




