from functools import cache

from .perm_lib import Permutation
from .plactic import Plactic


class NilPlactic(Plactic):
    def __init__(self, word=()):
        super().__init__(word)

    @classmethod
    @cache
    def all_skew_ed_tableaux(cls, outer_shape, bruhat_perm, inner_shape=None):
        """
        Return all skew Edelman–Greene tableaux of skew shape outer_shape/inner_shape
        (represented as row lengths) that are increasing (strictly increasing rows
        left-to-right and columns top-to-bottom) and whose row-reading word (E-G
        row word: read rows bottom-to-top, left-to-right, omitting zeros) has
        Permutation.ref_product(...) less than or equal to `bruhat_perm` in
        Bruhat order.

        Representation details:
        - outer_shape: sequence of row lengths (one int per row).
        - inner_shape: optional sequence of left offsets per row (defaults to zeros).
        - The returned tableaux are NilPlactic instances whose rows are left-aligned
          of length `outer_shape[r]`, with the first `inner_shape[r]` entries set
          to 0 to represent the inner (removed) cells of the skew shape.
        """
        # normalize shapes
        outer = tuple(int(x) for x in outer_shape)
        r = len(outer)
        if inner_shape is None:
            inner = (0,) * r
        else:
            inner = tuple(int(x) for x in inner_shape) + (0,) * (r - len(inner_shape))

        # determine a safe upper bound for entries (simple reflection indices).
        # Use permutation length as conservative bound (perm acts on 1..n).
        perm_len = len(bruhat_perm)
        max_entry = max(1, perm_len - 1)

        # build list of cell coordinates in row-major (top-to-bottom, left-to-right)
        cells = []
        for i in range(r):
            start = inner[i]
            end = outer[i]
            for j in range(start, end):
                cells.append((i, j))

        # prepare an array-of-rows initialized with zeros for inner cells
        template_rows = []
        for i in range(r):
            row_len = outer[i]
            row = [0] * row_len
            template_rows.append(row)

        results = set()

        # helper to get neighbors for increasing constraints
        def get_left_value(rows, i, j):
            if j - 1 < 0:
                return None
            return rows[i][j - 1]

        def get_top_value(rows, i, j):
            if i - 1 < 0 or j >= len(rows[i - 1]):
                return None
            return rows[i - 1][j]

        # backtracking over cells
        def backtrack(k, rows):
            if k == len(cells):
                # finished, construct NilPlactic and test Bruhat condition
                tableau = tuple(tuple([entry if entry is not None and entry != 0 else None for entry in rw]) for rw in rows)
                tpl = cls(tableau)
                if tpl.perm.bruhat_leq(bruhat_perm) and tpl.perm.inv == sum(outer_shape) - sum(inner_shape) and len(tpl.row_word) == sum(outer_shape) - sum(inner_shape):
                    results.add(tpl)
                return

            i, j = cells[k]
            # minimal value is 1, or one more than left/top neighbor if present
            left = get_left_value(rows, i, j)
            top = get_top_value(rows, i, j)
            min_val = 1
            if left is not None and left != 0:
                min_val = max(min_val, left + 1)
            if top is not None and top != 0:
                min_val = max(min_val, top + 1)

            # maximum choice is max_entry
            for val in range(min_val, max_entry + 1):
                rows[i][j] = val
                backtrack(k + 1, rows)
                # reset cell to 0 for cleanliness (not strictly necessary)
                rows[i][j] = None

        # start backtracking
        backtrack(0, template_rows)
        # cacheable return: convert to tuple for immutability
        return set(results)

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
                if up > left:
                    new_word[i][j] = up
                    i -= 1
                if up < left:
                    new_word[i][j] = left
                    j -= 1
                if up == left:
                    new_word[i][j] = up
                    new_word[i - 1][j] = new_word[i - 1][j - 1]
                    new_word[i][j - 1] = new_word[i - 1][j - 1]
                    i -= 1
                    j -= 1

        new_word[i][j] = 0

        return NilPlactic(tuple(tuple(r) for r in new_word))

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
                if down > right:
                    new_word[i][j] = right
                    j += 1
                if down == right:
                    if down == 0:
                        break
                    new_word[i][j] = down
                    new_word[i + 1][j] = new_word[i + 1][j + 1]
                    new_word[i][j + 1] = new_word[i + 1][j + 1]
                    i += 1
                    j += 1

        # Remove the last empty cell
        new_word[i].pop()
        # Remove any empty rows at the bottom
        while new_word and len(new_word[-1]) == 0:
            new_word.pop()

        return NilPlactic(tuple(tuple(r) for r in new_word))

    def bruhat_leq(self, other):
        if not isinstance(other, NilPlactic):
            return NotImplemented
        if len(self.row_word) > len(other.row_word):
            return False
        if len(self.row_word) == len(other.row_word):
            return self.row_word == other.row_word
        if not self.perm.bruhat_leq(other.perm):
            return False
        for i in range(len(other.row_word)):
            new_word = list(other.row_word[:i] + other.row_word[i + 1 :])
            new_perm = Permutation.ref_product(*new_word)
            if not self.perm.bruhat_leq(new_perm):
                continue
            return self.bruhat_leq(NilPlactic().ed_insert(*new_word))
        return True

    def hw_rc(self, length=None):
        """
        Return the highest-weight RCGraph corresponding to this NilPlactic
        tableau via Edelman–Greene insertion.
        """
        from schubmult.schub_lib.rc_graph import RCGraph

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

    def ed_insert(self, *letters):
        """Insert a letter/entry into this NilPlactic tableau and return a new Plactic."""
        new_word = self._word
        for letter in letters:
            new_word = NilPlactic._ed_insert(new_word, int(letter))
        return NilPlactic(new_word)

    @classmethod
    def ed_insert_rsk(cls, *letters):
        """Insert a letter/entry into this NilPlactic tableau and return a new Plactic."""
        word0, word2 = (), ()
        for idx, letter in enumerate(letters):
            word0, word2 = cls._ed_insert_rsk(word0, word2, int(letter), idx + 1)
        return cls(word0), Plactic(word2)

    @staticmethod
    def _ed_insert(word, letter, i=0):
        """
        Row insertion for NilPlactic. `word` is a tuple-of-rows (each a tuple).
        Returns a tuple-of-rows.
        """
        word = tuple(tuple(r) for r in word)

        # determine current rows (safe when word2 shorter)
        if i >= len(word):
            row_i = ()
        else:
            row_i = word[i]

        x0 = letter

        # append case (no bump in this row)
        if len(row_i) == 0 or x0 >= max(row_i):
            new_word = (*word[:i], (*row_i, x0), *word[i + 1 :])
            return new_word

        # find bump
        x1 = min(a for a in row_i if a > x0)

        # normal replace + recurse into next row
        if x1 != x0 + 1 or x0 not in row_i:
            new_first_row = list(row_i)
            new_first_row[new_first_row.index(x1)] = x0
            new_word = (*word[:i], tuple(new_first_row), *word[i + 1 :])
            return NilPlactic._ed_insert(new_word, x1, i=i + 1)

        # special case: continue bumping without changing current row
        return NilPlactic._ed_insert(word, x1, i=i + 1)

    @staticmethod
    def _ed_insert_rsk(word, word2, letter, letter2, i=0):
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
            return NilPlactic._ed_insert_rsk(new_word, word2, x1, letter2, i=i + 1)

        # special case: continue bumping without changing current row
        return NilPlactic._ed_insert_rsk(word, word2, x1, letter2, i=i + 1)

    # @classmethod
    # def ed_insert_rsk(cls, letters, letters2=()):
    #     """Edelman–Greene two-row insertion for NilPlactic tableaux."""
    #     word = ()
    #     word2 = ()
    #     for idx, letter in enumerate(letters):
    #         letter2 = letters2[idx] if idx < len(letters2) else None
    #         word, word2 = cls._ed_insert_rsk(word, word2, int(letter), int(letter2) if letter2 is not None else None)
    #     return cls(word), Plactic(word2)

    def __hash__(self):
        return hash(self._word)

    def __eq__(self, other):
        if isinstance(other, NilPlactic):
            return self._word == other._word
        return False
