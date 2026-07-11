
from .permutation import Permutation
from .plactic import Plactic


def _compute_inner_shape_from_grid(grid):
    """
    Compute a valid inner_shape (partition) from a grid by counting leading holes per row.
    Returns a tuple representing leading holes, ensuring it forms a valid partition
    (weakly decreasing sequence from top to bottom).
    """

    if grid.size == 0:
        return None

    # Count leading holes (None or 0) in each row
    leading_holes = []
    for row_idx in range(grid.shape[0]):
        count = 0
        for col_idx in range(grid.shape[1]):
            val = grid[row_idx, col_idx]
            if val is None or val == 0:
                count += 1
            else:
                break  # Stop at first non-hole
        leading_holes.append(count)

    # Ensure it forms a valid partition (weakly decreasing from top to bottom)
    # Process from top row (index 0) to bottom row
    # Each row must have at most as many leading holes as the row above it
    for row_idx in range(1, len(leading_holes)):
        # Current row should have <= holes than row above
        if leading_holes[row_idx] > leading_holes[row_idx - 1]:
            # If current row has more holes, reduce to match row above
            leading_holes[row_idx] = leading_holes[row_idx - 1]

    # Return None if all zeros (no inner shape)
    if all(h == 0 for h in leading_holes):
        return None

    # Strip trailing zeros
    while leading_holes and leading_holes[-1] == 0:
        leading_holes.pop()

    return tuple(leading_holes) if leading_holes else None


def _compact_grid(grid):
    """
    Compact a grid by shifting values down in each column to remove gaps.
    A gap is when there's a None/0 at position (i, j) but a value at (i-1, j).
    Returns a compacted grid with no column gaps.

    Note: Grid includes a border row/column (last row and last column), which
    should not be modified.
    """
    if grid.size == 0:
        return grid

    new_grid = grid.copy()

    # Grid has border at last row and last column - exclude them
    num_data_rows = new_grid.shape[0] - 1  # Exclude border row
    num_data_cols = new_grid.shape[1] - 1  # Exclude border column

    # For each column (excluding border), collect non-None/non-0 values from top to bottom
    for col in range(num_data_cols):
        values = []
        for row in range(num_data_rows):
            val = new_grid[row, col]
            if val is not None and val != 0:
                values.append(val)

        # Clear the column (data rows only)
        new_grid[:num_data_rows, col] = None

        # Place values from bottom up (fill from highest row index, staying below border)
        if values:
            for i, val in enumerate(reversed(values)):
                new_grid[num_data_rows - 1 - i, col] = val

    return new_grid


class IncreasingTableau(Plactic):
    def __init__(self, word=(), inner_shape=None):
        super().__init__(word, inner_shape=inner_shape)


    @property
    def perm(self):
        return Permutation.hecke_ref_product(*self.row_word)

    def hecke_insert(self, *letters):
        """Insert a letter/entry into this IncreasingTableau tableau and return a new Plactic."""
        # Convert grid to tuple of tuples format for _ed_insert (excluding border)
        word_tuples = []
        for i in range(self.rows):
            row_vals = []
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is not None and val != 0:
                    row_vals.append(val)
            if row_vals:  # Only add non-empty rows
                word_tuples.append(tuple(row_vals))
        new_word = tuple(word_tuples)

        for letter in letters:
            new_word = IncreasingTableau._hecke_insert(new_word, int(letter))
        return IncreasingTableau(new_word)

    @classmethod
    def ed_insert_rsk(cls, w1, w2):
        """Insert a letter/entry into this IncreasingTableau tableau and return a new Plactic."""
        word0, word2 = (), ()
        for a, b in zip(w1, w2):
            word0, word2 = cls._ed_insert_rsk(word0, word2, int(a), int(b))
        return cls(word0), Plactic(word2)

    @classmethod
    def ed_column_insert_rsk(cls, w1, w2):
        """Insert a letter/entry into this IncreasingTableau tableau and return a new Plactic."""
        word0, word2 = (), ()
        for a, b in zip(w1, w2):
            word0, word2 = cls._ed_insert_rsk(word0, word2, int(a), int(b))
        return cls(word0), Plactic(word2).transpose()

    @staticmethod
    def _hecke_insert(word, letter, i=0):
        """
        Row insertion for IncreasingTableau. `word` is a tuple-of-rows (each a tuple).
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
            return (*word[:i], (*row_i, x0), *word[i + 1 :])

        # find bump
        x1 = min(a for a in row_i if a > x0)

        # normal replace + recurse into next row
        if x1 != x0 + 1 or x0 not in row_i:
            new_first_row = list(row_i)
            new_first_row[new_first_row.index(x1)] = x0
            new_word = (*word[:i], tuple(new_first_row), *word[i + 1 :])
            return IncreasingTableau._hecke_insert(new_word, x1, i=i + 1)

        # special case: continue bumping without changing current row
        return IncreasingTableau._hecke_insert(word, x1, i=i + 1)

    def _to_pcells(self):
        """Return this tableau as a dict ``{(row, col): value}`` (border excluded)."""
        cells = {}
        for i in range(self.rows):
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is not None and val != 0:
                    cells[(i, j)] = int(val)
        return cells

    @classmethod
    def _from_pcells(cls, pcells):
        """Build an ``IncreasingTableau`` from a dict ``{(row, col): value}``."""
        if not pcells:
            return cls()
        max_r = max(r for (r, _) in pcells)
        rows = []
        for r in range(max_r + 1):
            cols = sorted(c for (rr, c) in pcells if rr == r)
            rows.append(tuple(pcells[(r, c)] for c in cols))
        return cls(tuple(rows))

    @staticmethod
    def _is_increasing_cells(cells):
        """Check that ``cells`` (dict ``{(row, col): value}``) is an increasing
        tableau: a valid Young diagram shape whose entries strictly increase
        along both rows and columns.
        """
        for (r, c), v in cells.items():
            # shape must remain a Young diagram
            if c > 0 and (r, c - 1) not in cells:
                return False
            if r > 0 and (r - 1, c) not in cells:
                return False
            # strictly increasing rows
            right = cells.get((r, c + 1))
            if right is not None and not v < right:
                return False
            # strictly increasing columns
            down = cells.get((r + 1, c))
            if down is not None and not v < down:
                return False
        return True

    @staticmethod
    def _hecke_column_insert_letter(pcells, x):
        """Column Hecke-insert the letter ``x`` into ``pcells``.

        ``pcells`` is a dict ``{(row, col): value}`` describing an increasing
        tableau. Returns ``(new_pcells, event)`` where ``event`` is either
        ``("grow", row, col)`` when a new box carrying the final bumped letter
        was appended at ``(row, col)``, or ``("absorb", row)`` when no box was
        added (the letter was absorbed) and the recording label belongs in the
        outer corner of ``row``.
        """
        pcells = dict(pcells)
        c = 0
        while True:
            col_rows = sorted(r for (r, cc) in pcells if cc == c)
            col_vals = [pcells[(r, c)] for r in col_rows]

            # Case 1: x is >= every entry of the current column (empty column too).
            if all(x >= y for y in col_vals):
                new_r = len(col_rows)
                cand = dict(pcells)
                cand[(new_r, c)] = x
                if IncreasingTableau._is_increasing_cells(cand):
                    # append x to the bottom of column c -> new box
                    return cand, ("grow", new_r, c)
                # cannot append: x is absorbed, tableau unchanged
                bottom_row = col_rows[-1] if col_rows else new_r
                return pcells, ("absorb", bottom_row)

            # Case 2: bump. y = min{ y' in C | y' > x }.
            yr = next(r for r in col_rows if pcells[(r, c)] > x)
            y = pcells[(yr, c)]
            cand = dict(pcells)
            cand[(yr, c)] = x
            if IncreasingTableau._is_increasing_cells(cand):
                # replace y with x only when the result stays increasing
                pcells = cand
            # proceed by inserting y into the next column
            x = y
            c += 1

    def hecke_column_insert(self, *letters):
        """Column Hecke-insert ``letters`` into a copy of this tableau.

        Returns the resulting :class:`IncreasingTableau` (the insertion tableau
        ``P``). Use :meth:`hecke_column_insert_rsk` when the set-valued
        recording tableau ``Q`` is also required.
        """
        pcells = self._to_pcells()
        for letter in letters:
            pcells, _ = IncreasingTableau._hecke_column_insert_letter(pcells, int(letter))
        return IncreasingTableau._from_pcells(pcells)

    @classmethod
    def hecke_column_insert_rsk(cls, recording, insertion):
        """Column Hecke (K-theoretic) insertion of a two-line array.

        The two-line array is given by two equal-length sequences: ``recording``
        holds the top-row labels ``k`` (weakly increasing, with the letters in
        each equal-label block strictly increasing) and ``insertion`` holds the
        bottom-row letters ``a`` that are column-inserted.

        Returns ``(P, Q)`` where ``P`` is an :class:`IncreasingTableau` (the
        insertion tableau) and ``Q`` is a
        :class:`~schubmult.combinatorics.set_valued_tableau.SetValuedTableau`,
        the semistandard *set-valued* recording tableau.
        """
        from .set_valued_tableau import SetValuedTableau

        pcells = {}
        qcells = {}
        for k, a in zip(recording, insertion):
            pcells, event = cls._hecke_column_insert_letter(pcells, int(a))
            if event[0] == "grow":
                _, r, c = event
                qcells[(r, c)] = (int(k),)
            else:
                _, r = event
                # place k in the outer (rightmost) corner of Q's row r
                corner_c = max(cc for (rr, cc) in qcells if rr == r)
                qcells[(r, corner_c)] = (*qcells[(r, corner_c)], int(k))
        return cls._from_pcells(pcells), SetValuedTableau(qcells)

    @staticmethod
    def _is_corner(cells, r, c):
        """True if ``(r, c)`` is an outer (removable) corner of the shape.

        Works for any dict keyed by ``(row, col)``; only the presence of keys is
        used, so it applies both to ``pcells`` (integer values) and to
        set-valued ``qcells`` (tuple values).
        """
        return (r, c) in cells and (r + 1, c) not in cells and (r, c + 1) not in cells

    @staticmethod
    def _reverse_hecke_column_insert(pcells, corner, alpha):
        """Reverse Hecke column insertion (Thomas–Yong / Buch–Samuel, §3.2).

        Given an increasing tableau ``pcells`` (dict ``{(row, col): value}``), a
        ``corner`` ``(r, c)`` of it, and ``alpha in {0, 1}``, produce ``(Y, x)``
        as follows. Let ``y`` be the entry in cell ``c``. If ``alpha == 1`` the
        entry ``y`` is removed (the box is deleted); if ``alpha == 0`` the box is
        kept. In either case ``y`` is reverse-inserted into the column to the
        left of ``c``.

        Whenever a value ``y`` is reverse-inserted into a column ``C``, let ``x``
        be the largest entry of ``C`` with ``x < y``. If replacing ``x`` with
        ``y`` keeps the tableau increasing this is done; in any case ``x`` is
        passed to the left. When the left-most column is reached, ``x`` becomes
        the output letter.

        Returns ``(new_pcells, x)``.
        """
        pcells = dict(pcells)
        r, c = corner
        v = pcells[(r, c)]
        if alpha == 1:
            del pcells[(r, c)]

        # reverse-insert v into the columns to the left of c, right-to-left
        cc = c - 1
        while cc >= 0:
            # x = largest entry of column cc such that x < v
            below = [(rr, pcells[(rr, cc)]) for (rr, ccc) in pcells if ccc == cc and pcells[(rr, cc)] < v]
            if not below:
                # No strictly-smaller entry in this column: the value passes
                # through unchanged (mirrors a forward step that did not bump).
                cc -= 1
                continue
            xr, x = max(below, key=lambda t: t[1])
            cand = dict(pcells)
            cand[(xr, cc)] = v
            if IncreasingTableau._is_increasing_cells(cand):
                # replace x with v only when the result stays increasing
                pcells = cand
            v = x  # pass x to the left
            cc -= 1

        return pcells, v

    def reverse_hecke_column_insert(self, corner, alpha=1):
        """Reverse Hecke-column-insert the box at ``corner`` ``(row, col)``.

        ``alpha`` is ``1`` when the corner box was created by the forward
        insertion (a ``"grow"`` step) and ``0`` when the forward step merely
        absorbed a letter at that corner. Returns ``(Y, x)`` where ``Y`` is the
        resulting :class:`IncreasingTableau` and ``x`` the reconstructed letter.
        """
        pcells = self._to_pcells()
        if not IncreasingTableau._is_corner(pcells, corner[0], corner[1]):
            raise ValueError(f"{corner} is not an outer corner of the tableau")
        new_cells, x = IncreasingTableau._reverse_hecke_column_insert(pcells, corner, int(alpha))
        return IncreasingTableau._from_pcells(new_cells), x

    @classmethod
    def hecke_column_uninsert_rsk(cls, P, Q):
        """Invert :meth:`hecke_column_insert_rsk`.

        Given an insertion tableau ``P`` (:class:`IncreasingTableau`) and a
        set-valued recording tableau ``Q`` of the same shape, reconstruct the
        two-line array as ``(recording, insertion)``.

        The recording labels are undone from largest to smallest. Within a box,
        the maximum label was the last placed there: a singleton box came from a
        ``"grow"`` step (reverse with ``alpha = 1``), while a box holding several
        labels had its largest label appended by an ``"absorb"`` step (reverse
        with ``alpha = 0``, keeping the box). Because insertion is by *columns*
        with the canonical convention that equal-label letters are inserted in
        strictly *decreasing* order, the corners sharing the maximal label are
        undone right-most column first (that being the most recently created
        box for the batch).
        """
        from .set_valued_tableau import SetValuedTableau

        pcells = P._to_pcells() if isinstance(P, IncreasingTableau) else dict(P)
        if isinstance(Q, SetValuedTableau):
            qcells = {rc: tuple(labels) for rc, labels in Q.cells.items()}
        else:
            qcells = {rc: tuple(v) for rc, v in Q.items()}

        recording = []
        insertion = []
        while qcells:
            # largest recording label still present anywhere in Q
            kmax = max(max(labels) for labels in qcells.values())
            # among removable corners containing kmax, undo the right-most
            # column first: with the canonical column-insertion convention the
            # equal-label letters were inserted in strictly decreasing order, so
            # the last box created for that batch sits in the right-most column
            corners = [rc for rc in qcells if cls._is_corner(qcells, rc[0], rc[1]) and kmax in qcells[rc]]
            r, c = max(corners, key=lambda rc: (rc[1], rc[0]))

            labels = qcells[(r, c)]
            if len(labels) == 1:
                alpha = 1
                del qcells[(r, c)]
            else:
                # kmax is the maximum, i.e. the last label placed in this box
                alpha = 0
                qcells[(r, c)] = tuple(v for v in labels if v != kmax) or labels[:-1]

            pcells, x = cls._reverse_hecke_column_insert(pcells, (r, c), alpha)
            recording.append(kmax)
            insertion.append(x)

        recording.reverse()
        insertion.reverse()
        return tuple(recording), tuple(insertion)

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
            return IncreasingTableau._ed_insert_rsk(new_word, word2, x1, letter2, i=i + 1)

        # special case: continue bumping without changing current row
        return IncreasingTableau._ed_insert_rsk(word, word2, x1, letter2, i=i + 1)

    def hw_rc(self, length):
        from schubmult.combinatorics.rc_graph import RCGraph

        seq = []
        last_spot = 0
        last_elem = -1000
        for a in self.column_word:
            if a > last_elem:
                last_elem = a
                last_spot += 1
            seq.append(last_spot)
            last_elem = a
        return RCGraph.from_reduced_compatible(self.column_word, seq).resize(length)

    @classmethod
    def from_word(cls, word):
        return cls().hecke_insert(*word)

    def right_zero_act(self, length):
        from schubmult import ASx, RCGraph, uncode
        from schubmult.utils.perm_utils import little_zero

        up_perms = ASx(self.perm, length) * ASx(uncode([0]), 1)

        word_set = set()
        rcc = self.hw_rc(length)
        # rcc = IncreasingTableau().ed_insert(*word)
        # shp = p_trans(rcc.shape)

        for perm1, _ in up_perms.keys():
            for rc in RCGraph.all_hw_rcs(perm1, length + 1, weight=(*rcc.length_vector, 0)):
                new_word = little_zero(rc.perm_word, length + 1)
                if new_word == self.column_word:
                    word_set.add(IncreasingTableau.from_word(rc.perm_word))
        return word_set

    def __mul__(self, other):
        """
        Plactic product: insert entries of `other` in row-reading order
        (top-to-bottom, left-to-right) into a copy of self.
        """
        if not isinstance(other, Plactic):
            return NotImplemented
        word = [*self.row_word, *other.row_word]
        if Permutation.ref_product(*word).inv != len(word):
            return None
        pl = IncreasingTableau()
        return pl.ed_insert(*word)
