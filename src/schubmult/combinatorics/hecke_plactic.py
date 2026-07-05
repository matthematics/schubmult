from .permutation import Permutation
from .plactic import Plactic


class HeckePlactic(Plactic):
    r"""Insertion tableau for *size-preserving Hecke column insertion*.

    This is a "made up" variant of Hecke (K-theoretic) column insertion that
    behaves like Edelman--Greene insertion in that every inserted letter adds
    exactly **one** box (the shape/size is preserved: ``#boxes == word length``),
    while still preserving the Hecke product / K-Knuth equivalence class of the
    reading word.

    Ordinary Hecke insertion is *not* size preserving: when a letter cannot
    extend the shape it is *absorbed* (no box is added). Here we never absorb --
    instead the letter is carried forward to the next column until it can be
    placed as a new box. The price is that the resulting tableau ``P`` is only
    **column-strict semistandard** (entries strictly increase down columns and
    weakly increase along rows) rather than a strict increasing tableau: a value
    may repeat within a row.

    The recording tableau ``Q`` is an ordinary (single-valued) semistandard
    :class:`~schubmult.combinatorics.plactic.Plactic` tableau.
    """

    _display_name = "HeckePlactic"

    def __init__(self, word=(), inner_shape=None):
        super().__init__(word, inner_shape=inner_shape)

    @property
    def perm(self):
        """Hecke product of the row-reading word (the preserved K-Knuth class)."""
        return Permutation.hecke_ref_product(*self.row_word)

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
        """Build a :class:`HeckePlactic` from a dict ``{(row, col): value}``."""
        if not pcells:
            return cls()
        max_r = max(r for (r, _) in pcells)
        rows = []
        for r in range(max_r + 1):
            cols = sorted(c for (rr, c) in pcells if rr == r)
            rows.append(tuple(pcells[(r, c)] for c in cols))
        return cls(tuple(rows))

    @staticmethod
    def _is_semistandard_cells(cells):
        """Check that ``cells`` (dict ``{(row, col): value}``) is a column-strict
        semistandard tableau: a valid Young diagram shape whose entries strictly
        increase down columns and weakly increase along rows.
        """
        for (r, c), v in cells.items():
            # shape must remain a Young diagram
            if c > 0 and (r, c - 1) not in cells:
                return False
            if r > 0 and (r - 1, c) not in cells:
                return False
            # rows weakly increasing
            right = cells.get((r, c + 1))
            if right is not None and not v <= right:
                return False
            # columns strictly increasing
            down = cells.get((r + 1, c))
            if down is not None and not v < down:
                return False
        return True

    @staticmethod
    def _column_insert_letter(pcells, x):
        """Size-preserving Hecke column-insert the letter ``x`` into ``pcells``.

        ``pcells`` is a dict ``{(row, col): value}`` describing a column-strict
        semistandard tableau. Exactly one new box is always added. Returns
        ``(new_pcells, (row, col))`` giving the position of that new box.
        """
        pcells = dict(pcells)
        c = 0
        while True:
            col_rows = sorted(r for (r, cc) in pcells if cc == c)
            col_vals = [pcells[(r, c)] for r in col_rows]

            # Case 1: x >= every entry of the current column (empty column too).
            if all(x >= y for y in col_vals):
                # Append only when strictly greater than the bottom entry, so the
                # column stays strictly increasing.
                if not col_vals or x > col_vals[-1]:
                    new_r = len(col_rows)
                    pcells[(new_r, c)] = x
                    return pcells, (new_r, c)
                # x equals the bottom entry: do NOT absorb -- carry x forward to
                # the next column so a box is still created (value may repeat).
                c += 1
                continue

            # Case 2: bump. y = min{ y' in C | y' > x }.
            yr = next(r for r in col_rows if pcells[(r, c)] > x)
            y = pcells[(yr, c)]
            cand = dict(pcells)
            cand[(yr, c)] = x
            if HeckePlactic._is_semistandard_cells(cand):
                # replace y with x only when the result stays semistandard
                pcells = cand
            # proceed by inserting y into the next column
            x = y
            c += 1

    def hecke_insert(self, *letters):
        """Size-preserving Hecke column-insert ``letters`` into a copy of self.

        Returns the resulting :class:`HeckePlactic` (the insertion tableau ``P``).
        Use :meth:`hecke_insert_rsk` when the recording tableau ``Q`` is required.
        """
        pcells = self._to_pcells()
        for letter in letters:
            pcells, _ = HeckePlactic._column_insert_letter(pcells, int(letter))
        return HeckePlactic._from_pcells(pcells)

    @classmethod
    def from_word(cls, word):
        return cls().hecke_insert(*word)

    @classmethod
    def hecke_insert_rsk(cls, recording, insertion):
        """Size-preserving Hecke column insertion of a two-line array.

        ``recording`` holds the top-row labels ``k`` (weakly increasing) and
        ``insertion`` holds the bottom-row letters ``a`` that are column-inserted.

        Returns ``(P, Q)`` where ``P`` is a :class:`HeckePlactic` insertion
        tableau and ``Q`` is a single-valued semistandard
        :class:`~schubmult.combinatorics.plactic.Plactic` recording tableau.
        Because insertion is size preserving, ``P`` and ``Q`` have the same shape
        and ``Q`` records the label ``k`` in each newly created box.
        """
        pcells = {}
        qcells = {}
        for k, a in zip(recording, insertion):
            pcells, (r, c) = cls._column_insert_letter(pcells, int(a))
            qcells[(r, c)] = int(k)
        return cls._from_pcells(pcells), cls._plactic_from_pcells(qcells)

    @staticmethod
    def _plactic_from_pcells(qcells):
        """Build a :class:`Plactic` from a dict ``{(row, col): value}``."""
        if not qcells:
            return Plactic()
        max_r = max(r for (r, _) in qcells)
        rows = []
        for r in range(max_r + 1):
            cols = sorted(c for (rr, c) in qcells if rr == r)
            rows.append(tuple(qcells[(r, c)] for c in cols))
        return Plactic(tuple(rows))

    @staticmethod
    def _is_corner(pcells, r, c):
        """True if ``(r, c)`` is an outer (removable) corner of the shape."""
        return (r, c) in pcells and (r + 1, c) not in pcells and (r, c + 1) not in pcells

    @staticmethod
    def _reverse_column_insert(pcells, corner):
        """Reverse one size-preserving Hecke column insertion.

        Given ``pcells`` (dict ``{(row, col): value}``) and an outer ``corner``
        ``(r, c)`` that was the box created by a forward insertion, remove that
        box and reverse-insert its value leftward. Returns ``(new_pcells, x)``
        where ``x`` is the reconstructed inserted letter.

        Because size-preserving insertion always creates a box, the corner entry
        is always removed (``alpha = 1`` in the reverse-Hecke terminology).
        """
        pcells = dict(pcells)
        r, c = corner
        v = pcells.pop((r, c))  # remove the corner box, carry its value left

        cc = c - 1
        while cc >= 0:
            # x = largest entry of column cc that is strictly less than v
            below = [(rr, pcells[(rr, cc)]) for (rr, ccc) in pcells if ccc == cc and pcells[(rr, cc)] < v]
            if not below:
                # No strictly-smaller entry: this column was a carry-through in the
                # forward pass (the value passed by unchanged). Pass v left as-is.
                cc -= 1
                continue
            xr, x = max(below, key=lambda t: t[1])
            cand = dict(pcells)
            cand[(xr, cc)] = v
            if HeckePlactic._is_semistandard_cells(cand):
                # replace x with v only when the result stays semistandard
                pcells = cand
            v = x  # pass x to the left
            cc -= 1

        return pcells, v

    def reverse_hecke_insert(self, corner):
        """Reverse-insert the box at ``corner`` ``(row, col)``.

        Returns ``(Y, x)`` where ``Y`` is the resulting :class:`HeckePlactic`
        and ``x`` is the reconstructed letter.
        """
        pcells = self._to_pcells()
        if not HeckePlactic._is_corner(pcells, corner[0], corner[1]):
            raise ValueError(f"{corner} is not an outer corner of the tableau")
        new_cells, x = HeckePlactic._reverse_column_insert(pcells, corner)
        return HeckePlactic._from_pcells(new_cells), x

    @classmethod
    def hecke_uninsert_rsk(cls, P, Q):
        """Invert :meth:`hecke_insert_rsk`.

        Given an insertion tableau ``P`` (:class:`HeckePlactic`) and a recording
        tableau ``Q`` (:class:`~schubmult.combinatorics.plactic.Plactic`) of the
        same shape, reconstruct the two-line array as ``(recording, insertion)``.
        """
        pcells = P._to_pcells() if isinstance(P, HeckePlactic) else dict(P)
        qcells = {}
        for i in range(Q.rows):
            for j in range(Q.cols):
                val = Q[i, j]
                if val is not None and val != 0:
                    qcells[(i, j)] = int(val)

        recording = []
        insertion = []
        while qcells:
            # Remove boxes in reverse insertion order. Forward insertion is
            # *column* insertion, so cells of equal recording label form a
            # horizontal strip and were created left-to-right; undo them
            # right-to-left. Among corners choose the maximal recording label,
            # then the right-most column.
            corners = [(r, c) for (r, c) in qcells if cls._is_corner(qcells, r, c)]
            kmax = max(qcells[rc] for rc in corners)
            r, c = max((rc for rc in corners if qcells[rc] == kmax), key=lambda rc: rc[1])

            pcells, x = cls._reverse_column_insert(pcells, (r, c))
            recording.append(qcells.pop((r, c)))
            insertion.append(x)

        recording.reverse()
        insertion.reverse()
        return tuple(recording), tuple(insertion)

