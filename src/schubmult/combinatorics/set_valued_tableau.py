from ..utils._grid_print import GridPrint
from .crystal_graph import CrystalGraph


class SetValuedTableau(GridPrint, CrystalGraph):
    """A semistandard set-valued tableau.

    Each box of a Young diagram holds a non-empty *set* of positive integers.
    The tableau is semistandard in the set-valued sense: reading each box as its
    set of labels, ``max`` of a box is strictly less than ``min`` of the box
    below it (columns strictly increase) and weakly less than ``min`` of the box
    to its right (rows weakly increase).

    Internally the tableau is stored as a dict ``{(row, col): tuple(labels)}``
    where each ``labels`` tuple is sorted ascending. Empty boxes are simply
    absent from the dict.

    The tableau carries the type ``A`` crystal structure of Monical--Pechenik--
    Scrimshaw (*Crystal structures for symmetric Grothendieck polynomials*,
    arXiv:1807.03294): :meth:`raising_operator` is ``e_i`` and
    :meth:`lowering_operator` is ``f_i``.
    """

    _display_name = "SetValuedTableau"

    def __init__(self, cells=None):
        """Create a set-valued tableau.

        ``cells`` may be:

        - a dict ``{(row, col): iterable of labels}``, or
        - a nested sequence of rows, where each entry is an iterable of labels
          (or a single integer for a singleton box); ``None`` marks an empty
          leading (skew) box.
        """
        parsed = {}
        if cells is None:
            cells = {}
        if isinstance(cells, dict):
            for (r, c), labels in cells.items():
                parsed[(int(r), int(c))] = self._normalize_labels(labels)
        else:
            for r, row in enumerate(cells):
                for c, labels in enumerate(row):
                    if labels is None:
                        continue
                    parsed[(r, c)] = self._normalize_labels(labels)
        self._cells = parsed

    @staticmethod
    def _normalize_labels(labels):
        if isinstance(labels, int):
            labels = (labels,)
        vals = tuple(sorted(int(v) for v in labels))
        if len(vals) == 0:
            raise ValueError("set-valued boxes must be non-empty")
        return vals

    @classmethod
    def from_cells(cls, cells):
        """Build a :class:`SetValuedTableau` from a dict ``{(row, col): labels}``."""
        return cls(cells)

    @property
    def cells(self):
        """Return the underlying ``{(row, col): tuple(labels)}`` dict (a copy)."""
        return dict(self._cells)

    @property
    def rows(self):
        if not self._cells:
            return 0
        return max(r for (r, _) in self._cells) + 1

    @property
    def cols(self):
        if not self._cells:
            return 0
        return max(c for (_, c) in self._cells) + 1

    def __getitem__(self, key):
        if isinstance(key, tuple):
            return self._cells.get((int(key[0]), int(key[1])))
        raise ValueError(f"Bad indexing {key=}")

    def __iter__(self):
        for key in sorted(self._cells):
            yield self._cells[key]

    @property
    def shape(self):
        shape_list = []
        for r in range(self.rows):
            count = sum(1 for (rr, _) in self._cells if rr == r)
            shape_list.append(count)
        while shape_list and shape_list[-1] == 0:
            shape_list.pop()
        return tuple(shape_list)

    @property
    def weight(self):
        """Return the content weight: ``weight[v - 1]`` counts occurrences of ``v``."""
        counts = {}
        for labels in self._cells.values():
            for v in labels:
                counts[v] = counts.get(v, 0) + 1
        if not counts:
            return ()
        top = max(counts)
        return tuple(counts.get(v, 0) for v in range(1, top + 1))

    def add_label(self, row, col, label):
        """Return a new tableau with ``label`` added to the box at ``(row, col)``."""
        new_cells = dict(self._cells)
        existing = new_cells.get((row, col), ())
        new_cells[(row, col)] = tuple(sorted({*existing, int(label)}))
        return SetValuedTableau(new_cells)

    def is_semistandard(self):
        """Check the semistandard set-valued conditions on rows and columns."""
        for (r, c), labels in self._cells.items():
            lo, hi = labels[0], labels[-1]
            if lo <= 0:
                return False
            right = self._cells.get((r, c + 1))
            if right is not None and not hi <= right[0]:
                return False
            down = self._cells.get((r + 1, c))
            if down is not None and not hi < down[0]:
                return False
            # shape must remain a Young diagram
            if c > 0 and (r, c - 1) not in self._cells:
                return False
            if r > 0 and (r - 1, c) not in self._cells:
                return False
        return True

    # ---- crystal structure (Monical--Pechenik--Scrimshaw) --------------
    def _reading_boxes(self):
        """Row-reading order over boxes: bottom-to-top rows, left-to-right."""
        boxes = []
        for r in range(self.rows - 1, -1, -1):
            for c in range(self.cols):
                if (r, c) in self._cells:
                    boxes.append((r, c))
        return boxes

    def _bracket(self, i):
        """Signature bracketing for the pair ``(i, i + 1)``.

        Within each box the ``i + 1`` sign (an opening ``-``) is read before the
        ``i`` sign (a closing ``+``), so a box holding both is self-cancelling.
        ``+`` closes the most recent unmatched ``-``. Returns the pair
        ``(leftmost_unmatched_minus_box, rightmost_unmatched_plus_box)`` (either
        may be ``None``).
        """
        open_stack = []       # boxes with an unmatched i+1 (opening '-')
        unmatched_plus = []   # boxes with an unmatched i (closing '+')
        for rc in self._reading_boxes():
            labels = self._cells[rc]
            if (i + 1) in labels:
                open_stack.append(rc)
            if i in labels:
                if open_stack:
                    open_stack.pop()
                else:
                    unmatched_plus.append(rc)
        left_minus = open_stack[0] if open_stack else None
        right_plus = unmatched_plus[-1] if unmatched_plus else None
        return left_minus, right_plus

    def lowering_operator(self, i):
        """Crystal lowering operator ``f_i`` (Definition 3.1 of MPS).

        Acts on the box ``b`` of the rightmost uncancelled ``+`` (an ``i``). If
        the box ``b_right`` immediately to the right of ``b`` contains ``i``,
        remove ``i`` from ``b_right`` and add ``i + 1`` to ``b``; otherwise
        replace ``i`` with ``i + 1`` in ``b``. Returns ``None`` if undefined.
        """
        if i < 1:
            return None
        _, b = self._bracket(i)
        if b is None:
            return None
        cells = {k: set(v) for k, v in self._cells.items()}
        r, c = b
        b_right = (r, c + 1)
        if b_right in cells and i in cells[b_right]:
            cells[b_right].discard(i)
            cells[b].add(i + 1)
        else:
            cells[b].discard(i)
            cells[b].add(i + 1)
        return SetValuedTableau({k: tuple(sorted(v)) for k, v in cells.items() if v})

    def raising_operator(self, i):
        """Crystal raising operator ``e_i`` (Definition 3.1 of MPS).

        Acts on the box ``b`` of the leftmost uncancelled ``-`` (an ``i + 1``).
        If the box ``b_left`` immediately to the left of ``b`` contains
        ``i + 1``, remove ``i + 1`` from ``b_left`` and add ``i`` to ``b``;
        otherwise replace ``i + 1`` with ``i`` in ``b``. Returns ``None`` if
        undefined.
        """
        if i < 1:
            return None
        b, _ = self._bracket(i)
        if b is None:
            return None
        cells = {k: set(v) for k, v in self._cells.items()}
        r, c = b
        b_left = (r, c - 1)
        if b_left in cells and (i + 1) in cells[b_left]:
            cells[b_left].discard(i + 1)
            cells[b].add(i)
        else:
            cells[b].discard(i + 1)
            cells[b].add(i)
        return SetValuedTableau({k: tuple(sorted(v)) for k, v in cells.items() if v})

    @property
    def crystal_weight(self):
        """Content weight ``(#1, #2, ...)`` (alias of :attr:`weight`)."""
        return self.weight

    def crystal_length(self):
        """Upper bound on crystal operator indices (matches ``Plactic``)."""
        return 100

    def __eq__(self, other):
        if isinstance(other, SetValuedTableau):
            return self._cells == other._cells
        return NotImplemented

    def __hash__(self):
        return hash(tuple(sorted(self._cells.items())))
