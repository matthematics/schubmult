from ..utils._grid_print import GridPrint


class SetValuedTableau(GridPrint):
    """A semistandard set-valued tableau.

    Each box of a Young diagram holds a non-empty *set* of positive integers.
    The tableau is semistandard in the set-valued sense: reading each box as its
    set of labels, ``max`` of a box is strictly less than ``min`` of the box
    below it (columns strictly increase) and weakly less than ``min`` of the box
    to its right (rows weakly increase).

    Internally the tableau is stored as a dict ``{(row, col): tuple(labels)}``
    where each ``labels`` tuple is sorted ascending. Empty boxes are simply
    absent from the dict.
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

    def __eq__(self, other):
        if isinstance(other, SetValuedTableau):
            return self._cells == other._cells
        return NotImplemented

    def __hash__(self):
        return hash(tuple(sorted(self._cells.items())))
