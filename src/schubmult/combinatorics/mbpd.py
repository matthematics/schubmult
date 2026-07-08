r"""Implementation of the MBPD <-> RCP <-> WCGraph bijection from
``writing/mbpd.solve.tex`` (Theorem "T: main").

Geometry (paper convention):
  * n x n grid, rows top->bottom (1..n), columns left->right (1..n).
  * Pipes ENTER from the RIGHT border and EXIT to the SOUTH border.
  * No pipe enters from the top, no pipe leaves to the left.

Seven tiles, described by the set of directions in which the pipe connects
(N=up, E=right, S=down, W=left) together with a "marked" bit:

    B  blank        {}            heavy
    H  horizontal   {E, W}
    V  vertical     {N, S}
    P  plus/cross   {N, E, S, W}
    R  R-elbow      {E, S}
    J  J-elbow      {N, W}
    M  marked J     {N, W}        heavy

A tile is *heavy* iff it is B or M.

The bijection ``Phi = MBPD.phi`` and its inverse ``Psi = RCP.psi`` preserve the
associated permutation and weight.  ``WCGraph.to_mbpd`` / ``WCGraph.from_mbpd``
expose the composition with the trivial ``RCP <-> WCGraph`` repackaging.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

# Directions as single-character strings.
N, E, S, W = "N", "E", "S", "W"
_OPP = {N: S, S: N, E: W, W: E}

# Canonical tile name  ->  frozenset of connection directions.
_CONN = {
    "B": frozenset(),
    "H": frozenset({E, W}),
    "V": frozenset({N, S}),
    "P": frozenset({N, E, S, W}),
    "R": frozenset({E, S}),
    "J": frozenset({N, W}),
    "M": frozenset({N, W}),  # marked J
}
# frozenset -> unmarked tile name (M shares connections with J).
_NAME = {v: k for k, v in _CONN.items() if k != "M"}

_HEAVY = {"B", "M"}


def tile_name(conn: frozenset, marked: bool) -> str:
    """Return the canonical tile name for a connection set / mark."""
    base = _NAME[conn]
    if marked:
        if base != "J":
            raise ValueError(f"only J tiles may be marked, got {base}")
        return "M"
    return base


def is_heavy(conn: frozenset, marked: bool) -> bool:
    return marked or len(conn) == 0


# ASCII glyphs for pretty printing (approximate the pipe drawing).
_GLYPH = {
    "B": "·",
    "H": "─",
    "V": "│",
    "P": "┼",
    "R": "┌",
    "J": "┘",
    "M": "◹",  # marked J
}


class MBPD:
    """A marked bumpless pipedream on an ``n x n`` grid.

    Internally we store, for every cell ``(i, j)`` (1-indexed), the connection
    set ``conn[i][j]`` (a frozenset of ``N/E/S/W``) and a boolean ``marked``.
    """

    __slots__ = ("_conn", "_marked", "n")

    def __init__(self, n: int, conn, marked):
        self.n = n
        # store as tuples-of-tuples for hashability; index [i-1][j-1]
        self._conn = tuple(tuple(row) for row in conn)
        self._marked = tuple(tuple(bool(x) for x in row) for row in marked)

    # ---- basic accessors -------------------------------------------------
    def conn(self, i: int, j: int) -> frozenset:
        return self._conn[i - 1][j - 1]

    def marked(self, i: int, j: int) -> bool:
        return self._marked[i - 1][j - 1]

    def tile(self, i: int, j: int) -> str:
        return tile_name(self.conn(i, j), self.marked(i, j))

    def heavy(self, i: int, j: int) -> bool:
        return is_heavy(self.conn(i, j), self.marked(i, j))

    def connects(self, i: int, j: int, d: str) -> bool:
        return d in self._conn[i - 1][j - 1]

    # ---- hashing / equality ---------------------------------------------
    def _key(self):
        return (self.n, self._conn, self._marked)

    def __eq__(self, other):
        return isinstance(other, MBPD) and self._key() == other._key()

    def __hash__(self):
        return hash(self._key())

    # ---- construction ----------------------------------------------------
    @classmethod
    def from_tiles(cls, grid) -> MBPD:
        """Build from an ``n x n`` grid of tile-name strings."""
        n = len(grid)
        conn = [[_CONN[grid[i][j]] for j in range(n)] for i in range(n)]
        marked = [[grid[i][j] == "M" for j in range(n)] for i in range(n)]
        return cls(n, conn, marked)

    def to_tiles(self):
        return [[self.tile(i, j) for j in range(1, self.n + 1)] for i in range(1, self.n + 1)]

    def with_tile(self, i: int, j: int, name: str) -> MBPD:
        """Return a copy with cell ``(i, j)`` set to tile ``name``."""
        conn = [list(row) for row in self._conn]
        marked = [list(row) for row in self._marked]
        conn[i - 1][j - 1] = _CONN[name]
        marked[i - 1][j - 1] = name == "M"
        return MBPD(self.n, conn, marked)

    # ---- Rothe diagram ---------------------------------------------------
    @classmethod
    def rothe(cls, perm, n=None) -> MBPD:
        r"""The Rothe MBPD ``D_w``: pipe ``i`` turns only at ``(i, w(i))``.

        Connections (derived and verified against the paper's conventions):
          E : j >= w(i)                 (horizontal present / border on right)
          W : j > w(i)
          S : w^{-1}(j) <= i            (vertical present going down)
          N : w^{-1}(j) < i
        """
        w = list(perm)
        n = max(len(w), n or 0)
        w = w + list(range(len(w) + 1, n + 1))  # pad identity to size n
        winv = [0] * (n + 1)
        for i in range(1, n + 1):
            winv[w[i - 1]] = i
        conn = [[None] * n for _ in range(n)]
        marked = [[False] * n for _ in range(n)]
        for i in range(1, n + 1):
            wi = w[i - 1]
            for j in range(1, n + 1):
                dirs = set()
                if j >= wi:
                    dirs.add(E)
                if j > wi:
                    dirs.add(W)
                if winv[j] <= i:
                    dirs.add(S)
                if winv[j] < i:
                    dirs.add(N)
                conn[i - 1][j - 1] = frozenset(dirs)
        return cls(n, conn, marked)

    # ---- validity --------------------------------------------------------
    def validity_errors(self):
        """Return a list of human-readable validity problems (empty if valid)."""
        errs = []
        n = self.n
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                c = self._conn[i - 1][j - 1]
                # tile must be one of the 7 recognized connection sets
                if c not in _NAME:
                    errs.append(f"({i},{j}): unknown connection set {set(c)}")
                    continue
                if self._marked[i - 1][j - 1] and c != _CONN["J"]:
                    errs.append(f"({i},{j}): marked non-J tile")
                # internal adjacency: E of (i,j) matches W of (i,j+1)
                if j < n:
                    right = self._conn[i - 1][j]
                    if (E in c) != (W in right):
                        errs.append(f"({i},{j})-({i},{j + 1}): E/W mismatch")
                # S of (i,j) matches N of (i+1,j)
                if i < n:
                    below = self._conn[i][j - 1]
                    if (S in c) != (N in below):
                        errs.append(f"({i},{j})-({i + 1},{j}): S/N mismatch")
        # border conditions
        for i in range(1, n + 1):
            if E not in self._conn[i - 1][n - 1]:
                errs.append(f"row {i}: no pipe entering from the right")
            if W in self._conn[i - 1][0]:
                errs.append(f"row {i}: pipe leaves to the left")
        for j in range(1, n + 1):
            if S not in self._conn[n - 1][j - 1]:
                errs.append(f"col {j}: no pipe leaving to the south")
            if N in self._conn[0][j - 1]:
                errs.append(f"col {j}: pipe enters from the top")
        return errs

    def is_valid(self) -> bool:
        return not self.validity_errors()

    # ---- associated permutation -----------------------------------------
    def perm(self):
        r"""Associated permutation ``w`` where the pipe entering the right of
        row ``i`` leaves the bottom of column ``w(i)``.

        Double crossings between an already-crossed pair are ignored (the
        Demazure convention): the *first* time a pair of pipes meets at a ``P``
        tile they cross, every later meeting is a *bump* (the ``P`` acts as
        superimposed elbows ``R``+``J``, i.e. ``E<->S`` and ``N<->W``).

        We compute the routing by a fixed point.  Start assuming every ``P`` is
        a genuine crossing (straight through) and trace all pipes.  Whenever a
        pair of pipes meets at more than one ``P`` tile we mark all but their
        first meeting as bumps, then re-trace.  This is repeated until the set
        of bump tiles stabilises; the exit columns then give ``w``.
        """
        from schubmult.combinatorics.permutation import Permutation

        n = self.n
        bump = set()  # P-cells (i, j) that behave as bumps instead of crosses
        routes = {}
        for _ in range(n * n + 5):
            routes = {p: self._trace_route(p, bump) for p in range(1, n + 1)}
            # For every P cell record the (pipe, orientation, arc-length) of the
            # (at most two) strands that pass through it.
            cellpass = {}
            for p, route in routes.items():
                for arc, (i, j, orient) in enumerate(route["cells"]):
                    cellpass.setdefault((i, j), []).append((p, orient, arc))
            # A genuine meeting: two DISTINCT pipes cross (one H, one V) at a P.
            pair_meets = {}  # frozenset({p, q}) -> list of (cell, arc_of_min_pipe)
            for cell, lst in cellpass.items():
                if self.tile(cell[0], cell[1]) != "P" or len(lst) != 2:
                    continue
                (p1, o1, a1), (p2, o2, a2) = lst
                if p1 == p2 or o1 == o2:
                    continue
                pr = frozenset({p1, p2})
                arc = a1 if p1 < p2 else a2
                pair_meets.setdefault(pr, []).append((cell, arc))
            new_bump = set()
            for pr, meets in pair_meets.items():
                if len(meets) > 1:
                    meets.sort(key=lambda x: x[1])
                    for cell, _arc in meets[1:]:
                        new_bump.add(cell)
            if new_bump == bump:
                break
            bump = new_bump
        w = [0] * n
        for p in range(1, n + 1):
            ei, ej = routes[p]["exit"]
            w[p - 1] = ej if ei == n + 1 else (ej or 0)
        return Permutation(w)

    def _trace_route(self, p, bump):
        """Trace pipe ``p`` (entering the right border of row ``p``).

        ``bump`` is a set of ``P`` cells that act as bumps (``E<->S``, ``N<->W``)
        rather than crossings (straight through).  Returns a dict with:
          ``cells``: list of ``(i, j, orient)`` visited, orient in {"H", "V"};
          ``exit``:  ``(i, j)`` of the border cell it leaves through.
        """
        n = self.n
        i, j, came = p, n, E  # entered cell (p, n) from the East border
        cells = []
        for _ in range(4 * n * n + 10):
            orient = "H" if came in (E, W) else "V"
            cells.append((i, j, orient))
            sides = set(self._conn[i - 1][j - 1])
            if came not in sides:
                sides = sides | {came}
            others = sides - {came}
            if len(others) == 1:
                leave = next(iter(others))
            elif (i, j) in bump:
                # bump geometry: E<->S and N<->W
                leave = {E: S, S: E, N: W, W: N}[came]
            else:
                # straight-through crossing
                leave = _OPP[came]
            if leave == E:
                next_i, next_j, next_came = i, j + 1, W
            elif leave == W:
                next_i, next_j, next_came = i, j - 1, E
            elif leave == S:
                next_i, next_j, next_came = i + 1, j, N
            else:
                next_i, next_j, next_came = i - 1, j, S
            i, j, came = next_i, next_j, next_came
            if i < 1 or i > n or j < 1 or j > n:
                return {"cells": cells, "exit": (i, j)}
        return {"cells": cells, "exit": (i, j)}

    # ---- weight ----------------------------------------------------------
    def weight(self):
        """``wt(D) = (m_1, ..., m_n)`` with ``m_i`` = #heavy tiles in row ``i``."""
        return tuple(
            sum(1 for j in range(1, self.n + 1) if self.heavy(i, j))
            for i in range(1, self.n + 1)
        )

    def num_heavy(self) -> int:
        return sum(self.weight())

    def heavy_cells(self):
        return [
            (i, j)
            for i in range(1, self.n + 1)
            for j in range(1, self.n + 1)
            if self.heavy(i, j)
        ]

    # ---- row / segment predicates ---------------------------------------
    def is_pipe_segment(self, r: int, b: int, c: int) -> bool:
        r"""``D_{r,[b,c]}`` is a pipe segment: a connected horizontal run.

        For ``b == c`` a single non-blank tile qualifies.  For ``b < c`` the
        interior ``D_{r,[b+1,c-1]}`` must be all ``H``/``P`` tiles and the run
        must actually connect: ``(r,b)`` connects east and ``(r,c)`` connects
        west (so the paper's inference that the endpoints connect toward the
        interior holds).
        """
        if b > c:
            return False
        if b == c:
            return self.tile(r, b) != "B"
        if not self.connects(r, b, E) or not self.connects(r, c, W):
            return False
        for j in range(b + 1, c):
            if self.tile(r, j) not in ("H", "P"):
                return False
        return True

    def rj_subsequence(self, r: int, lo: int, hi: int):
        """The RJ subsequence (list of 'R'/'J') of light tiles in ``[lo,hi]``.

        Returns ``None`` if any tile in the range is heavy (not a light seq).
        """
        if lo > hi:
            return []
        out = []
        for j in range(lo, hi + 1):
            t = self.tile(r, j)
            if t in ("B", "M"):
                return None  # not a light sequence
            if t in ("R", "J"):
                out.append(t)
        return out

    def light_seq_type(self, r: int, lo: int, hi: int):
        r"""Classify the light sequence ``D_{r,[lo,hi]}`` as one of
        ``'paired'``, ``'J'``, ``'R'``, ``'JR'`` (or ``None`` if it is not a
        light sequence).  Empty range counts as ``'paired'``."""
        seq = self.rj_subsequence(r, lo, hi)
        if seq is None:
            return None
        lead = 1 if seq and seq[0] == "J" else 0
        tail = 1 if seq and seq[-1] == "R" else 0
        core = seq[lead: len(seq) - tail] if tail else seq[lead:]
        # core must be (R J)*
        if len(core) % 2 != 0:
            return None
        for k in range(0, len(core), 2):
            if core[k] != "R" or core[k + 1] != "J":
                return None
        return {(0, 0): "paired", (1, 0): "J", (0, 1): "R", (1, 1): "JR"}[(lead, tail)]

    def is_doublecross(self, r: int, b: int, d: int) -> bool:
        r"""``D_{[r,r+1],[b,d]}`` is a doublecross: both rows are pipe segments,
        ``D_{r,b}=R`` and ``D_{r+1,d}=J``."""
        if not (1 <= r < self.n and 1 <= b < d <= self.n):
            return False
        return (
            self.is_pipe_segment(r, b, d)
            and self.is_pipe_segment(r + 1, b, d)
            and self.tile(r, b) == "R"
            and self.tile(r + 1, d) == "J"
        )

    # ---- droop / undroop -------------------------------------------------
    def admits_droop(self, r: int, b: int, d: int) -> bool:
        if not (1 <= r < self.n and 1 <= b < d <= self.n):
            return False
        for i in (r, r + 1):
            for j in range(b, d + 1):
                if self.heavy(i, j) and not (i == r + 1 and j == d and self.tile(i, j) == "B"):
                    return False
        if not self.is_pipe_segment(r, b, d):
            return False
        if self.light_seq_type(r + 1, b + 1, d - 1) != "paired":
            return False
        # The four corner shapes forced by an admissible droop (paper's Remark):
        #   (r,b) connects E and S        [R or P]
        #   (r,d) connects W, not S       [J or H]
        #   (r+1,b) connects N, not E     [J or V]
        #   (r+1,d) not N, not W          [B or R]
        if not (self.connects(r, b, E) and self.connects(r, b, S)):
            return False
        if not (self.connects(r, d, W) and not self.connects(r, d, S)):
            return False
        if not (self.connects(r + 1, b, N) and not self.connects(r + 1, b, E)):
            return False
        if self.connects(r + 1, d, N) or self.connects(r + 1, d, W):
            return False
        return True

    def admits_undroop(self, r: int, b: int, d: int) -> bool:
        if not (1 <= r < self.n and 1 <= b < d <= self.n):
            return False
        for i in (r, r + 1):
            for j in range(b, d + 1):
                if self.heavy(i, j) and not (i == r and j == b and self.tile(i, j) == "B"):
                    return False
        if not self.is_pipe_segment(r + 1, b, d):
            return False
        if self.light_seq_type(r, b + 1, d - 1) != "paired":
            return False
        # Corner shapes for an admissible undroop (inverse of the droop Remark):
        #   (r,b) not E, not S            [B or J]
        #   (r,d) connects S, not W       [V or R]
        #   (r+1,b) connects E, not N     [H or R]
        #   (r+1,d) connects N and W      [J or P]
        if self.connects(r, b, E) or self.connects(r, b, S):
            return False
        if not (self.connects(r, d, S) and not self.connects(r, d, W)):
            return False
        if not (self.connects(r + 1, b, E) and not self.connects(r + 1, b, N)):
            return False
        if not (self.connects(r + 1, d, N) and self.connects(r + 1, d, W)):
            return False
        return True

    def _swap_rows_core(self, r: int, b: int, d: int) -> MBPD:
        r"""The common involution underlying both droop and undroop on the
        window ``[r,r+1] x [b,d]``.

        Swap the *interior* horizontal edges between the two rows, keep all
        external edges fixed, and re-solve the vertical mid-edges so that every
        tile in the window is valid.  This is an involution, so
        ``undroop = droop^{-1}`` is computed by the very same routine.
        """
        cols = list(range(b, d + 1))
        VT = {j: self.connects(r, j, N) for j in cols}       # external, above row r
        VB = {j: self.connects(r + 1, j, S) for j in cols}   # external, below row r+1
        # external horizontals on the sides of the window (preserved by droop)
        W_r_b = self.connects(r, b, W)
        E_r_d = self.connects(r, d, E)
        W_r1_b = self.connects(r + 1, b, W)
        E_r1_d = self.connects(r + 1, d, E)
        # interior horizontal edges (edge between column j and j+1), j in [b, d-1]
        HR = {j: self.connects(r, j, E) for j in range(b, d)}
        HL = {j: self.connects(r + 1, j, E) for j in range(b, d)}
        newHR = dict(HL)   # row r interior horizontals <- old row r+1
        newHL = dict(HR)   # row r+1 interior horizontals <- old row r

        conn = [list(row) for row in self._conn]
        marked = [list(row) for row in self._marked]
        for j in cols:
            Wr = W_r_b if j == b else newHR[j - 1]
            Er = E_r_d if j == d else newHR[j]
            Wr1 = W_r1_b if j == b else newHL[j - 1]
            Er1 = E_r1_d if j == d else newHL[j]
            chosen = None
            for x in (0, 1):  # x = vertical mid-edge between (r,j) and (r+1,j)
                top = set()
                if VT[j]:
                    top.add(N)
                if x:
                    top.add(S)
                if Er:
                    top.add(E)
                if Wr:
                    top.add(W)
                bot = set()
                if x:
                    bot.add(N)
                if VB[j]:
                    bot.add(S)
                if Er1:
                    bot.add(E)
                if Wr1:
                    bot.add(W)
                if frozenset(top) in _NAME and frozenset(bot) in _NAME:
                    chosen = (frozenset(top), frozenset(bot))
                    break
            if chosen is None:
                raise ValueError(f"_swap_rows_core: no valid vertical at column {j} of window ({r},[{b},{d}])")
            conn[r - 1][j - 1] = chosen[0]
            conn[r][j - 1] = chosen[1]
            marked[r - 1][j - 1] = False
            marked[r][j - 1] = False
        return MBPD(self.n, conn, marked)

    def droop(self, r: int, b: int, d: int) -> MBPD:
        if not self.admits_droop(r, b, d):
            raise ValueError(f"({r},[{b},{d}])-droop not admitted")
        return self._swap_rows_core(r, b, d)

    def undroop(self, r: int, b: int, d: int) -> MBPD:
        if not self.admits_undroop(r, b, d):
            raise ValueError(f"({r},[{b},{d}])-undroop not admitted")
        return self._swap_rows_core(r, b, d)

    # ---- F-targets and F-moves ------------------------------------------
    def _f1(self, r: int, c: int) -> bool:
        """(f1): ``(r,c)`` is the rightmost heavy tile in row ``r``."""
        if not self.heavy(r, c):
            return False
        for j in range(c + 1, self.n + 1):
            if self.heavy(r, j):
                return False
        return True

    def is_f_target(self, r: int, c: int) -> bool:
        if not (1 <= r < self.n):
            return False
        if not self._f1(r, c):
            return False
        cp = None  # (f2) min c' > c with D_{r+1,c'} = J
        for j in range(c + 1, self.n + 1):
            if self.tile(r + 1, j) == "J":
                cp = j
                break
        if cp is None:
            return False
        for j in range(1, cp):  # (f3)
            if self.heavy(r + 1, j):
                return False
        return True

    def is_fstar_target(self, r: int, c: int) -> bool:
        if not (1 <= r < self.n):
            return False
        if not self._f1(r, c):
            return False
        # (f*2): no J or M in row r+1 right of c
        for j in range(c + 1, self.n + 1):
            if self.tile(r + 1, j) in ("J", "M"):
                return False
        cp = None  # c' = max col with D_{r,c'} = R
        for j in range(self.n, 0, -1):
            if self.tile(r, j) == "R":
                cp = j
                break
        if cp is None or cp <= c:
            return False
        for j in range(1, cp):  # (f3)
            if self.heavy(r + 1, j):
                return False
        return True

    def is_F_target(self, r: int, c: int) -> bool:
        return self.is_f_target(r, c) or self.is_fstar_target(r, c)

    def max_F_target(self):
        """Bottommost then rightmost heavy tile, or ``None`` for ``D_id``."""
        heavy = self.heavy_cells()
        if not heavy:
            return None
        r = max(i for (i, _j) in heavy)
        c = max(j for (i, j) in heavy if i == r)
        return (r, c)

    def is_F_terminal(self) -> bool:
        t = self.max_F_target()
        if t is None:
            return True
        return self.is_fstar_target(*t)

    def F_target_info(self, r: int, c: int) -> dict:
        fstar = self.is_fstar_target(r, c)
        assert fstar or self.is_f_target(r, c), f"({r},{c}) is not an F-target"
        kind = "fstar" if fstar else "f"
        if kind == "f":
            cp = next(j for j in range(c + 1, self.n + 1) if self.tile(r + 1, j) == "J")
        else:
            cp = max(j for j in range(1, self.n + 1) if self.tile(r, j) == "R")
        tc = self.tile(r, c)
        if tc == "B":
            b = c
        else:  # M
            b = max(j for j in range(1, c) if self.tile(r, j) == "R")
        # left trichotomy
        if tc == "B":
            left = "Blank"
        elif self.is_pipe_segment(r + 1, b, c):
            left = "Cross"
        else:
            left = "NonCross"
        # right trichotomy
        if kind == "fstar":
            right, rho = "Terminal", cp
        else:
            d = None
            for dd in range(c + 1, cp):
                if self.is_doublecross(r, dd, cp):
                    d = dd
                    break
            if d is not None:
                right, rho = "DCross", d
            else:
                right, rho = "Ordinary", cp
        return {"kind": kind, "cp": cp, "b": b, "left": left, "right": right, "rho": rho, "lam": c}

    def f_move(self, r: int, c: int) -> MBPD:
        r"""Apply the ``F``-move at the ``F``-target ``(r,c)`` (paper
        \S "do F move")."""
        info = self.F_target_info(r, c)
        D = self
        if D.tile(r, c) == "M":
            D = D.with_tile(r, c, "J")
        if info["left"] in ("Blank", "Cross"):
            D = D._swap_rows_core(r, c, info["rho"])  # (r,[c,rho])-undroop
        if D.tile(r + 1, info["cp"]) == "J":
            D = D.with_tile(r + 1, info["cp"], "M")
        return D

    # ---- E-targets and E-moves (inverse of F-moves) ---------------------
    def _e_target_cprime(self, r: int, c: int):
        r"""(e2): the maximum ``c' < c`` with ``D_{r,c'}=R`` such that
        ``D_{r+1,[c',c]}`` is *not* a pipe segment; then check (e3): no heavy
        tiles strictly right of ``(r,c')`` in row ``r``.  Returns ``c'`` or
        ``None`` if (e2)/(e3) fail."""
        cprime = None
        for j in range(c - 1, 0, -1):
            if self.tile(r, j) == "R" and not self.is_pipe_segment(r + 1, j, c):
                cprime = j
                break
        if cprime is None:
            return None
        for j in range(cprime + 1, self.n + 1):  # (e3)
            if self.heavy(r, j):
                return None
        return cprime

    def e_target(self, r: int):
        r"""Determine the ``E``-target ``(r+1,c)`` in the two-row window with top
        row ``r`` (``1 <= r < n``).  Returns a dict with keys
        ``kind`` (``'e'`` or ``'estar'``), ``c``, ``cprime`` or ``None`` if no
        ``E``-target exists for this row."""
        if not (1 <= r < self.n):
            return None
        # e-target: (e1) leftmost heavy tile on row r+1.
        heavy_cols = [j for j in range(1, self.n + 1) if self.heavy(r + 1, j)]
        if heavy_cols:
            c = min(heavy_cols)
            cprime = self._e_target_cprime(r, c)
            if cprime is not None:
                return {"kind": "e", "c": c, "cprime": cprime}
            return None
        # e*-target: (e*1) no heavy in row r+1; c = largest with D_{r,c}=R or
        # D_{r+1,c}=R.
        c = None
        for j in range(self.n, 0, -1):
            if self.tile(r, j) == "R" or self.tile(r + 1, j) == "R":
                c = j
                break
        if c is None:
            return None
        cprime = self._e_target_cprime(r, c)
        if cprime is None:
            return None
        return {"kind": "estar", "c": c, "cprime": cprime}

    def E_target_info(self, r: int) -> dict:
        r"""Full case analysis for the ``E``-move at row ``r``: the right
        trichotomy (Initial/Plus/NPlus) with ``rho``, and the left trichotomy
        (Straight/DCross/Leftturn) with ``lambda``."""
        t = self.e_target(r)
        assert t is not None, f"row {r} has no E-target"
        c, cprime = t["c"], t["cprime"]
        # right trichotomy + rho
        if t["kind"] == "estar":
            right, rho = "Initial", c
        elif self.tile(r, c) == "P":
            right = "Plus"
            rho = max(j for j in range(cprime, c) if self.tile(r + 1, j) == "R")
        else:
            right, rho = "NPlus", c
        # left trichotomy + lambda
        if not self.is_pipe_segment(r, cprime, c):
            left = "Leftturn"
            lam = min(j for j in range(cprime + 1, self.n + 1) if self.tile(r, j) == "J")
        else:
            dd = next((j for j in range(cprime + 1, self.n + 1) if self.tile(r + 1, j) == "J"), None)
            if dd is not None and self.is_doublecross(r, cprime, dd):
                left, lam = "DCross", dd
            else:
                left, lam = "Straight", cprime
        return {"kind": t["kind"], "c": c, "cprime": cprime, "right": right, "rho": rho, "left": left, "lam": lam}

    def e_move(self, r: int) -> MBPD:
        r"""Apply the ``E``-move ``E_r`` at row ``r`` (paper \S "the E moves"):
          * if ``(r+1,c)`` is ``M``, unmark it;
          * in Cases Straight and DCross do the ``(r,[lambda,rho])``-droop;
          * if ``(r,lambda)`` is ``J``, mark it."""
        info = self.E_target_info(r)
        c, lam, rho = info["c"], info["lam"], info["rho"]
        D = self
        if D.tile(r + 1, c) == "M":
            D = D.with_tile(r + 1, c, "J")
        if info["left"] in ("Straight", "DCross"):
            D = D._swap_rows_core(r, lam, rho)  # (r,[lambda,rho])-droop
        if D.tile(r, lam) == "J":
            D = D.with_tile(r, lam, "M")
        return D

    # ---- Phi : MBPD -> RCP (row-pop recursion) --------------------------
    def phi(self) -> RCP:
        r"""The bijection ``Phi`` of the paper (Theorem "T: main"), realized by
        the row-pop recursion:

          * ``Phi(D_id) = ()``
          * F-terminal:    ``Phi(D) = (i,i) . Phi(f*_i(D))``
          * F-nonterminal: ``Phi(D) = down( Phi(f_i(D)) )``

        where ``i`` is the row of the maximum ``F``-target and ``down`` replaces
        the first biletter ``(i',a)`` by ``(i'-1,a)`` (RRCP2)."""
        biletters = self._phi_biletters()
        return RCP(tuple(biletters), self.n)

    def _phi_biletters(self):
        t = self.max_F_target()
        if t is None:
            return []
        r, c = t
        if self.is_fstar_target(r, c):  # F-terminal
            rest = self.f_move(r, c)._phi_biletters()
            return [(r, r), *rest]
        # F-nonterminal
        rest = self.f_move(r, c)._phi_biletters()
        (i0, a0) = rest[0]
        return [(i0 - 1, a0)] + rest[1:]

    # ---- BPD bridge ------------------------------------------------------
    # MBPD tiles map onto ordinary BPD tiles identically; the extra mark ``M``
    # is a marked up-elbow (``J``) recorded separately.  This lets us reuse the
    # (unreduced-capable) ``BPD.resize`` machinery for zeroing out rows.
    _MBPD_TO_BPD: ClassVar[dict[str, str]] = {
        "B": "BLANK", "P": "CROSS", "H": "HORIZ",
        "V": "VERT", "J": "ELBOW_NW", "R": "ELBOW_SE", "M": "ELBOW_NW",
    }
    _BPD_TO_MBPD: ClassVar[dict[str, str]] = {
        "BLANK": "B", "CROSS": "P", "HORIZ": "H",
        "VERT": "V", "ELBOW_NW": "J", "ELBOW_SE": "R",
    }

    def to_bpd(self):
        """Return ``(bpd, marks)``: the underlying :class:`BPD` and the frozenset
        of ``(i, j)`` cells (1-indexed) that carry a mark (``M`` tiles)."""
        import numpy as np

        from schubmult.combinatorics.bpd import BPD, TileType

        grid = np.array(
            [[TileType[self._MBPD_TO_BPD[self.tile(i, j)]] for j in range(1, self.n + 1)] for i in range(1, self.n + 1)],
            dtype=TileType,
        )
        marks = frozenset((i, j) for i in range(1, self.n + 1) for j in range(1, self.n + 1) if self.tile(i, j) == "M")
        return BPD(grid), marks

    @classmethod
    def from_bpd(cls, bpd, marks=frozenset()) -> MBPD:
        """Rebuild an :class:`MBPD` from a :class:`BPD` and a set of marked
        ``(i, j)`` cells.  The BPD is squared up to ``max(rows, len(perm))`` first
        (Rothe completion), and a mark is only reinstated where the tile is a
        ``J`` up-elbow."""
        n = max(bpd.rows, len(bpd.perm))
        bpd = bpd.resize(n)
        grid = [[cls._BPD_TO_MBPD[bpd._grid[i, j].name] for j in range(n)] for i in range(n)]
        D = cls.from_tiles(grid)
        for (i, j) in marks:
            if i <= n and j <= n and D.tile(i, j) == "J":
                D = D.with_tile(i, j, "M")
        return D

    def zero_out_last_row(self, rows: int) -> MBPD:
        """Zero out below ``rows`` by resizing the underlying BPD down to ``rows``
        rows (Weigandt/Lascoux transition) and re-deducing the marks that survive.

        Weight is preserved; the permutation changes according to the transition."""
        bpd, marks = self.to_bpd()
        bpd = bpd.resize(rows)
        return MBPD.from_bpd(bpd, frozenset((i, j) for (i, j) in marks if i <= rows))

    def __str__(self):
        rows = []
        for i in range(1, self.n + 1):
            rows.append("".join(_GLYPH[self.tile(i, j)] for j in range(1, self.n + 1)))
        return "\n".join(rows)

    def __repr__(self):
        return f"MBPD(n={self.n})\n{self}"


# =====================================================================
#  Reverse compatible pairs (RCP)
# =====================================================================
@dataclass(frozen=True)
class RCP:
    r"""A reverse compatible pair: a tuple of biletters ``(i, a)`` with
    ``1 <= i <= a < n`` (``n`` is the ambient size), strictly *decreasing* in
    the order ``(i1,a1) > (i2,a2)`` iff ``i1 > i2`` or (``i1 == i2`` and
    ``a1 < a2``).
    """

    biletters: tuple  # tuple of (i, a)
    n: int

    @staticmethod
    def _gt(x, y) -> bool:
        (i1, a1), (i2, a2) = x, y
        return i1 > i2 or (i1 == i2 and a1 < a2)

    def is_valid(self) -> bool:
        for (i, a) in self.biletters:
            if not (1 <= i <= a < self.n):
                return False
        for k in range(len(self.biletters) - 1):
            if not self._gt(self.biletters[k], self.biletters[k + 1]):
                return False
        return True

    def weight(self):
        m = [0] * self.n
        for (i, _a) in self.biletters:
            m[i - 1] += 1
        return tuple(m)

    def __len__(self):
        return len(self.biletters)

    def perm(self):
        r"""``w(B) = s_{a_l} * ... * s_{a_1}`` (Demazure product, subscripts
        decreasing = reversed biletter order)."""
        from schubmult.combinatorics.permutation import Permutation

        perm = Permutation([])
        # a_l, a_{l-1}, ..., a_1  == reversed order of the biletter list.
        for (_i, a) in reversed(self.biletters):
            perm = perm @ Permutation.ref_product(a)
        return perm

    # ---- bridges ---------------------------------------------------------
    def to_wcgraph(self):
        r"""Package the RCP as a :class:`WCGraph`.

        ``WCGraph.from_word_compatible`` wants a *weakly increasing* compatible
        sequence.  Reversing the biletter list turns the RCP's decreasing order
        into weakly increasing ``i`` with strictly decreasing ``a`` within a
        block -- exactly the WCGraph compatibility condition.
        """
        from schubmult.combinatorics.wc_graph import WCGraph

        rev = list(reversed(self.biletters))
        seq = [i for (i, _a) in rev]
        word = [a for (_i, a) in rev]
        return WCGraph.from_word_compatible(word, seq, length=self.n)

    @classmethod
    def from_wcgraph(cls, wc, n=None):
        """Inverse of :meth:`to_wcgraph`."""
        word = wc.perm_word
        seq = wc.compatible_sequence
        rev = list(zip(seq, word))  # (i, a) in weakly-increasing-i order
        biletters = tuple(reversed(rev))  # back to RCP (decreasing) order
        if n is None:
            n = (max((a for (_i, a) in biletters), default=0) + 1) if biletters else 1
            n = max(n, max((i for (i, _a) in biletters), default=0) + 1, len(wc) + 1)
        return cls(biletters, n)

    def to_pd_crossings(self):
        r"""Remark "R:CP and PD": biletter ``(i, a)`` -> crossing at
        ``(i, a - i + 1)`` in the ordinary pipedream picture."""
        return [(i, a - i + 1) for (i, a) in self.biletters]

    # ---- Psi : RCP -> MBPD (row-unpop, inverse of Phi) ------------------
    def psi(self) -> MBPD:
        r"""The inverse bijection ``Psi = Phi^{-1}`` (Theorem "T: main"),
        realized by row-unpop.

        ``Phi(D) = (i,a) . Phi(nabla_r D)`` where the maximum ``F``-target of
        ``D`` is in row ``i`` and ``nabla_r D = f*_a f_{a-1} ... f_i (D)`` is
        obtained by climbing rows ``i, i+1, ..., a`` (the last, at row ``a``,
        being the ``f*``-move).  Inverting: given the first biletter ``(i,a)``,
        first rebuild ``N = nabla_r D = Psi(rest)``, then apply
        ``D = e_i(e_{i+1}( ... e_{a-1}( e*_a(N) ) ... ))``.
        """
        from schubmult.combinatorics.permutation import Permutation

        if len(self.biletters) == 0:
            ident = Permutation(list(range(1, self.n + 1)))
            return MBPD.rothe(ident, n=self.n)
        (i, a) = self.biletters[0]
        rest = RCP(self.biletters[1:], self.n)
        N = rest.psi()
        D = N.e_move(a)  # e*_a  (terminal inverse)
        for r in range(a - 1, i - 1, -1):
            D = D.e_move(r)
        return D

