from schubmult.utils._grid_print import GridPrint


class ChuteMoveElement(GridPrint):
    def __init__(self, rc_graph, rows):
        super().__init__()
        self._rc = rc_graph
        self._rows = set(rows)
        if len(self._rows) == 0:
            self._pair = (self._rc, self._rc)
            return
        new_rc = self._rc
        if any(abs(r1 - r2) <= 1 for r1 in self._rows for r2 in self._rows if r1 != r2):
            raise ValueError(f"Invalid chute move: rows {self._rows} are not separated by at least one row")

        for r in rows:
            row = r
            if row >= len(self._rc) - 1:
                raise ValueError(f"Invalid chute move: row {row} is the last row of the RC graph")
            found = False
            for index in range(len(self._rc[row - 1])):
                col = self._rc[row - 1][index] - row + 1
                if not self._rc.has_element(row + 1, col):
                    col2 = col - 1
                    while col2 > 0 and self._rc.has_element(row, col2) and not self._rc.has_element(row+1, col2):
                        col2 -= 1
                    initial_col2 = col2 + 1
                    while col2 > 0 and self._rc.has_element(row, col2) and self._rc.has_element(row + 1, col2):
                        col2 -= 1
                    if col2 == 0:
                        raise ValueError(f"Invalid chute move: no empty space below row {row} in column {col}")
                    # col2 += 1
                    new_rc = new_rc.toggle_ref_at(row, initial_col2).toggle_ref_at(row+1, col2)
                    found = True
                    break
            if not found:
                raise ValueError(f"Invalid chute move: no empty space below row {row} in any column")
        self._pair = (self._rc, new_rc)

    @property
    def chute_degree(self):
        return len(self._rows)

    @property
    def chute_move_rows(self):
        return self._rows

    @property
    def cols(self):
        return self._rc.cols

    @property
    def rows(self):
        return len(self._rc)

    @property
    def print_element(self):
        class PrintElement(GridPrint):

            _display_name = self.__class__.__name__

            @property
            def rows(inner_self) -> int:
                return self.rows

            @property
            def cols(inner_self) -> int:
                return self.cols

            @property
            def grey_background_rows(inner_self) -> set:
                """Returns set of 0-indexed row numbers that should have grey background padding."""
                grey_rows = set()
                for row_idx in range(self.rows):
                    # Row is amorphous if it's not in self._rows (1-indexed) and not row below a chute move
                    if row_idx + 1 not in self._rows and row_idx not in self._rows:
                        grey_rows.add(row_idx)
                return grey_rows

            def __getitem__(inner_self, key):
                if isinstance(key, tuple):
                    row, col = key
                    # if row + 1 not in self._rows and row not in self._rows:
                    #     return "\\e[9m" + str(self._rc[row, col]) + "\\e[0m"
                    #item1, item2 = self._pair[0][key], self._pair[1][key]
                    if row + 1 in self._rows:
                        if not self._rc.has_element(row + 1, col + 1):
                            return " "
                        if self._rc.has_element(row + 1, col + 1) and not self._pair[1].has_element(row + 1, col + 1):
                            return f"\033[48;5;242m{col + row + 1}\033[0m"
                        return col + row + 1
                    if (row - 1) + 1 in self._rows:
                        if not self._pair[1].has_element(row + 1, col + 1):
                            return " "
                        if self._pair[1].has_element(row + 1, col + 1) and not self._rc.has_element(row + 1, col + 1):
                            return f"\033[48;5;242m{col + row + 1}\033[0m"
                        return col + row + 1
                    # Amorphous rows - grey out everything
                    if self._rc.has_element(row + 1, col + 1):
                        return f"\033[48;5;242m{row + col + 1}\033[0m"
                    return "\033[48;5;242m \033[0m"
                if isinstance(key, int):
                    return self._rc[key]
                return "FATBACON"
        return PrintElement()

    def __hash__(self):
        return hash(tuple([self._rc.perm] + [len(self._rc)] + [(self._rc[r - 1], self._rc[r]) for r in self._rows]))

    def __eq__(self, other):
        if not isinstance(other, ChuteMoveElement):
            return False
        if self._rc.perm != other._rc.perm or self._rows != other._rows:
            return False
        if len(self._rc) != len(other._rc):
            return False
        for r in self._rows:
            if self._rc[r - 1] != other._rc[r - 1]:
                return False
            if self._rc[r] != other._rc[r]:
                return False
        return True
