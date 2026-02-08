from sympy.printing.defaults import Printable


class GridPrint(Printable):

    @property
    def rows(self) -> int: ...

    @property
    def cols(self) -> int: ...

    def __getitem__(self, key): ...

    def _format_str(self, printer=None) -> str:
        from sympy.printing import StrPrinter

        if not printer:
            printer = StrPrinter()
        # Handle zero dimensions:
        if self.rows == 1:
            return self._display_name + "([%s])" % self.table(printer, rowsep=",\n")  # noqa: UP031
        return self._display_name + "([\n%s])" % self.table(printer, rowsep=",\n")  # noqa: UP031

    def table(self, printer, rowstart="|", rowend="|", rowsep="\n", colsep=" ", align="right"):
        table: list[list[str]] = []
        # Track per-column max lengths for pretty alignment
        maxlen = [0] * self.cols
        for i in range(self.rows):
            table.append([])
            for j in range(self.cols):
                s = printer._print(self[i, j]) if self[i, j] is not None else " "
                table[-1].append(s)
                maxlen[j] = max(len(s), maxlen[j])
        # Patch strings together
        align = {
            "left": "ljust",
            "right": "rjust",
            "center": "center",
            "<": "ljust",
            ">": "rjust",
            "^": "center",
        }[align]
        res = [""] * len(table)
        for i, row in enumerate(table):
            for j, elem in enumerate(row):
                row[j] = getattr(elem, align)(maxlen[j])
            res[i] = rowstart + colsep.join(row) + rowend
        return rowsep.join(res)

    def _sympystr(self, printer=None):
        from sympy.printing.str import StrPrinter

        if not printer:
            printer = StrPrinter()
        return printer._print_MatrixBase(self.print_element)

    def __repr__(self):
        return self._sympyrepr()

    def _sympyrepr(self):
        return self._sympystr()

    def _pretty(self, printer):
        return printer._print_MatrixBase(self.print_element)

    def _latex(self, printer):
        return printer._print_MatrixBase(self.print_element)

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

            def __getitem__(inner_self, key):
                item = self[key]
                if item is None:
                    return " "
                if isinstance(item, tuple):
                    return tuple([i if i is not None else " " for i in item])
                return item
        return PrintElement()
