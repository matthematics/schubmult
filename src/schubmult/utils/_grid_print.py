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
        import re
        table: list[list[str]] = []
        # Track per-column max lengths for pretty alignment
        maxlen = [0] * self.cols
        ansi_escape = re.compile(r'\x1b\[[0-9;]*m')
        for i in range(self.rows):
            table.append([])
            for j in range(self.cols):
                s = printer._print(self[i, j]) if self[i, j] is not None else " "
                table[-1].append(s)
                # Calculate visual length by removing ANSI escape codes
                visual_len = len(ansi_escape.sub('', s))
                maxlen[j] = max(visual_len, maxlen[j])
        # Patch strings together
        align = {
            "left": "ljust",
            "right": "rjust",
            "center": "center",
            "<": "ljust",
            ">": "rjust",
            "^": "center",
        }[align]

        # Check if this object has grey_background_rows property
        grey_rows = getattr(self, 'grey_background_rows', None)
        if grey_rows is None:
            grey_rows = set()

        res = [""] * len(table)
        for i, row in enumerate(table):
            for j, elem in enumerate(row):
                # Calculate how much padding is needed based on visual length
                visual_len = len(ansi_escape.sub('', elem))
                padding_needed = maxlen[j] - visual_len
                # Check if this row should have grey background padding
                should_grey_pad = i in grey_rows

                # If this is a grey row, wrap everything in one continuous grey background
                if should_grey_pad:
                    # Strip any existing ANSI codes to get clean content
                    clean_elem = ansi_escape.sub('', elem)
                    # Apply padding first
                    if align == "ljust":
                        padded_content = clean_elem + " " * padding_needed
                    elif align == "rjust":
                        padded_content = " " * padding_needed + clean_elem
                    else:  # center
                        left_pad = padding_needed // 2
                        right_pad = padding_needed - left_pad
                        padded_content = " " * left_pad + clean_elem + " " * right_pad
                    # Wrap entire cell in grey background
                    row[j] = '\033[48;5;242m' + padded_content + '\033[0m'
                else:
                    # Normal padding without grey background
                    if align == "ljust":
                        row[j] = elem + " " * padding_needed
                    elif align == "rjust":
                        row[j] = " " * padding_needed + elem
                    else:  # center
                        left_pad = padding_needed // 2
                        right_pad = padding_needed - left_pad
                        row[j] = " " * left_pad + elem + " " * right_pad
            # Join columns with appropriate separator
            if i in grey_rows:
                # For grey rows, make colsep also grey
                grey_colsep = '\033[48;5;242m' + colsep + '\033[0m'
                res[i] = rowstart + grey_colsep.join(row) + rowend
            else:
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
