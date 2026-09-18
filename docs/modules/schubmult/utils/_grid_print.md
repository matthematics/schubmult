<a id="schubmult.utils._grid_print"></a>

# schubmult.utils.\_grid\_print

`GridPrint`: a SymPy ``Printable`` mixin that renders a 2-D grid (``rows``, ``cols``, ``self[i, j]``)
as an aligned table for str/pretty/LaTeX output; used by RC graphs, BPDs and similar diagrams.

<a id="schubmult.utils._grid_print.GridPrint"></a>

## GridPrint Objects

```python
class GridPrint(Printable)
```

Mixin: subclasses provide ``rows``, ``cols``, ``__getitem__((i, j))`` and ``_display_name``.

