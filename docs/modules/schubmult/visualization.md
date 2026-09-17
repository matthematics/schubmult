<a id="schubmult.visualization"></a>

# schubmult.visualization

Visualization tools for Schubert calculus objects.

This module provides functions to visualize combinatorial structures like
pipe dreams from RC graphs.

<a id="schubmult.visualization.draw_pipe_dream"></a>

#### draw\_pipe\_dream

```python
def draw_pipe_dream(rc,
                    max_size=None,
                    title=None,
                    ax=None,
                    flip_horizontal=True,
                    top_labeled=False,
                    show_refs=True)
```

Draw a pipe dream visualization of an RC graph.

In a pipe dream:
- Positions with elements: strands CROSS (go straight through)
- Empty positions: strands AVOID (make 90-degree elbow turns)

Label positioning depends on flip_horizontal and top_labeled:

- flip_horizontal=True, top_labeled=False (default):
Top shows column numbers, right shows permutation output
- flip_horizontal=True, top_labeled=True:
Top shows permutation output, right shows row numbers (1,2,3,...)
- flip_horizontal=False, top_labeled=False:
Left shows permutation output, top shows column numbers
- flip_horizontal=False, top_labeled=True:
Left shows row numbers (1,2,3,...), top shows permutation output

**Arguments**:

- `rc` - RCGraph object to visualize
- `max_size` - Maximum grid size to display (default: determined from permutation)
- `title` - Optional title for the plot
- `ax` - Optional matplotlib axes to draw on (creates new figure if None)
- `flip_horizontal` - If True, reflect horizontally (default: True)
- `top_labeled` - If True, swap which side shows input vs output labels (default: False)
- `show_refs` - If True, display reflection numbers at crossings in green (default: False)
  

**Returns**:

- `tuple` - (fig, ax) matplotlib figure and axes objects
  

**Examples**:

  >>> from schubmult import Permutation, RCGraph
  >>> from schubmult.visualization import draw_pipe_dream
  >>> perm = Permutation([2, 1, 3])
  >>> rc = list(RCGraph.all_rc_graphs(perm, 2))[0]
  >>> fig, ax = draw_pipe_dream(rc, title=f"Pipe Dream for {perm}")
  >>> plt.show()

<a id="schubmult.visualization.draw_pipe_dream_tikz"></a>

#### draw\_pipe\_dream\_tikz

```python
def draw_pipe_dream_tikz(rc,
                         max_size=None,
                         flip_horizontal=True,
                         top_labeled=False,
                         show_refs=False,
                         scale=1.0,
                         outline_rows=None,
                         clip_at_outline=True)
```

Generate TikZ code for a pipe dream visualization of an RC graph.

**Arguments**:

- `rc` - RCGraph object to visualize
- `max_size` - Maximum grid size to display (default: determined from permutation)
- `flip_horizontal` - If True, reflect horizontally (default: True)
- `top_labeled` - If True, swap which side shows input vs output labels (default: False)
- `show_refs` - If True, display reflection numbers at crossings (default: True)
- `scale` - Scale factor for the TikZ picture (default: 1.0)
- `outline_rows` - If provided, draw a thick black outline around this many rows from bottom (default: None)
- `clip_at_outline` - If True and outline_rows is set, clip strands at outline_rows + 1 (default: False)
  

**Returns**:

- `str` - TikZ code as a string
  

**Examples**:

  >>> from schubmult import Permutation, RCGraph
  >>> from schubmult.visualization import draw_pipe_dream_tikz
  >>> perm = Permutation([2, 1, 3])
  >>> rc = list(RCGraph.all_rc_graphs(perm, 2))[0]
  >>> tikz_code = draw_pipe_dream_tikz(rc)
  >>> print(tikz_code)

