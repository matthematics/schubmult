"""
Visualization tools for Schubert calculus objects.

This module provides functions to visualize combinatorial structures like
pipe dreams from RC graphs.
"""

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.patches import Arc


def draw_pipe_dream(rc, max_size=None, title=None, ax=None, flip_horizontal=True, top_labeled=False, show_refs=True):
    """
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

    Args:
        rc: RCGraph object to visualize
        max_size: Maximum grid size to display (default: determined from permutation)
        title: Optional title for the plot
        ax: Optional matplotlib axes to draw on (creates new figure if None)
        flip_horizontal: If True, reflect horizontally (default: True)
        top_labeled: If True, swap which side shows input vs output labels (default: False)
        show_refs: If True, display reflection numbers at crossings in green (default: False)

    Returns:
        tuple: (fig, ax) matplotlib figure and axes objects

    Examples:
        >>> from schubmult import Permutation, RCGraph
        >>> from schubmult.visualization import draw_pipe_dream
        >>> perm = Permutation([2, 1, 3])
        >>> rc = list(RCGraph.all_rc_graphs(perm, 2))[0]
        >>> fig, ax = draw_pipe_dream(rc, title=f"Pipe Dream for {perm}")
        >>> plt.show()
    """
    # Get the permutation size from the RC graph's permutation
    perm = rc.perm
    perm_size = len(perm) if len(perm) > 0 else 0

    # Calculate the actual grid size needed from RC graph elements (for unreduced graphs)
    if max_size is None:
        max_size = perm_size if perm_size > 0 else 1
        # Check all elements in the RC graph to ensure grid is large enough
        for row_idx, row in enumerate(rc, start=1):
            if len(row) > 0:
                max_element = max(row)
                max_size = max(max_size, max_element, row_idx)

    # Create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(max_size * 1.2, max_size * 1.2))
    else:
        fig = ax.get_figure()

    # Track which positions have elements (crossings) and their reflection values
    crossings = set()
    crossing_values = {}  # Maps (row, col) to reflection value
    for row_idx, row in enumerate(rc, start=1):
        for element in row:
            col = element - row_idx + 1
            if 1 <= col <= max_size and row_idx <= max_size:
                crossings.add((row_idx, col))
                crossing_values[(row_idx, col)] = element

    # Draw grid lines (light)
    for i in range(max_size + 1):
        ax.axhline(i, color="lightgray", linewidth=0.3, alpha=0.3)
        ax.axvline(i, color="lightgray", linewidth=0.3, alpha=0.3)

    # Draw all pipe segments - only within the actual grid
    for row in range(1, max_size + 1):
        for col in range(1, max_size + 1):
            # Skip cells beyond the anti-diagonal
            if row + col > perm_size + 1:
                continue

            # Check if we're on the anti-diagonal (boundary cells)
            on_diagonal = row + col == perm_size + 1

            y_pos = max_size - row + 1

            if flip_horizontal:
                # Reflect x-coordinates horizontally
                x_pos = max_size - col + 1
            else:
                x_pos = col

            if (row, col) in crossings:
                # CROSSING: strands go straight through
                if flip_horizontal:
                    # Horizontal strand (blue, from right to left)
                    ax.plot([x_pos, x_pos - 1], [y_pos - 0.5, y_pos - 0.5], "b-", linewidth=2.5, zorder=10, solid_capstyle="round")
                else:
                    # Horizontal strand (blue, from left to right)
                    ax.plot([x_pos - 1, x_pos], [y_pos - 0.5, y_pos - 0.5], "b-", linewidth=2.5, zorder=10, solid_capstyle="round")
                # Vertical strand (red, from bottom) - same in both orientations
                ax.plot([x_pos - 0.5, x_pos - 0.5], [y_pos - 1, y_pos], "r-", linewidth=2.5, zorder=10, solid_capstyle="round")
            else:
                # ELBOW: strands avoid each other with 90-degree curves
                if flip_horizontal:
                    # Blue arc: from right to up, centered at lower-right corner
                    arc1 = Arc((x_pos, y_pos), 1, 1, angle=0, theta1=180, theta2=270, color="blue", linewidth=2.5, zorder=10)
                    ax.add_patch(arc1)

                    # Red arc: from bottom to left, centered at upper-left corner
                    if not on_diagonal:
                        arc2 = Arc((x_pos - 1, y_pos - 1), 1, 1, angle=0, theta1=0, theta2=90, color="red", linewidth=2.5, zorder=10)
                        ax.add_patch(arc2)
                else:
                    # Blue arc: from left to up, centered at lower-left corner
                    arc1 = Arc((x_pos - 1, y_pos), 1, 1, angle=0, theta1=270, theta2=360, color="blue", linewidth=2.5, zorder=10)
                    ax.add_patch(arc1)

                    # Red arc: from bottom to right, centered at upper-right corner
                    if not on_diagonal:
                        arc2 = Arc((x_pos, y_pos - 1), 1, 1, angle=0, theta1=90, theta2=180, color="red", linewidth=2.5, zorder=10)
                        ax.add_patch(arc2)

    # Draw reflection numbers at crossings if requested
    if show_refs:
        # Scale fontsize with grid size for consistency
        ref_fontsize = max(8, min(20, max_size * 2.5))
        for row, col in crossings:
            y_pos = max_size - row + 1
            if flip_horizontal:
                x_pos = max_size - col + 1
            else:
                x_pos = col

            # Display the reflection value (simple reflection number)
            ref_value = crossing_values[(row, col)]
            text = ax.text(x_pos - 0.5, y_pos - 0.5, str(ref_value), ha="center", va="center", fontsize=ref_fontsize, fontweight="bold", color="forestgreen", zorder=20)
            text.set_path_effects([path_effects.Stroke(linewidth=3, foreground="white"), path_effects.Normal()])

    if flip_horizontal:
        if top_labeled:
            # Top shows permutation output, right shows row numbers
            for col in range(1, max_size + 1):
                # Find which row ends up in this column
                row_label = None
                for row in range(1, max_size + 1):
                    if perm[row - 1] == col:
                        row_label = str(row)
                        break
                if row_label:
                    x_pos = max_size - col + 0.5
                    ax.text(x_pos, max_size + 0.3, row_label, ha="center", va="bottom", fontsize=12, fontweight="bold", color="red")

            # Right border shows row numbers (input)
            for row in range(1, max_size + 1):
                y_pos = max_size - row + 0.5
                ax.text(max_size + 0.3, y_pos, str(row), ha="left", va="center", fontsize=12, fontweight="bold", color="blue")
        else:
            # Add labels for column positions (top border) - fixed labels counting down
            for col in range(1, max_size + 1):
                # Reflect column position: rightmost is column 1, leftmost is column max_size
                x_pos = max_size - col + 0.5
                ax.text(x_pos, max_size + 0.3, str(col), ha="center", va="bottom", fontsize=12, fontweight="bold", color="blue")

            # Add labels for output positions (right border) - show permutation values
            for row in range(1, max_size + 1):
                y_pos = max_size - row + 0.5
                # Row i enters from right and exits at column perm[i-1]
                output_col = perm[row - 1]
                ax.text(max_size + 0.3, y_pos, str(output_col), ha="left", va="center", fontsize=12, fontweight="bold", color="red")

        ax.set_xlim(-0.5, max_size + 0.8)
    else:
        if top_labeled:
            # Top shows permutation output (which input row ends up in each output column), left shows row numbers
            # Need inverse permutation: for each output column, find which input row goes there
            inverse_perm = [0] * max_size
            for row in range(1, max_size + 1):
                output_col = perm[row - 1]
                inverse_perm[output_col - 1] = row

            for col in range(1, max_size + 1):
                x_pos = col - 0.5
                ax.text(x_pos, max_size + 0.3, str(inverse_perm[col - 1]), ha="center", va="bottom", fontsize=12, fontweight="bold", color="red")

            # Left border shows row numbers (input)
            for row in range(1, max_size + 1):
                y_pos = max_size - row + 0.5
                ax.text(-0.3, y_pos, str(row), ha="right", va="center", fontsize=12, fontweight="bold", color="blue")
        else:
            # Add labels for output positions (left border) - show permutation values (same vertical as flipped mode)
            for row in range(1, max_size + 1):
                y_pos = max_size - row + 0.5
                # Row i outputs to column perm[i-1]
                output_col = perm[row - 1]
                ax.text(-0.3, y_pos, str(output_col), ha="right", va="center", fontsize=12, fontweight="bold", color="red")

            # Add labels for column positions (top border) - show column numbers in increasing order
            for col in range(1, max_size + 1):
                # Normal position: leftmost is 1, rightmost is max_size
                ax.text(col - 0.5, max_size + 0.3, str(col), ha="center", va="bottom", fontsize=12, fontweight="bold", color="blue")

        ax.set_xlim(-0.8, max_size + 0.5)

    # Title
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold", pad=20)

    ax.set_ylim(-0.5, max_size + 0.8)
    ax.set_aspect("equal")
    ax.axis("off")

    plt.tight_layout()
    return fig, ax


def draw_pipe_dream_tikz(rc, max_size=None, flip_horizontal=True, top_labeled=False, show_refs=True, scale=1.0, outline_rows=None, clip_at_outline=True):
    """
    Generate TikZ code for a pipe dream visualization of an RC graph.

    Args:
        rc: RCGraph object to visualize
        max_size: Maximum grid size to display (default: determined from permutation)
        flip_horizontal: If True, reflect horizontally (default: True)
        top_labeled: If True, swap which side shows input vs output labels (default: False)
        show_refs: If True, display reflection numbers at crossings (default: True)
        scale: Scale factor for the TikZ picture (default: 1.0)
        outline_rows: If provided, draw a thick black outline around this many rows from bottom (default: None)
        clip_at_outline: If True and outline_rows is set, clip strands at outline_rows + 1 (default: False)

    Returns:
        str: TikZ code as a string

    Examples:
        >>> from schubmult import Permutation, RCGraph
        >>> from schubmult.visualization import draw_pipe_dream_tikz
        >>> perm = Permutation([2, 1, 3])
        >>> rc = list(RCGraph.all_rc_graphs(perm, 2))[0]
        >>> tikz_code = draw_pipe_dream_tikz(rc)
        >>> print(tikz_code)
    """
    # Get the permutation size from the RC graph's permutation
    perm = rc.perm
    perm_size = len(perm) if len(perm) > 0 else 0

    # Calculate the actual grid size needed from RC graph elements (for unreduced graphs)
    if max_size is None:
        max_size = perm_size if perm_size > 0 else 1
        # Check all elements in the RC graph to ensure grid is large enough
        for row_idx, row in enumerate(rc, start=1):
            if len(row) > 0:
                max_element = max(row)
                max_size = max(max_size, max_element, row_idx)

    # Track which positions have elements (crossings) and their reflection values
    crossings = set()
    crossing_values = {}
    for row_idx, row in enumerate(rc, start=1):
        for element in row:
            col = element - row_idx + 1
            if 1 <= col <= max_size and row_idx <= max_size:
                crossings.add((row_idx, col))
                crossing_values[(row_idx, col)] = element

    # Start building TikZ code
    lines = []
    lines.append("\\begin{tikzpicture}[scale=" + str(scale) + "]")

    # Scale line width with the scale parameter
    line_width = 2.5 * scale

    # Draw grid lines (light gray)
    if clip_at_outline and outline_rows is not None:
        # Only draw grid in the clipped region
        y_min = max_size - outline_rows
        for i in range(outline_rows + 1):
            y = y_min + i
            lines.append(f"  \\draw[lightgray, very thin, opacity=0.3] (0,{y}) -- ({max_size},{y});")
        for i in range(max_size + 1):
            lines.append(f"  \\draw[lightgray, very thin, opacity=0.3] ({i},{y_min}) -- ({i},{max_size});")
    else:
        # Draw full grid
        for i in range(max_size + 1):
            lines.append(f"  \\draw[lightgray, very thin, opacity=0.3] (0,{i}) -- ({max_size},{i});")
            lines.append(f"  \\draw[lightgray, very thin, opacity=0.3] ({i},0) -- ({i},{max_size});")

    # Draw all pipe segments
    row_limit = outline_rows + 1 if (outline_rows is not None and clip_at_outline) else max_size + 1
    for row in range(1, row_limit):
        for col in range(1, max_size + 1):
            # Skip cells beyond the anti-diagonal
            if row + col > max_size + 1:
                continue

            # Check if we're on the anti-diagonal (boundary cells)
            on_diagonal = row + col == max_size + 1

            # TikZ coordinates (y increases upward, so we flip row)
            y_pos = max_size - row

            if flip_horizontal:
                x_pos = max_size - col
            else:
                x_pos = col - 1

            if (row, col) in crossings:
                # CROSSING: strands go straight through
                if flip_horizontal:
                    # Horizontal strand (blue, from right to left)
                    lines.append(f"  \\draw[blue, line width={line_width}pt, line cap=round] ({x_pos + 1},{y_pos + 0.5}) -- ({x_pos},{y_pos + 0.5});")
                else:
                    # Horizontal strand (blue, from left to right)
                    lines.append(f"  \\draw[blue, line width={line_width}pt, line cap=round] ({x_pos},{y_pos + 0.5}) -- ({x_pos + 1},{y_pos + 0.5});")
                # Vertical strand (red, from bottom to top)
                lines.append(f"  \\draw[red, line width={line_width}pt, line cap=round] ({x_pos + 0.5},{y_pos}) -- ({x_pos + 0.5},{y_pos + 1});")
            else:
                # ELBOW: strands avoid each other with 90-degree curves
                if flip_horizontal:
                    # Blue arc: from right to up (starts horizontal from right, ends vertical going up)
                    # Control points: first keeps horizontal tangent, second keeps vertical tangent
                    lines.append(
                        f"  \\draw[blue, line width={line_width}pt] ({x_pos + 1},{y_pos + 0.5}) .. controls ({x_pos + 0.7},{y_pos + 0.5}) and ({x_pos + 0.5},{y_pos + 0.7}) .. ({x_pos + 0.5},{y_pos + 1});",
                    )
                    # Red arc: from bottom to left (starts vertical from bottom, ends horizontal going left)
                    if not on_diagonal:
                        lines.append(f"  \\draw[red, line width={line_width}pt] ({x_pos + 0.5},{y_pos}) .. controls ({x_pos + 0.5},{y_pos + 0.3}) and ({x_pos + 0.3},{y_pos + 0.5}) .. ({x_pos},{y_pos + 0.5});")
                else:
                    # Blue arc: from left to up (starts horizontal from left, ends vertical going up)
                    lines.append(
                        f"  \\draw[blue, line width={line_width}pt] ({x_pos},{y_pos + 0.5}) .. controls ({x_pos + 0.3},{y_pos + 0.5}) and ({x_pos + 0.5},{y_pos + 0.7}) .. ({x_pos + 0.5},{y_pos + 1});",
                    )
                    # Red arc: from bottom to right (starts vertical from bottom, ends horizontal going right)
                    if not on_diagonal:
                        lines.append(
                            f"  \\draw[red, line width={line_width}pt] ({x_pos + 0.5},{y_pos}) .. controls ({x_pos + 0.5},{y_pos + 0.3}) and ({x_pos + 0.7},{y_pos + 0.5}) .. ({x_pos + 1},{y_pos + 0.5});",
                        )

    # Draw reflection numbers at crossings if requested
    if show_refs:
        for row, col in crossings:
            y_pos = max_size - row
            if flip_horizontal:
                x_pos = max_size - col
            else:
                x_pos = col - 1

            ref_value = crossing_values[(row, col)]
            lines.append(f"  \\node[font=\\Large\\bfseries, text=green!60!black, fill=white, inner sep=0.5pt, circle, transform shape] at ({x_pos + 0.5},{y_pos + 0.5}) {{{ref_value}}};")

    # Add labels
    if flip_horizontal:
        if top_labeled:
            # Top shows permutation output
            for col in range(1, max_size + 1):
                row_label = None
                for row in range(1, max_size + 1):
                    if perm[row - 1] == col:
                        row_label = str(row)
                        break
                if row_label:
                    x_pos = max_size - col
                    lines.append(f"  \\node[font=\\bfseries, text=red, anchor=south] at ({x_pos + 0.5},{max_size + 0.3}) {{{row_label}}};")

            # Right shows row numbers
            row_range = range(1, outline_rows + 1) if (clip_at_outline and outline_rows is not None) else range(1, max_size + 1)
            for row in row_range:
                y_pos = max_size - row
                lines.append(f"  \\node[font=\\bfseries, text=blue, anchor=west] at ({max_size + 0.3},{y_pos + 0.5}) {{{row}}};")
        else:
            # Top shows column numbers
            for col in range(1, max_size + 1):
                x_pos = max_size - col
                lines.append(f"  \\node[font=\\bfseries, text=blue, anchor=south] at ({x_pos + 0.5},{max_size + 0.3}) {{{col}}};")

            # Right shows permutation output
            row_range = range(1, outline_rows + 1) if (clip_at_outline and outline_rows is not None) else range(1, max_size + 1)
            for row in row_range:
                y_pos = max_size - row
                output_col = perm[row - 1]
                lines.append(f"  \\node[font=\\bfseries, text=red, anchor=west] at ({max_size + 0.3},{y_pos + 0.5}) {{{output_col}}};")
    else:
        if top_labeled:
            # Top shows permutation output (inverse permutation)
            inverse_perm = [0] * max_size
            for row in range(1, max_size + 1):
                output_col = perm[row - 1]
                inverse_perm[output_col - 1] = row

            for col in range(1, max_size + 1):
                x_pos = col - 1
                lines.append(f"  \\node[font=\\bfseries, text=red, anchor=south] at ({x_pos + 0.5},{max_size + 0.3}) {{{inverse_perm[col - 1]}}};")

            # Left shows row numbers
            row_range = range(1, outline_rows + 1) if (clip_at_outline and outline_rows is not None) else range(1, max_size + 1)
            for row in row_range:
                y_pos = max_size - row
                lines.append(f"  \\node[font=\\bfseries, text=blue, anchor=east] at ({-0.3},{y_pos + 0.5}) {{{row}}};")
        else:
            # Left shows permutation output
            row_range = range(1, outline_rows + 1) if (clip_at_outline and outline_rows is not None) else range(1, max_size + 1)
            for row in row_range:
                y_pos = max_size - row
                output_col = perm[row - 1]
                lines.append(f"  \\node[font=\\bfseries, text=red, anchor=east] at ({-0.3},{y_pos + 0.5}) {{{output_col}}};")

            # Top shows column numbers
            for col in range(1, max_size + 1):
                x_pos = col - 1
                lines.append(f"  \\node[font=\\bfseries, text=blue, anchor=south] at ({x_pos + 0.5},{max_size + 0.3}) {{{col}}};")

    # Draw thick black outline if requested
    if outline_rows is not None:
        # Outline the TOP outline_rows rows (from y = max_size - outline_rows to y = max_size)
        outline_width = 3 * scale
        lines.append(f"  \\draw[black, line width={outline_width}pt] (0,{max_size - outline_rows}) rectangle ({max_size},{max_size});")

    lines.append("\\end{tikzpicture}")

    return "\n".join(lines)
