"""
Visualization tools for Schubert calculus objects.

This module provides functions to visualize combinatorial structures like
pipe dreams from RC graphs.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Arc


def draw_pipe_dream(rc, max_size=None, title=None, ax=None):
    """
    Draw a pipe dream visualization of an RC graph.

    In a pipe dream:
    - Positions with elements: strands CROSS (go straight through)
    - Empty positions: strands AVOID (make 90-degree elbow turns)

    Strands start at the left border (indices 1, 2, ...) and end at the top
    border showing the permutation values (which row ends up in each column).

    Args:
        rc: RCGraph object to visualize
        max_size: Maximum grid size to display (default: determined from permutation)
        title: Optional title for the plot
        ax: Optional matplotlib axes to draw on (creates new figure if None)

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
    perm_size = len(perm) if len(perm) > 0 else 2

    # Use perm_size as the grid size if not specified
    if max_size is None:
        max_size = perm_size

    # Create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(max_size * 1.2, max_size * 1.2))
    else:
        fig = ax.get_figure()

    # Track which positions have elements (crossings)
    crossings = set()
    for row_idx, row in enumerate(rc, start=1):
        for element in row:
            col = element - row_idx + 1
            if 1 <= col <= max_size:
                crossings.add((row_idx, col))

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

            if (row, col) in crossings:
                # CROSSING: strands go straight through
                # Horizontal strand (blue, from left)
                ax.plot([col - 1, col], [y_pos - 0.5, y_pos - 0.5], "b-", linewidth=2.5, zorder=10, solid_capstyle="round")
                # Vertical strand (red, from bottom)
                ax.plot([col - 0.5, col - 0.5], [y_pos - 1, y_pos], "r-", linewidth=2.5, zorder=10, solid_capstyle="round")
            else:
                # ELBOW: strands avoid each other with 90-degree curves
                # Blue arc: from left to up, centered at lower-left corner
                arc1 = Arc((col - 1, y_pos), 1, 1, angle=0, theta1=270, theta2=360, color="blue", linewidth=2.5, zorder=10)
                ax.add_patch(arc1)

                # Red arc: from bottom to right, centered at upper-right corner
                # Only draw if NOT on the diagonal (boundary cells only need upper arc)
                if not on_diagonal:
                    arc2 = Arc((col, y_pos - 1), 1, 1, angle=0, theta1=90, theta2=180, color="red", linewidth=2.5, zorder=10)
                    ax.add_patch(arc2)

    # Add labels for initial positions (left border)
    for i in range(1, max_size + 1):
        y_pos = max_size - i + 0.5
        ax.text(-0.3, y_pos, str(i), ha="right", va="center", fontsize=12, fontweight="bold", color="blue")

    # Add labels for final positions (top border) - show which row ends up in each column
    for col in range(1, max_size + 1):
        # Find which row's strand ends up in this column by looking at the permutation
        # perm[i-1] gives the column where row i ends up, so we need the inverse
        row_label = None
        for row in range(1, max_size + 1):
            if perm[row - 1] == col:
                row_label = str(row)
                break
        if row_label:
            ax.text(col - 0.5, max_size + 0.3, row_label, ha="center", va="bottom", fontsize=12, fontweight="bold", color="red")

    # Title
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold", pad=20)

    ax.set_xlim(-0.8, max_size + 0.5)
    ax.set_ylim(-0.5, max_size + 0.8)
    ax.set_aspect("equal")
    ax.axis("off")

    plt.tight_layout()
    return fig, ax
