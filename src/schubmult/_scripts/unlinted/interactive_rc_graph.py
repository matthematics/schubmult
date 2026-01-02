#!/usr/bin/env python3
"""
Interactive RC Graph Editor

Click on any grid cell to toggle a crossing at that position.
The pipe dream visualization updates in real time.

Usage:
    python interactive_rc_graph.py [permutation]

Examples:
    python interactive_rc_graph.py 3 1 2
    python interactive_rc_graph.py 4 2 1 3
    python interactive_rc_graph.py  # defaults to [2, 1, 3]
"""

import sys
import os
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import TextBox, Button
from schubmult import Permutation, RCGraph
from schubmult.visualization import draw_pipe_dream, draw_pipe_dream_tikz


class InteractiveRCGraph:
    """Interactive RC graph editor with clickable crossings."""
    
    def __init__(self, initial_perm=None, max_size=None):
        """
        Initialize the interactive RC graph editor.
        
        Args:
            initial_perm: Starting permutation (defaults to [2, 1, 3])
            max_size: Maximum grid size (defaults to 4)
        """
        if initial_perm is None:
            initial_perm = Permutation([2, 1, 3])
        elif not isinstance(initial_perm, Permutation):
            initial_perm = Permutation(initial_perm)
            
        self.current_perm = initial_perm
        
        # Start with the first RC graph for this permutation
        all_graphs = list(RCGraph.all_rc_graphs(initial_perm))
        if all_graphs:
            self.rc_graph = all_graphs[0]
        else:
            # Create empty RC graph
            self.rc_graph = RCGraph([tuple() for _ in range(len(initial_perm))])
        
        # Set fixed grid size (default 4x4)
        if max_size is None:
            self.max_size = 4
        else:
            self.max_size = max_size
            
        # Create the figure and connect events
        self.fig, self.ax = plt.subplots(figsize=(10, 10))
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        
        # Draw initial state
        self.redraw()
        
    def get_cell_from_click(self, x, y):
        """
        Convert click coordinates to grid cell (row, col).
        
        Returns (row, col) in 1-indexed coordinates, or None if outside grid.
        """
        # The grid goes from 0 to max_size in both dimensions
        # We need to convert to 1-indexed row and col
        
        # y coordinate: higher y means lower row index
        # row 1 is at top (y = max_size to max_size-1)
        # row max_size is at bottom (y = 1 to 0)
        row = self.max_size - int(y)
        
        # x coordinate: column increases with x
        col = int(x) + 1
        
        # Check if within valid grid
        if row < 1 or row > self.max_size or col < 1 or col > self.max_size:
            return None
            
        return (row, col)
    
    def toggle_crossing(self, row, col):
        """
        Toggle the crossing at position (row, col).
        
        Args:
            row: 1-indexed row
            col: 1-indexed column
        """
        try:
            # Extend RC graph if necessary to accommodate the new crossing
            current_length = len(self.rc_graph)
            if row > current_length:
                # Need to extend the RC graph to have enough rows
                self.rc_graph = self.rc_graph.extend(row - current_length)
            
            # Toggle the reflection at this position
            self.rc_graph = self.rc_graph.toggle_ref_at(row, col)
            
            # Update current permutation based on new RC graph
            self.current_perm = self.rc_graph.perm
            
            # Resize display if permutation grew
            self.update_display_size()
            
            return True
        except Exception as e:
            print(f"Error toggling crossing at ({row}, {col}): {e}")
            import traceback
            traceback.print_exc()
            return False
    
    def update_display_size(self):
        """Update the display size based on the current RC graph and permutation."""
        # Calculate required size from permutation
        perm_size = len(self.current_perm) if len(self.current_perm) > 0 else 1
        new_size = max(perm_size, 4)  # At least 4x4
        
        # Check all elements in the RC graph
        for row_idx, row in enumerate(self.rc_graph, start=1):
            if len(row) > 0:
                max_element = max(row)
                new_size = max(new_size, max_element, row_idx)
        
        # Update if size changed
        if new_size != self.max_size:
            self.max_size = new_size
            print(f"Display resized to {self.max_size}x{self.max_size}")
    
    def on_click(self, event):
        """Handle mouse click events."""
        if event.inaxes != self.ax:
            return
            
        # Get the cell coordinates
        cell = self.get_cell_from_click(event.xdata, event.ydata)
        if cell is None:
            return
            
        row, col = cell
        
        # Toggle the crossing
        if self.toggle_crossing(row, col):
            # Redraw if successful
            self.redraw()
            print(f"Toggled crossing at row={row}, col={col}")
            print(f"New permutation: {self.current_perm}")
            print(f"RC graph: {self.rc_graph}")
    
    def on_key_press(self, event):
        """Handle keyboard events."""
        if event.key == 's':
            self.save_to_latex()
        elif event.key == 'h':
            self.print_help()
    
    def print_help(self):
        """Print keyboard shortcuts."""
        print("\n=== Keyboard Shortcuts ===")
        print("s - Save current RC graph to LaTeX/TikZ document")
        print("h - Show this help message")
        print("==========================\n")
    
    def save_to_latex(self):
        """Export the current RC graph to a LaTeX document."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"rc_graph_{self.current_perm}_{timestamp}.tex"
        
        # Generate TikZ code for the pipe dream
        tikz_code = draw_pipe_dream_tikz(
            self.rc_graph,
            max_size=self.max_size,
            flip_horizontal=False,
            top_labeled=True,
            show_refs=True,
            scale=0.8
        )
        
        # Create a complete LaTeX document
        latex_document = f"""\\documentclass[border=2mm]{{standalone}}
\\usepackage{{tikz}}
\\usetikzlibrary{{arrows.meta}}

\\begin{{document}}

% RC Graph for permutation {self.current_perm}
% Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

{tikz_code}

\\end{{document}}
"""
        
        # Save to file
        try:
            with open(filename, 'w') as f:
                f.write(latex_document)
            print(f"\n✓ Saved LaTeX document to: {filename}")
            print(f"  Compile with: pdflatex {filename}")
            print(f"  Current permutation: {self.current_perm}")
            print(f"  RC graph: {self.rc_graph}\n")
            return filename
        except Exception as e:
            print(f"\n✗ Error saving LaTeX document: {e}\n")
            return None
    
    def redraw(self):
        """Redraw the entire visualization."""
        self.ax.clear()
        
        # Update figure size to match grid size
        self.fig.set_size_inches(self.max_size * 1.2, self.max_size * 1.2)
        
        # Draw the pipe dream
        draw_pipe_dream(
            self.rc_graph,
            max_size=self.max_size,
            title=f"RC Graph for {self.current_perm}\n(Click to toggle | Press 's' to save LaTeX | 'h' for help)",
            ax=self.ax,
            flip_horizontal=False,
            top_labeled=True,
            show_refs=True
        )
        
        # Add grid highlighting to make all cells visible and clickable
        for row in range(1, self.max_size + 1):
            for col in range(1, self.max_size + 1):
                # Add a very subtle rectangle to show clickable area
                y_pos = self.max_size - row
                x_pos = col - 1
                
                rect = Rectangle(
                    (x_pos, y_pos), 1, 1,
                    fill=False,
                    edgecolor='lightblue',
                    linewidth=0.5,
                    alpha=0.3,
                    linestyle=':'
                )
                self.ax.add_patch(rect)
        
        # Important: use draw_idle() instead of draw() to avoid blocking
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
    
    def show(self):
        """Display the interactive plot."""
        plt.show()


def main():
    """Main entry point for the interactive RC graph editor."""
    # Parse command line arguments
    if len(sys.argv) > 1:
        try:
            perm_list = [int(x) for x in sys.argv[1:]]
            initial_perm = Permutation(perm_list)
        except ValueError:
            print(f"Error: Invalid permutation. Please provide space-separated integers.")
            print(f"Example: python {sys.argv[0]} 3 1 2")
            return 1
    else:
        initial_perm = Permutation([2, 1, 3])
        print("Using default permutation [2, 1, 3]")
        print("Usage: python interactive_rc_graph.py [permutation]")
        print("Example: python interactive_rc_graph.py 3 1 2")
        print()
    
    # Create and show the interactive editor
    print(f"Starting interactive RC graph editor for permutation {initial_perm}")
    print("Click on any grid cell to toggle a crossing at that position.")
    print("Press 's' to save the current RC graph to a LaTeX document.")
    print("Press 'h' for help with keyboard shortcuts.")
    print()
    
    editor = InteractiveRCGraph(initial_perm)
    editor.show()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
