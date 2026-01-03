"""Generate LaTeX equation for RC graph multiplication."""

from schubmult import RCGraph
from schubmult.visualization import draw_pipe_dream_tikz


def generate_multiplication_equation(rc1, rc2, result_rcs, terms_per_line=3, scale=0.4):
    """
    Generate LaTeX code for an RC graph multiplication equation.
    
    Args:
        rc1: First RC graph (left multiplicand)
        rc2: Second RC graph (right multiplicand)
        result_rcs: List of RC graphs on the right side of equation
        terms_per_line: Number of terms to display per line (default: 3)
        scale: Scale factor for TikZ drawings (default: 0.4)
    
    Returns:
        str: LaTeX code for the equation
    """
    lines = []
    lines.append("\\begin{align*}")
    
    # Left side: rc1 * rc2 (2 rows outlined)
    lines.append("  &")
    lines.append(f"  \\vcenter{{\\hbox{{{draw_pipe_dream_tikz(rc1, flip_horizontal=False, top_labeled=True, show_refs=True, scale=scale, outline_rows=2)}}}}}")
    lines.append("  \\cdot")
    lines.append(f"  \\vcenter{{\\hbox{{{draw_pipe_dream_tikz(rc2, flip_horizontal=False, top_labeled=True, show_refs=True, scale=scale, outline_rows=2)}}}}}")
    lines.append("  \\\\[1em]")
    
    # Right side: sum of RC graphs
    lines.append("  ={}& ")
    
    for i, rc in enumerate(result_rcs):
        # Add the RC graph (4 rows outlined)
        tikz_code = draw_pipe_dream_tikz(rc, flip_horizontal=False, top_labeled=True, show_refs=True, scale=scale, outline_rows=4)
        lines.append(f"  \\vcenter{{\\hbox{{{tikz_code}}}}}")
        
        # Add + sign if not the last term
        if i < len(result_rcs) - 1:
            # Check if we need to break to a new line
            if (i + 1) % terms_per_line == 0:
                lines.append("  \\\\[1em]")
                if i > 0:
                    lines.append("  + & ")
                else:
                    lines.append("  & ")
            else:
                lines.append("  +")
    
    lines.append("\\end{align*}")
    
    return "\n".join(lines)


if __name__ == "__main__":
    # Define the RC graphs
    rc1 = RCGraph([(3, 1), (2,)])
    rc2 = RCGraph([(1,), (2,)])

    result_rcs = [
        RCGraph([(3, 2), (5,), (3,), (4,)]),
        RCGraph([(3, 1), (2,), (3,), (4,)]),
        RCGraph([(3, 2), (4,), (3,), (4,)]),
        RCGraph([(4, 1), (2,), (3,), (4,)]),
        RCGraph([(4, 2), (5,), (3,), (4,)]),
        RCGraph([(5, 1), (2,), (3,), (4,)]),
        RCGraph([(5, 1), (4,), (3,), (4,)]),
        RCGraph([(5, 2), (4,), (3,), (4,)]),
        RCGraph([(5, 4), (5,), (3,), (4,)]),
        RCGraph([(6, 1), (5,), (3,), (4,)]),
        RCGraph([(6, 2), (5,), (3,), (4,)]),
        RCGraph([(6, 4), (5,), (3,), (4,)]),
    ]

    # Generate the equation
    latex_code = generate_multiplication_equation(rc1, rc2, result_rcs, terms_per_line=3, scale=0.4)

    print(latex_code)
    print("\n\n% To use this in a LaTeX document:")
    print("% \\documentclass{article}")
    print("% \\usepackage{amsmath}")
    print("% \\usepackage{tikz}")
    print("% \\begin{document}")
    print("% [paste the generated code here]")
    print("% \\end{document}")
