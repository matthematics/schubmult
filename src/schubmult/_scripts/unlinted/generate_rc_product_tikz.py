"""Generate TikZ code for a specific RC-graph product equation.

Requested equation:
  RCGraph([(2,),(),(e,)]) \star RCGraph([(3,2,1),(),(3,)])
= RCGraph([(4,3,2,1),(),(4,3)])
+ RCGraph([(4,3,2,1),(4,),(3,)])
+ RCGraph([(4,3,2,1),(4,3),()])
- RCGraph([(4,3,2,1),(2,),(3,)])

Notes:
- We render all diagrams at fixed width = 5.
- Cells above the antidiagonal are explicitly filled with bumps (elbows), so
  every square has either a crossing or a bump.
"""

from __future__ import annotations

from schubmult import RCGraph
from schubmult.visualization import draw_pipe_dream_tikz


def _draw_single_rc_tikz(
    rc: RCGraph,
    *,
    width: int = 5,
    scale: float = 0.58,
    flip_horizontal: bool = True,
    show_refs: bool = False,
) -> str:
    return draw_pipe_dream_tikz(
            rc=rc,
            max_size=width,
            flip_horizontal=False,
            top_labeled=True,
            show_refs=show_refs,
            scale=scale
        )


def generate_requested_equation(width: int = 5, scale: float = 0.58) -> str:
    # Interpreting the provided (e,) row as (3,) for the 3rd row crossing.
    # rc1 = RCGraph([(2,), (), (3,)])
    # rc2 = RCGraph([(3, 2, 1), (), (3,)])
    rc1, rc2 = (RCGraph([
( ),
(2,),
(3,)]), RCGraph([
(3, 2,1,),
(     ),
(    3,)]))

    positives = [
        # RCGraph([(4, 3, 2, 1), (), (4, 3)]),
        # RCGraph([(4, 3, 2, 1), (4,), (3,)]),
        # RCGraph([(4, 3, 2, 1), (4, 3), ()]),
        RCGraph([
(3,2,1,),
(    2,),
(  4,3,)]) , RCGraph([
(3,2,1,),
(4,   2,),
(    3,)])
    ]
    negative = RCGraph([(4, 3, 2, 1), (2,), (3,)])

    lhs1 = _draw_single_rc_tikz(rc1, width=width, scale=scale)
    lhs2 = _draw_single_rc_tikz(rc2, width=width, scale=scale)
    rhs = [_draw_single_rc_tikz(rc, width=width, scale=scale) for rc in positives]
    rhs_neg = _draw_single_rc_tikz(negative, width=width, scale=scale)

    lines: list[str] = []
    lines.append("\\[")
    lines.append(f"\\vcenter{{\\hbox{{{lhs1}}}}}")
    lines.append("\\;\\star\\;")
    lines.append(f"\\vcenter{{\\hbox{{{lhs2}}}}}")
    lines.append("\\;=\\;")
    lines.append(f"\\vcenter{{\\hbox{{{rhs[0]}}}}}")
    lines.append("\\;+\\;")
    lines.append(f"\\vcenter{{\\hbox{{{rhs[1]}}}}}")
    # lines.append("\\;+\\;")
    # lines.append(f"\\vcenter{{\\hbox{{{rhs[2]}}}}}")
    # lines.append("\\;-\\;")
    # lines.append(f"\\vcenter{{\\hbox{{{rhs_neg}}}}}")
    lines.append("\\]")

    return "\n".join(lines)


if __name__ == "__main__":
    print("% Paste this into your TeX file (requires \\usepackage{tikz})")
    print(generate_requested_equation(width=5, scale=0.52))
