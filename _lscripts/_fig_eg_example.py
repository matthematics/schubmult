"""TikZ for the Edelman-Greene example: R = {(1,2),(1,3),(2,2)} and f_1(R), with P and the negated recording tableau."""
import sys

sys.path.insert(0, "_lscripts")
from _eg_convention_check import eg  # noqa: E402

from schubmult import RCGraph  # noqa: E402
from schubmult.visualization import draw_pipe_dream_tikz  # noqa: E402


def tab(T, x0=0, y0=0, s=0.55):
    out = []
    for ri, row in enumerate(T):
        for ci, v in enumerate(row):
            x, y = x0 + ci * s, y0 - ri * s
            out.append(f"  \\draw ({x:.2f},{y:.2f}) rectangle ({x + s:.2f},{y - s:.2f});")
            out.append(f"  \\node[font=\\small] at ({x + s / 2:.2f},{y - s / 2:.2f}) {{${v}$}};")
    return "\n".join(out)


def data(R):
    rows = [tuple(r) for r in R]
    word = [l for r in rows for l in r]
    ridx = [k + 1 for k, r in enumerate(rows) for _ in r]
    return eg(word, [-x for x in ridx])


R = RCGraph([(3, 2), (3,), ()])
F = R.lowering_operator(1)
for name, X in (("R", R), ("f1R", F)):
    P, T = data(X)
    print("%", name, [tuple(r) for r in X], "word", [l for r in X for l in r], "P", P, "T", T)
    print(draw_pipe_dream_tikz(X, max_size=4, scale=0.6, outline_rows=3, show_refs=True, flip_horizontal=False))
    print("% P")
    print("\\begin{tikzpicture}[scale=1]\n" + tab(P) + "\n\\end{tikzpicture}")
    print("% T")
    print("\\begin{tikzpicture}[scale=1]\n" + tab(T) + "\n\\end{tikzpicture}")
