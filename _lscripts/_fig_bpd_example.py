"""Generate TikZ for Phi(R) and Phi(Z(R)) for the paper's Figure example R = {(1,2),(1,3),(2,2)}."""
from schubmult import RCGraph
from schubmult.combinatorics.bpd import BPD, TileType


def tikz(B):
    n = B.rows
    out = ["\\begin{tikzpicture}[scale=0.7]"]
    for i in range(n + 1):
        out.append(f"  \\draw[lightgray, very thin] (0,{i}) -- ({n},{i});")
        out.append(f"  \\draw[lightgray, very thin] ({i},0) -- ({i},{n});")
    for i in range(n):
        for j in range(n):
            t = B._grid[i, j]
            x, y = j, n - 1 - i
            cx, cy = x + 0.5, y + 0.5
            if t in (TileType.HORIZ, TileType.CROSS):
                out.append(f"  \\draw[blue, line width=1.0pt] ({x},{cy}) -- ({x + 1},{cy});")
            if t in (TileType.VERT, TileType.CROSS):
                out.append(f"  \\draw[blue, line width=1.0pt] ({cx},{y}) -- ({cx},{y + 1});")
            if t == TileType.ELBOW_SE:
                out.append(f"  \\draw[blue, line width=1.0pt] ({cx},{y}) .. controls ({cx},{y + 0.3}) and ({x + 0.7},{cy}) .. ({x + 1},{cy});")
            if t == TileType.ELBOW_NW:
                out.append(f"  \\draw[blue, line width=1.0pt] ({x},{cy}) .. controls ({x + 0.3},{cy}) and ({cx},{y + 0.7}) .. ({cx},{y + 1});")
            if t == TileType.BLANK:
                out.append(f"  \\fill[gray!20] ({x + 0.1},{y + 0.1}) rectangle ({x + 0.9},{y + 0.9});")
    for j in range(n):
        out.append(f"  \\node[font=\\scriptsize] at ({j + 0.5},-0.35) {{{j + 1}}};")
    for i in range(n):
        out.append(f"  \\node[font=\\scriptsize, anchor=west] at ({n + 0.15},{n - 1 - i + 0.5}) {{{B.perm[i]}}};")
    out.append(f"  \\draw[black, line width=1pt] (0,0) rectangle ({n},{n});")
    out.append("\\end{tikzpicture}")
    return "\n".join(out)


R = RCGraph([(3, 2), (3,), (), ()])
Z = RCGraph([(3, 1), (2,), (), ()])
B = BPD.from_rc_graph(R).resize(4)
C = BPD.from_rc_graph(Z).resize(4)
with open("/tmp/fig_bpd.tex", "w") as f:
    f.write(tikz(B) + "\n%%%%\n" + tikz(C) + "\n")
print(B.perm, C.perm)
