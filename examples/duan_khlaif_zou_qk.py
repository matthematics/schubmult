#!/usr/bin/env python3
"""Compute the explicit equivariant QK products in arXiv:2606.15261v2.

Run from the repository root with the source tree on the import path::

    PYTHONPATH=src python examples/duan_khlaif_zou_qk.py --list
    PYTHONPATH=src python examples/duan_khlaif_zou_qk.py --example fl4-04
    PYTHONPATH=src python examples/duan_khlaif_zou_qk.py --section fl134
    PYTHONPATH=src python examples/duan_khlaif_zou_qk.py --tex --output products.tex

Inventory (printed page numbers, not PDF page indices):

* Table 3.2, p.14: 15 products in Fl(3).
* Appendix A.1, p.34: 10 selected products in Fl(4).
* Appendix A.2, pp.37--43: 66 distinct products in Fl(1,3;4).
  The product 4123 * 3142 is printed twice on p.42; compute it only once.
* Appendix A.3, pp.47--53: 65 products in Fl(1,2;4).
  The equivariant list omits 3214 * 3214.

These 156 distinct products cover the explicit equivariant Schubert-class QK
products, including ones whose answers contain no quantum parameters. Metrics,
correlators, dual-class decompositions, Whitney relations, and the quantum
cohomology products in Section 4 are different calculations and are not included.
The repeated non-equivariant lists are specializations of these calculations;
--non-equivariant specializes the selected products, not the missing A.3 square.

This generates schubmult's answers, not a transcription or an automated
verification of the paper's right-hand sides. In particular, fl4-04 is the
example discussed with the author. The equivariant quantum K-Pieri
multiplication rule is empirically tested but conjectural.
"""

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from schubmult import ChangeOfVariables, QuantumDoubleGrothendieckRing
from schubmult.abc import x, y, z
from schubmult.mult.groth_quantum_double import apply_kato
from schubmult.symbolic import latex, sstr, sympify_sympy


@dataclass(frozen=True)
class Flag:
    key: str
    name: str
    dimensions: tuple[int, ...]
    n: int
    source: str
    pages: tuple[tuple[int, str], ...]

    @property
    def parabolic(self):
        """Simple reflections inside blocks, rather than the block boundaries."""
        return tuple(i for i in range(1, self.n) if i not in self.dimensions)


# Only the factor permutations and their source pages are recorded. Every
# right-hand side is recomputed symbolically; there are no fitted coefficients.
# Keep the paper's order, even where it is not lexicographic.
FLAGS = (
    Flag("fl3", "Fl(3)", (1, 2), 3, "Table 3.2", (
        (14, "213*213 213*132 213*231 213*312 213*321 "
             "132*132 132*231 132*312 132*321 "
             "231*231 231*312 231*321 312*312 312*321 321*321"),
    )),
    Flag("fl4", "Fl(4)", (1, 2, 3), 4, "Appendix A.1", (
        (34, "1243*1324 1243*3241 1243*2341 2341*3214 1243*1243 "
             "1243*1342 1324*1324 1324*2143 4123*4231 1324*1432"),
    )),
    Flag("fl134", "Fl(1,3;4)", (1, 3), 4, "Appendix A.2", (
        (37, "1243*1243 1243*2134 1243*1342 1243*2143 1243*3124 "
             "1243*2341 1243*3241 1243*4132 1243*4123 1243*3142 "
             "1243*4231 2134*2134 2134*1342 2134*2143"),
        (38, "2134*3124 2134*2341 2134*3241 2134*4132 2134*4123 "
             "2134*3142 2134*4231 1342*1342 1342*2143 1342*3124 "
             "1342*2341 1342*3241 1342*4132 1342*4123 1342*3142 1342*4231"),
        (39, "2143*3124 2143*2143 2143*2341 2143*3241 2143*4132 "
             "2143*4123 2143*3142 2143*4231 3124*3124 3124*2341"),
        (40, "3124*3241 3124*4132 3124*4123 3124*3142 3124*4231 "
             "2341*2341 2341*3241 2341*4132 2341*4123 2341*3142 "
             "2341*4231 3241*3241"),
        (41, "3241*4132 3241*4123 3241*3142 4132*4123 3241*4231"),
        (42, "4132*4132 4132*3142 4123*3142 4132*4231 4123*4123 4123*3142"),
        (43, "4123*4231 3142*3142 3142*4231 4231*4231"),
    )),
    Flag("fl124", "Fl(1,2;4)", (1, 2), 4, "Appendix A.3", (
        (47, "1324*1324 1324*1423 1324*2134 1324*2314 1324*2413 "
             "1324*3124 1324*3214 1324*3412 1324*4123 1324*4213 "
             "1324*4312 1423*1423"),
        (48, "1423*2134 1423*2314 1423*3124 1423*2413 1423*3214 "
             "1423*3412 1423*4123 1423*4213 1423*4312 2134*2134 "
             "2134*2314 2134*2413 2134*3124 2134*3214 2134*3412 "
             "2134*4123 2134*4213 2134*4312"),
        (49, "2314*2314 2314*2413 2314*3124 2314*3214 2314*3412 "
             "2314*4123 2314*4213 2314*4312 2413*2413 2413*3124 "
             "2413*3214 2413*3412"),
        (50, "2413*4123 2413*4213 2413*4312 3124*3124 3124*3214 "
             "3124*3412 3124*4123"),
        (51, "3124*4213 3124*4312 3214*3412 3214*4123 3214*4213 "
             "3214*4312 3412*3412"),
        (52, "3412*4123 3412*4213 3412*4312 4123*4123 4123*4213 4123*4312"),
        (53, "4213*4213 4213*4312 4312*4312"),
    )),
)


@dataclass(frozen=True)
class Example:
    identifier: str
    flag: Flag
    page: int
    left: tuple[int, ...]
    right: tuple[int, ...]
    occurrences: int = 1

    @property
    def reference(self):
        repeated = f"; printed {self.occurrences} times" if self.occurrences > 1 else ""
        return f"{self.flag.source}, p.{self.page}{repeated}"


def examples():
    """Expand the page inventory, retaining a note for repeated products."""
    result = []
    for flag in FLAGS:
        unique = {}
        for page, pairs in flag.pages:
            for pair in pairs.split():
                left, right = (tuple(map(int, word)) for word in pair.split("*"))
                key = tuple(sorted((left, right)))
                if key in unique:
                    old = unique[key]
                    unique[key] = Example(old.identifier, flag, old.page, old.left, old.right, old.occurrences + 1)
                else:
                    unique[key] = Example(f"{flag.key}-{len(unique) + 1:02}", flag, page, left, right)
        result.extend(unique.values())
    return tuple(result)


def make_ring():
    """Use multiplicative torus characters y_i, as in the paper."""
    # The package normally uses additive K-coordinates z_i. At beta=-1,
    # 1-z_i becomes 1/y_i. No substitution of quantum parameters is made here.
    characters = ChangeOfVariables(z, y, lambda i, target: 1 - 1 / target[i])
    return QuantumDoubleGrothendieckRing(x, characters, beta=-1)


def compute(example, ring, non_equivariant=False):
    """Multiply in the quantum ring and take the appropriate finite-flag image."""
    # Do not call .expand(): we want structure constants in the quantum basis,
    # not the explicit x-polynomials. The stable ring also retains w outside S_n.
    product = ring(example.left) * ring(example.right)

    # apply_kato drops w outside S_n, takes minimal right-coset representatives,
    # sets q_j=1 inside parabolic blocks (beta=-1), and renumbers the surviving
    # q's consecutively. Thus for Fl(1,3;4), full-flag q_3 becomes paper q_2;
    # for Fl(1,2;4), q_3 is specialized to 1, not 0.
    # With an empty parabolic set this is just the finite full-flag projection.
    coefficients = apply_kato(product, example.flag.parabolic, n=example.flag.n, beta=-1)
    answer = ring.from_dict(coefficients).simplify()
    if non_equivariant:
        # Multiplicative characters specialize to 1, not 0. Substitute only
        # after coefficient cancellation, leaving all surviving q's symbolic.
        substitution = {y[i]: 1 for i in range(1, example.flag.n + 1)}
        answer = ring.from_dict({w: c.subs(substitution) for w, c in answer.items()}).simplify()
    return answer


def padded_permutation(w, n):
    """Permutation objects trim fixed tails; display the full S_n labels."""
    window = tuple(w)
    # len(Permutation(())) is 2 even though its stored window is empty.
    return window + tuple(range(len(window) + 1, n + 1))


def basis(w, n, tex=False):
    label = ", ".join(map(str, padded_permutation(w, n)))
    if tex:
        return rf"\widetilde{{\mathfrak{{G}}}}_{{({label})}}(x;y)"
    return f"G~_({label})(x;y)"


def terms(answer):
    """Deterministic output, with the longest permutations first."""
    return sorted(answer.items(), key=lambda item: (-item[0].inv, tuple(item[0])))


def render_text(example, answer):
    n = example.flag.n
    lines = [
        f"{example.identifier}: {example.flag.name} ({example.reference})",
        f"{basis(example.left, n)} * {basis(example.right, n)} =",
    ]
    # Explicit coefficient parentheses avoid ambiguities in plain-text email.
    lines.extend(f"  {'  ' if i == 0 else '+ '}({sstr(sympify_sympy(c))}) {basis(w, n)}" for i, (w, c) in enumerate(terms(answer)))
    if not answer:
        lines.append("  0")
    return "\n".join(lines)


def render_tex(example, answer):
    n = example.flag.n
    lines = [
        rf"\subsection*{{{example.identifier}: {example.flag.name}}}",
        example.reference + ".",
        r"\begin{gather*}",
        basis(example.left, n, tex=True) + r"\star_q " + basis(example.right, n, tex=True) + r" = \\",
    ]
    # One basis term per row permits page breaks in long products. The document
    # macro shrinks only rows that would otherwise run outside the text width.
    for i, (w, c) in enumerate(terms(answer)):
        coefficient = latex(sympify_sympy(c))
        sign = "" if i == 0 else "{}+"
        lines.append(r"\fitmath{" + sign + rf"\left({coefficient}\right)" + basis(w, n, tex=True) + r"} \\")
    if not answer:
        lines.append("0")
    lines.append(r"\end{gather*}")
    return "\n".join(lines)


TEX_HEADER = r"""\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,graphicx}
\allowdisplaybreaks
\newsavebox{\equationbox}
\newcommand{\fitmath}[1]{%
  \begingroup
  \sbox{\equationbox}{$\displaystyle #1$}%
  \ifdim\wd\equationbox>\linewidth
    \resizebox{\linewidth}{!}{\usebox{\equationbox}}%
  \else
    \usebox{\equationbox}%
  \fi
  \endgroup}
\begin{document}
\section*{Quantum K-theory examples from arXiv:2606.15261v2}
Computed with \texttt{QuantumDoubleGrothendieckRing}, $\beta=-1$, and
\texttt{ChangeOfVariables} implementing $z_i\mapsto 1-1/y_i$.
The notation $\widetilde{\mathfrak{G}}_w(x;y)$ denotes the package basis,
viewed in the indicated finite flag variety; $\star_q$ denotes its quantum
product. Partial flags use Kato's projection and consecutive quantum indices.
The multiplication rule is empirically tested but conjectural.
These are generated answers, not assertions that the printed formulas agree.
Printed page numbers are used below.

Appendix A.2 repeats $4123\star_q3142$ on p.42 (computed once).
Appendix A.3 omits the equivariant square of $3214$ (not added here).
"""


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--section", choices=[flag.key for flag in FLAGS], help="restrict to one flag variety")
    parser.add_argument("--example", help="restrict to an identifier from --list, e.g. fl4-04")
    parser.add_argument("--list", action="store_true", help="list the selected products without multiplying")
    parser.add_argument("--tex", action="store_true", help="emit a standalone LaTeX document with typeset equations")
    parser.add_argument("--non-equivariant", action="store_true", help="specialize y_i=1 in the selected products")
    parser.add_argument("--output", type=Path, help="write output to this file instead of stdout")
    args = parser.parse_args(argv)
    inventory = examples()
    if args.example and args.example not in {example.identifier for example in inventory}:
        parser.error(f"unknown example: {args.example}; use --list for identifiers")
    selected = [example for example in inventory if (not args.section or example.flag.key == args.section) and (not args.example or example.identifier == args.example)]
    if not selected:
        parser.error("the selected example does not belong to the selected section")
    if args.list and (args.tex or args.non_equivariant):
        parser.error("--list cannot be combined with --tex or --non-equivariant")

    if args.list:
        output = "\n".join(
            f"{e.identifier}: {e.flag.name}, {e.reference}: {''.join(map(str, e.left))} * {''.join(map(str, e.right))}" for e in selected
        ) + "\n"
    else:
        ring = make_ring()
        render = render_tex if args.tex else render_text
        blocks = []
        for i, example in enumerate(selected, start=1):
            print(f"[{i}/{len(selected)}] {example.identifier}: {example.flag.name}", file=sys.stderr, flush=True)
            blocks.append(render(example, compute(example, ring, args.non_equivariant)))
        specialization = "Non-equivariant specialization: y_i=1." if args.non_equivariant else "Equivariant coefficients: z_i -> 1-1/y_i; beta=-1."
        if args.tex:
            specialization = r"Non-equivariant specialization: $y_i=1$." if args.non_equivariant else r"Equivariant coefficients in the multiplicative characters $y_i$."
            output = TEX_HEADER + "\n" + specialization + "\n\n" + "\n\n".join(blocks) + "\n\\end{document}\n"
        else:
            output = (
                "schubmult calculations for arXiv:2606.15261v2\n"
                + specialization
                + "\nQuantumDoubleGrothendieckRing; conjectural quantum K-Pieri rule.\n"
                + "G~_w(x;y) denotes the package basis in the indicated finite flag variety.\n"
                + "Generated answers; no automatic comparison with the printed right-hand sides.\n\n"
                + "\n\n".join(blocks) + "\n"
            )
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    else:
        sys.stdout.write(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
