from __future__ import annotations

import argparse
import shlex
import sys

from schubmult.combinatorics.anti_rc_graph import AntiRCGraph
from schubmult.combinatorics.rc_graph import RCGraph


def _parse_int_csv(text: str) -> tuple[int, ...]:
    text = text.strip()
    if text == "":
        return ()
    return tuple(int(piece.strip()) for piece in text.split(",") if piece.strip())


def _parse_rows_spec(spec: str) -> tuple[tuple[int, ...], ...]:
    rows = []
    for row_text in spec.split(";"):
        row_text = row_text.strip()
        if row_text == "":
            rows.append(())
        else:
            rows.append(_parse_int_csv(row_text))
    return tuple(rows)


def _grid_values(obj) -> list[list[int | None]]:
    return [[obj[i, j] for j in range(obj.cols)] for i in range(obj.rows)]


def _print_summary(anti: AntiRCGraph, *, title: str) -> None:
    print("=" * 78)
    print(title)
    print("-" * 78)
    print("stored rows:", tuple(anti))
    print("anti display:")
    print(anti)
    print("anti display matrix values:")
    print(_grid_values(anti))
    print("reflection view (ordinary RCGraph display):")
    print(anti.reflection_view)
    print("anti reduced word:", anti.anti_reduced_word)
    print("anti compatible sequence:", anti.anti_compatible_sequence)
    print("anti permutation:", anti.anti_permutation)
    print("anti valid:", anti.anti_is_valid)
    print("=" * 78)


def _demo_quick() -> None:
    anti = AntiRCGraph([(4, 2), (3,), (1,)])
    _print_summary(anti, title="Demo 1: Build from stored rows")

    word = (3, 2, 2, 1)
    seq = (3, 3, 2, 1)
    anti2 = AntiRCGraph.from_reduced_anticompatible(word, seq)
    _print_summary(anti2, title="Demo 2: Build from anti reduced+compatible data")


def _demo_all() -> None:
    _demo_quick()

    rc = RCGraph([(3, 1), (2,), ()])
    anti = AntiRCGraph.from_rc_graph(rc)
    _print_summary(anti, title="Demo 3: Convert from ordinary RCGraph")
    print("converted back to RCGraph rows:", tuple(anti.to_rc_graph()))


def _run_one(rows: str | None, word: str | None, seq: str | None, length: int | None) -> None:
    if rows is not None:
        anti = AntiRCGraph(_parse_rows_spec(rows))
        _print_summary(anti, title="Custom run from --rows")
        return

    if word is not None or seq is not None:
        if word is None or seq is None:
            raise ValueError("Provide both --word and --seq together.")
        anti = AntiRCGraph.from_reduced_anticompatible(_parse_int_csv(word), _parse_int_csv(seq), length=length)
        _print_summary(anti, title="Custom run from --word/--seq")
        return

    _demo_quick()


def _run_repl() -> None:
    print("AntiRCGraph showcase REPL")
    print("Commands:")
    print("  rows <spec>                     e.g. rows 4,2;3;1")
    print("  ws <word_csv> <seq_csv> [len]  e.g. ws 3,2,2,1 3,3,2,1 4")
    print("  quick                           run quick demo")
    print("  all                             run all demos")
    print("  quit                            exit")

    while True:
        try:
            line = input("anti> ").strip()
        except EOFError:
            print()
            return
        if not line:
            continue

        parts = shlex.split(line)
        cmd = parts[0].lower()

        try:
            if cmd in {"quit", "q", "exit"}:
                return
            if cmd == "quick":
                _demo_quick()
                continue
            if cmd == "all":
                _demo_all()
                continue
            if cmd == "rows":
                if len(parts) < 2:
                    print("usage: rows <spec>")
                    continue
                anti = AntiRCGraph(_parse_rows_spec(parts[1]))
                _print_summary(anti, title="REPL run from rows")
                continue
            if cmd == "ws":
                if len(parts) < 3:
                    print("usage: ws <word_csv> <seq_csv> [length]")
                    continue
                length = int(parts[3]) if len(parts) >= 4 else None
                anti = AntiRCGraph.from_reduced_anticompatible(_parse_int_csv(parts[1]), _parse_int_csv(parts[2]), length=length)
                _print_summary(anti, title="REPL run from word/seq")
                continue

            print(f"Unknown command: {cmd}")
        except Exception as exc:
            print(f"error: {exc}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Showcase and playground for AntiRCGraph drawing and encoding.")
    parser.add_argument("--rows", help='Stored row spec, e.g. "4,2;3;1" (semicolon separates rows).')
    parser.add_argument("--word", help='Anti reduced word CSV, e.g. "3,2,2,1".')
    parser.add_argument("--seq", help='Anti compatible sequence CSV, e.g. "3,3,2,1".')
    parser.add_argument("--length", type=int, default=None, help="Optional target number of rows when using --word/--seq.")
    parser.add_argument("--demo", choices=["quick", "all"], default="quick", help="Which built-in demo set to run when no custom input is provided.")
    parser.add_argument("--repl", action="store_true", help="Launch interactive mode after initial run.")
    return parser


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv
    args = build_parser().parse_args(argv[1:])

    if args.rows is not None or args.word is not None or args.seq is not None:
        _run_one(args.rows, args.word, args.seq, args.length)
    else:
        if args.demo == "quick":
            _demo_quick()
        else:
            _demo_all()

    if args.repl:
        _run_repl()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
