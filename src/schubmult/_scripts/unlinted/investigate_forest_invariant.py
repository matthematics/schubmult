from schubmult import *


def _html_escape(text):
        import html

        return html.escape(str(text))


def _render_html(report_rows, n):
        import datetime

        generated_at = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        total_perms = len(report_rows)
        total_invariants = sum(row["num_invariants"] for row in report_rows)
        total_rc_graphs = sum(row["num_rc_graphs"] for row in report_rows)

        section_html = []
        for row in report_rows:
                invariant_rows = []
                for invariant, rcs in row["invariants"]:
                        rc_list = "\n".join(_html_escape(rc) for rc in rcs)
                        invariant_rows.append(
                                """
                                <tr>
                                    <td><code>{invariant}</code></td>
                                    <td>{count}</td>
                                    <td>
                                        <details>
                                            <summary>Show RC graphs</summary>
                                            <pre>{rc_list}</pre>
                                        </details>
                                    </td>
                                </tr>
                                """.format(invariant=_html_escape(invariant), count=len(rcs), rc_list=rc_list)
                        )

                section_html.append(
                        """
                        <section class="perm-block">
                            <h2>Permutation {perm}</h2>
                            <div class="meta">
                                <span><strong>trimcode:</strong> <code>{trimcode}</code></span>
                                <span><strong>RC graphs:</strong> {num_rc_graphs}</span>
                                <span><strong>Invariants:</strong> {num_invariants}</span>
                            </div>
                            <table>
                                <thead>
                                    <tr>
                                        <th>Invariant</th>
                                        <th>Count</th>
                                        <th>RC Graphs</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {invariant_rows}
                                </tbody>
                            </table>
                        </section>
                        """.format(
                                perm=_html_escape(row["perm"]),
                                trimcode=_html_escape(row["trimcode"]),
                                num_rc_graphs=row["num_rc_graphs"],
                                num_invariants=row["num_invariants"],
                                invariant_rows="\n".join(invariant_rows),
                        )
                )

        return """
<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>Forest Invariant Report (n={n})</title>
    <style>
        :root {{
            --bg: #0f172a;
            --surface: #111827;
            --muted: #94a3b8;
            --text: #e2e8f0;
            --accent: #38bdf8;
            --border: #1f2937;
        }}
        body {{
            margin: 0;
            background: linear-gradient(180deg, #020617 0%, #0f172a 100%);
            color: var(--text);
            font-family: Inter, ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif;
            line-height: 1.45;
        }}
        .container {{ max-width: 1200px; margin: 0 auto; padding: 2rem 1rem 3rem; }}
        h1 {{ margin: 0 0 0.4rem; font-size: 1.8rem; }}
        .subtitle {{ color: var(--muted); margin-bottom: 1.2rem; }}
        .summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 0.8rem;
            margin-bottom: 1.2rem;
        }}
        .card {{
            background: rgba(17, 24, 39, 0.9);
            border: 1px solid var(--border);
            border-radius: 12px;
            padding: 0.8rem 0.9rem;
        }}
        .card .label {{ color: var(--muted); font-size: 0.85rem; }}
        .card .value {{ font-size: 1.2rem; font-weight: 700; margin-top: 0.2rem; }}
        .perm-block {{
            background: rgba(17, 24, 39, 0.95);
            border: 1px solid var(--border);
            border-radius: 14px;
            padding: 1rem;
            margin-bottom: 1rem;
        }}
        .perm-block h2 {{ margin: 0 0 0.6rem; font-size: 1.15rem; color: var(--accent); }}
        .meta {{ display: flex; flex-wrap: wrap; gap: 0.9rem; color: var(--muted); margin-bottom: 0.8rem; }}
        table {{ width: 100%; border-collapse: collapse; overflow: hidden; border-radius: 10px; }}
        th, td {{ border: 1px solid var(--border); padding: 0.55rem; vertical-align: top; }}
        th {{ background: rgba(30, 41, 59, 0.7); text-align: left; }}
        code {{ background: rgba(15, 23, 42, 0.9); border: 1px solid var(--border); padding: 0.1rem 0.35rem; border-radius: 6px; }}
        details summary {{ cursor: pointer; color: var(--accent); }}
        pre {{
            margin: 0.5rem 0 0;
            max-height: 280px;
            overflow: auto;
            white-space: pre-wrap;
            background: rgba(2, 6, 23, 0.9);
            border: 1px solid var(--border);
            padding: 0.6rem;
            border-radius: 8px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Forest Invariant Report</h1>
        <div class="subtitle">n={n} · generated at {generated_at}</div>

        <div class="summary">
            <div class="card"><div class="label">Permutations</div><div class="value">{total_perms}</div></div>
            <div class="card"><div class="label">RC Graphs (total)</div><div class="value">{total_rc_graphs}</div></div>
            <div class="card"><div class="label">Invariants (total)</div><div class="value">{total_invariants}</div></div>
        </div>

        {sections}
    </div>
</body>
</html>
""".format(
                n=n,
                generated_at=generated_at,
                total_perms=total_perms,
                total_rc_graphs=total_rc_graphs,
                total_invariants=total_invariants,
                sections="\n".join(section_html),
        )


if __name__ == "__main__":
        import argparse
        from pathlib import Path

        parser = argparse.ArgumentParser(description="Generate forest-invariant report as HTML.")
        parser.add_argument("n", type=int, help="Permutation size")
        parser.add_argument("--out", type=str, default=None, help="Output HTML file path")
        args = parser.parse_args()

        n = args.n
        out_path = Path(args.out) if args.out else Path(f"forest_invariant_report_n{n}.html")

        perms = Permutation.all_permutations(n)
        report_rows = []
        for perm in perms:
                invariant_to_rcs = {}
                rc_graphs = list(RCGraph.all_rc_graphs(perm))
                for rc in rc_graphs:
                        invariant = rc.forest_invariant
                        invariant_to_rcs[invariant] = invariant_to_rcs.get(invariant, set())
                        invariant_to_rcs[invariant].add(rc)

                ordered_invariants = sorted(
                        ((inv, sorted(list(rcs), key=lambda x: str(x))) for inv, rcs in invariant_to_rcs.items()),
                        key=lambda x: str(x[0]),
                )

                report_rows.append(
                        {
                                "perm": perm,
                                "trimcode": perm.trimcode,
                                "num_rc_graphs": len(rc_graphs),
                                "num_invariants": len(invariant_to_rcs),
                                "invariants": ordered_invariants,
                        }
                )

        html_report = _render_html(report_rows, n)
        out_path.write_text(html_report, encoding="utf-8")
        print(f"Wrote report to {out_path.resolve()}")