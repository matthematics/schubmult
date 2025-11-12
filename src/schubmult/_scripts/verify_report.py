import argparse
import base64
import html
import json
import os
import pickle
import sys
from datetime import datetime
from typing import Any, Dict, Tuple


def json_key_load(s: str) -> Any:
    """Inverse of the json_key_rep used by the verification saver."""
    if not isinstance(s, str):
        raise TypeError("json_key_load expects a string")
    if s.startswith("PICKLE:"):
        payload = s[len("PICKLE:") :]
        raw = base64.b64decode(payload.encode("ascii"))
        return pickle.loads(raw)
    if s.startswith("REPR:"):
        rep = s[len("REPR:") :]
        try:
            return eval(rep)
        except Exception:
            return rep
    try:
        return eval(s)
    except Exception:
        return s


def load_verification_json(path: str) -> Tuple[Dict[Any, Any], Dict[str, Any]]:
    """
    Load verification file saved by lr_rule_verify:
      - new format: {"meta": {...}, "records": {json_key: bool}}
      - legacy: flat mapping {json_key: bool}
    Returns (records: dict[decoded_key -> value], meta: dict)
    """
    with open(path, "r") as f:
        obj = json.load(f)
    if isinstance(obj, dict) and "records" in obj:
        raw = obj.get("records", {})
        meta = obj.get("meta", {}) or {}
    else:
        raw = obj
        meta = {}
    records = {}
    for k, v in raw.items():
        try:
            dk = json_key_load(k)
        except Exception:
            dk = f"<unreadable key: {k!r}>"
        records[dk] = v
    return records, meta


def key_to_columns(key: Any) -> Tuple[str, ...]:
    """
    Try to produce useful columns for typical key shapes used by the verifier.
    E.g. keys often are (u, v, w, n) where u,v,w are permutations/objects and n is int.
    If unpackable to 4 items, return their reprs as columns, else return single repr.
    """
    if isinstance(key, tuple) and len(key) == 4:
        u, v, w, n = key
        return (repr(u), repr(v), repr(w), str(n))
    # common alt: (hw_tab, perm, n) -> 3-tuple
    if isinstance(key, tuple) and len(key) == 3:
        a, b, c = key
        return (repr(a), repr(b), str(c))
    # fallback: single column
    return (repr(key),)


def make_html_report(records: Dict[Any, Any], meta: Dict[str, Any], title: str = "LR verification report") -> str:
    total = len(records)
    num_true = sum(1 for v in records.values() if v is True)
    num_false = sum(1 for v in records.values() if v is False)
    num_other = total - num_true - num_false
    timestamp = datetime.utcnow().isoformat() + "Z"

    css = """
    body{font-family: Arial, Helvetica, sans-serif; margin:20px}
    table{border-collapse:collapse;width:100%}
    th,td{border:1px solid #ddd;padding:8px}
    th{background:#f2f2f2;text-align:left}
    tr.fail td{background:#ffecec}
    tr.ok td{background:#eaffea}
    pre{white-space:pre-wrap;word-break:break-word}
    .meta{margin-bottom:16px;padding:8px;border:1px solid #ddd;background:#fafafa}
    """

    html_parts = []
    html_parts.append(f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title>")
    html_parts.append(f"<style>{css}</style></head><body>")
    html_parts.append(f"<h1>{html.escape(title)}</h1>")
    html_parts.append(f"<p>Generated: {html.escape(timestamp)}</p>")

    # Meta
    html_parts.append("<div class='meta'><h2>Test conditions (meta)</h2>")
    if meta:
        html_parts.append("<dl>")
        for k, v in sorted(meta.items()):
            html_parts.append(f"<dt><strong>{html.escape(str(k))}</strong></dt><dd><pre>{html.escape(repr(v))}</pre></dd>")
        html_parts.append("</dl>")
    else:
        html_parts.append("<p><em>No meta found in file.</em></p>")
    html_parts.append("</div>")

    # Summary
    html_parts.append("<h2>Summary</h2>")
    html_parts.append("<ul>")
    html_parts.append(f"<li>Total recorded entries: {total}</li>")
    html_parts.append(f"<li>Verified OK (True): {num_true}</li>")
    html_parts.append(f"<li>Failures (False): {num_false}</li>")
    html_parts.append(f"<li>Other / undecoded: {num_other}</li>")
    html_parts.append("</ul>")

    # Failures table
    html_parts.append("<h2>Failures</h2>")
    if num_false == 0:
        html_parts.append("<p><strong>No failures found.</strong></p>")
    else:
        html_parts.append("<table>")
        # table headers: try to detect common 4-col keys
        sample_keys = [k for k, v in records.items()][:20]
        sample_cols = None
        for sk in sample_keys:
            if isinstance(sk, tuple) and len(sk) == 4:
                sample_cols = ("u", "v", "w", "n")
                break
            if isinstance(sk, tuple) and len(sk) == 3:
                sample_cols = ("hw_tab", "perm", "n")
                break
        if sample_cols:
            html_parts.append("<tr>" + "".join(f"<th>{html.escape(c)}</th>" for c in sample_cols) + "<th>value</th></tr>")
            for k, v in records.items():
                if v is False:
                    cols = key_to_columns(k)
                    cls = "fail"
                    html_parts.append(f"<tr class='{cls}'>" + "".join(f"<td><pre>{html.escape(c)}</pre></td>" for c in cols) + f"<td><pre>{html.escape(repr(v))}</pre></td></tr>")
        else:
            html_parts.append("<tr><th>key</th><th>value</th></tr>")
            for k, v in records.items():
                if v is False:
                    html_parts.append("<tr class='fail'><td><pre>" + html.escape(repr(k)) + "</pre></td><td><pre>" + html.escape(repr(v)) + "</pre></td></tr>")
        html_parts.append("</table>")

    # Optionally show successes (collapsible)
    html_parts.append("<h2>Verified OK (sample)</h2>")
    html_parts.append("<details><summary>Show up to 200 OK entries</summary><pre>")
    ok_count = 0
    for k, v in records.items():
        if v is True:
            html_parts.append(html.escape(repr(k)) + "\n")
            ok_count += 1
            if ok_count >= 200:
                html_parts.append("...truncated...\n")
                break
    html_parts.append("</pre></details>")

    # Provide raw JSON download link (embed as data URL if small)
    try:
        raw_text = json.dumps({"meta": meta, "records": {str(k): v for k, v in records.items()}}, indent=2)
        if len(raw_text) < 200000:
            data_url = "data:application/json;charset=utf-8," + html.escape(raw_text)
            html_parts.append(f"<p><a href='{data_url}' download='verification_export.json'>Download raw (embedded)</a></p>")
    except Exception:
        pass

    html_parts.append("</body></html>")
    return "\n".join(html_parts)


def main(argv=None):
    p = argparse.ArgumentParser(description="Render verification JSON to an HTML report")
    p.add_argument("json_path", help="Path to verification JSON (filename.json)")
    p.add_argument("-o", "--output", help="HTML output path (defaults to stdout)")
    args = p.parse_args(argv)

    path = args.json_path
    if not os.path.exists(path):
        print(f"file not found: {path}", file=sys.stderr)
        return 2

    records, meta = load_verification_json(path)
    report = make_html_report(records, meta, title=f"LR verification report: {os.path.basename(path)}")

    if args.output:
        with open(args.output, "w", encoding="utf8") as f:
            f.write(report)
        print(f"Wrote HTML report to {args.output}")
    else:
        sys.stdout.write(report)


if __name__ == "__main__":
    sys.exit(main())