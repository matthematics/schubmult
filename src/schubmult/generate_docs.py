"""
Minimal, robust docs generator for this repo.
Writes docs/API.md and per-module files to docs/modules/.
"""
from pathlib import Path
import ast
import sys

THIS = Path(__file__).resolve()
PROJECT_ROOT = THIS.parents[2]  # repo_root
SRC = PROJECT_ROOT / "src"
OUT_DIR = PROJECT_ROOT / "docs"
MODULES_DIR = OUT_DIR / "modules"
OUT_DIR.mkdir(exist_ok=True)
MODULES_DIR.mkdir(exist_ok=True)

def process_file(path):
    try:
        src = path.read_text(encoding="utf8")
    except Exception as e:
        return {"error": f"read_error: {e}"}
    try:
        tree = ast.parse(src)
    except Exception as e:
        return {"error": f"parse_error: {e}"}
    mod_doc = ast.get_docstring(tree) or ""
    items = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            doc = ast.get_docstring(node) or ""
            items.append({
                "type": node.__class__.__name__,
                "name": getattr(node, "name", "<anon>"),
                "lineno": getattr(node, "lineno", None),
                "doc_first": (doc.splitlines()[0] if doc else ""),
            })
    return {"module_doc": mod_doc, "items": items}

def write_module_md(relpath, info):
    outpath = MODULES_DIR / (relpath.with_suffix(".md").name)
    lines = [f"<!-- filepath: {relpath} -->", ""]
    if "error" in info:
        lines.append("ERROR: " + info["error"])
    else:
        if info["module_doc"]:
            lines.append(info["module_doc"].strip())
            lines.append("")
        for it in info.get("items", []):
            lines.append(f"- `{it['type']}` — `{it['name']}` (line {it['lineno']})")
            if it["doc_first"]:
                lines.append(f"    - {it['doc_first']}")
    outpath.write_text("\n".join(lines), encoding="utf8")

def main():
    api_lines = ["# schubmult API index", ""]
    py_files = sorted(SRC.rglob("*.py"))
    for p in py_files:
        if any(part in ("venv", ".venv", "__pycache__", "build", "docs") for part in p.parts):
            continue
        rel = p.relative_to(PROJECT_ROOT)
        info = process_file(p)
        api_lines.append(f"## {rel}")
        api_lines.append("")
        if "error" in info:
            api_lines.append(f"- parse error: {info['error']}")
            api_lines.append("")
        else:
            for it in info.get("items", []):
                api_lines.append(f"- `{it['type']}` — `{it['name']}` (line {it['lineno']})")
            api_lines.append("")
        write_module_md(rel, info)
    (OUT_DIR / "API.md").write_text("\n".join(api_lines), encoding="utf8")
    print("Docs written to", OUT_DIR)

if __name__ == "__main__":
    main()