"""
Flatten all Python modules under src/schubmult into the package root.

Behavior:
- For each Python file src/schubmult/.../foo.py (not already at src/schubmult/foo.py),
  create a flattened filename by joining path components with '__', e.g.
  sub/inner/foo.py -> src/schubmult/sub__inner__foo.py and module name
  schubmult.sub__inner__foo.
- Build a mapping original_module_name -> new_module_name and update all
  import statements across the repository (Import and ImportFrom nodes)
  to use the new module names.
- Move files into src/schubmult with the new names.
- Supports --dry-run to preview planned changes without writing.

Caveats:
- Relative imports (from .sub import ...) are left unchanged.
- The script rewrites only syntactic import nodes (AST). String-based dynamic
  imports won't be updated.
- Back up your repo or run on a branch. Review diffs before committing.
"""
import argparse
import ast
import os
import shutil
from pathlib import Path
from typing import Dict, Tuple

REPO_ROOT = Path(__file__).resolve().parent
PKG_ROOT = REPO_ROOT / "schubmult"
if not PKG_ROOT.exists():
    raise SystemExit(f"Package root not found: {PKG_ROOT}")

def collect_nested_modules(pkg_root: Path) -> Dict[str, Path]:
    mapping = {}
    for p in pkg_root.rglob("*.py"):
        # skip files already at package root
        if p.parent == pkg_root:
            continue
        # skip package __init__ files
        if p.name == "__init__.py":
            continue
        # skip anything under the scripts directory entirely
        # (leave src/schubmult/scripts untouched)
        rel_parts = p.relative_to(pkg_root).parts
        if "scripts" in rel_parts:
            continue

        rel = p.relative_to(pkg_root)
        mod = "schubmult." + ".".join(rel.with_suffix("").parts)
        flattened_name = "__".join(rel.with_suffix("").parts)
        new_mod = "schubmult." + flattened_name
        mapping[mod] = (p, flattened_name + ".py")
    return mapping

class ImportRewriter(ast.NodeTransformer):
    def __init__(self, modmap: Dict[str,str]):
        super().__init__()
        # sort keys by length desc so longest matches are replaced first
        self.modmap = {k: v for k,v in modmap.items()}
        self.sorted_keys = sorted(self.modmap.keys(), key=len, reverse=True)

    def _remap(self, name: str) -> str:
        if not name:
            return name
        for old in self.sorted_keys:
            if name == old or name.startswith(old + "."):
                newbase = self.modmap[old]
                suffix = name[len(old):]  # includes leading dot if any
                return newbase + suffix
        return name

    def visit_Import(self, node: ast.Import):
        changed = False
        for alias in node.names:
            newname = self._remap(alias.name)
            if newname != alias.name:
                alias.name = newname
                changed = True
        return node

    def visit_ImportFrom(self, node: ast.ImportFrom):
        # leave relative imports (level > 0) alone
        if node.level and node.level > 0:
            return node
        node_module = node.module or ""
        newmod = self._remap(node_module)
        if newmod != node_module:
            node.module = newmod
        return node

def rewrite_imports_in_file(path: Path, modmap: Dict[str,str], dry_run: bool) -> Tuple[bool,str]:
    src = path.read_text(encoding="utf8")
    try:
        tree = ast.parse(src)
    except SyntaxError as e:
        return False, f"SKIP_SYNTAX_ERR: {path}: {e}"
    rewriter = ImportRewriter(modmap)
    new_tree = rewriter.visit(tree)
    ast.fix_missing_locations(new_tree)
    try:
        new_src = ast.unparse(new_tree)
    except Exception:
        # fallback: if ast.unparse not available or fails, skip
        return False, f"SKIP_UNPARSE: {path}"
    if new_src != src:
        if dry_run:
            return True, f"WOULD_REWRITE: {path}"
        path.write_text(new_src, encoding="utf8")
        return True, f"REWRITTEN: {path}"
    return False, f"UNCHANGED: {path}"

def main(dry_run: bool):
    mappinginfo = collect_nested_modules(PKG_ROOT)
    if not mappinginfo:
        print("No nested modules found; nothing to do.")
        return

    # build modmap original -> new_module_name
    modmap = {}
    moves = []  # tuples (src_path, dest_path)
    for orig_mod, (src_path, new_filename) in mappinginfo.items():
        new_path = PKG_ROOT / new_filename
        # if dest exists, add numeric suffix
        idx = 1
        base = new_path.with_suffix("").name
        while new_path.exists():
            new_path = PKG_ROOT / f"{base}_{idx}.py"
            idx += 1
        new_mod = "schubmult." + new_path.with_suffix("").name
        modmap[orig_mod] = new_mod
        moves.append((src_path, new_path))

    print("Planned module remapping (original -> flattened):")
    for k, v in modmap.items():
        print(f"  {k}  ->  {v}")

    # Rewrite imports across all .py files under src, but exclude scripts directory
    all_pyfiles = list((REPO_ROOT / "src").rglob("*.py"))
    pyfiles = []
    scripts_prefix = (PKG_ROOT / "scripts").resolve()
    for f in all_pyfiles:
        try:
            # if f is under src/schubmult/scripts, skip it
            if scripts_prefix in f.resolve().parents or f.resolve().parent == scripts_prefix:
                continue
        except Exception:
            # on resolution errors, conservatively include file
            pass
        pyfiles.append(f)

    print(f"\nScanning {len(pyfiles)} python files to rewrite imports (scripts/ excluded)...")
    any_changes = []
    for f in pyfiles:
        changed, msg = rewrite_imports_in_file(f, modmap, dry_run=dry_run)
        if changed:
            any_changes.append(msg)
            print(msg)
    if not any_changes:
        print("No import rewrites necessary/needed.")

    print("\nPlanned file moves:")
    for s, d in moves:
        print(f"  {s} -> {d}")

    if dry_run:
        print("\nDry-run complete. No files were moved or rewritten.")
        return

    # perform moves
    for s, d in moves:
        d_parent = d.parent
        d_parent.mkdir(parents=True, exist_ok=True)
        print(f"Moving {s} -> {d}")
        shutil.move(str(s), str(d))

    print("\nFlattening complete. You should inspect changes, run tests and update any packaging metadata if needed.")
    print("Suggested next steps:")
    print("  git status")
    print("  git add -A")
    print('  git commit -m "Flatten schubmult nested modules into package root (flatten tool)"')
    print("  pytest -q")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true", help="Preview changes without writing")
    args = p.parse_args()
    main(dry_run=args.dry_run)