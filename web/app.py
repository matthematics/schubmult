"""Embeddable web wrapper around schubmult_py and schubmult_double.

Run with:

    conda activate schubmult_312
    python web/app.py

Then open http://localhost:5000/

The page at /embed renders only the form + results (no <html> chrome) so it
can be dropped into an <iframe> from another site.

Configuration via environment variables:

    SCHUBMULT_ALLOWED_ORIGINS   Comma-separated list of origins allowed to
                                embed via <iframe> (frame-ancestors).
                                Default: 'self' only.
                                Example: "https://your-site.com,https://www.your-site.com"
    SCHUBMULT_COMPUTE_TIMEOUT   Per-request compute timeout in seconds.
                                Default: 8.
    SCHUBMULT_MAX_PERM_LENGTH   Max integers across all permutations in a
                                single request. Default: 24.
    SCHUBMULT_ENABLE_MULT       Set to "1" to enable the --mult polynomial
                                factor (passes through sympify; treat as
                                untrusted). Default: disabled.
"""

import io
import multiprocessing as mp
import os
import re
import shlex
import traceback
from contextlib import redirect_stderr, redirect_stdout

from flask import Flask, jsonify, render_template, request

app = Flask(__name__)

PERM_TOKEN_RE = re.compile(r"^-?\d+$")

# ---------- config ----------
ALLOWED_ORIGINS = [
    o.strip() for o in os.environ.get("SCHUBMULT_ALLOWED_ORIGINS", "").split(",")
    if o.strip()
]
COMPUTE_TIMEOUT = float(os.environ.get("SCHUBMULT_COMPUTE_TIMEOUT", "8"))
MAX_PERM_LENGTH = int(os.environ.get("SCHUBMULT_MAX_PERM_LENGTH", "24"))
ENABLE_MULT = os.environ.get("SCHUBMULT_ENABLE_MULT", "0") == "1"
MAX_INT_VALUE = 64  # reject permutation entries above this


def _tokenize(s: str) -> list[str]:
    """Split user input into shell-like tokens. Accepts spaces, commas, '-'."""
    s = s.strip()
    if not s:
        return []
    s = s.replace(",", " ")
    return shlex.split(s)


def _build_argv(prog: str, perms_raw: str, *, ascode: bool, coprod: bool,
                display_positive: bool, mult: str | None) -> list[str]:
    """Build sys.argv-like list to pass to the script's main(). Validates input."""
    argv: list[str] = [prog]
    tokens = _tokenize(perms_raw)
    if not tokens:
        raise ValueError("No permutation supplied.")

    int_count = 0
    for tok in tokens:
        if tok == "-":
            continue
        if not PERM_TOKEN_RE.match(tok):
            raise ValueError(f"Invalid token {tok!r}: expected integers and '-' separators.")
        v = int(tok)
        if not (-MAX_INT_VALUE <= v <= MAX_INT_VALUE):
            raise ValueError(f"Permutation entry {v} out of allowed range "
                             f"[-{MAX_INT_VALUE}, {MAX_INT_VALUE}].")
        int_count += 1

    if int_count > MAX_PERM_LENGTH:
        raise ValueError(f"Too many integers ({int_count}); limit is "
                         f"{MAX_PERM_LENGTH} per request.")

    if "-" not in tokens and not coprod:
        raise ValueError("Need at least two permutations separated by '-'. "
                         "Example:  3 1 2 - 2 1 3")

    if coprod:
        argv.append("--coprod")
    if ascode:
        argv.append("--code")
    if display_positive:
        argv.append("--display-positive")
    if mult:
        if not ENABLE_MULT:
            raise ValueError("The --mult option is disabled on this server.")
        if len(mult) > 200:
            raise ValueError("--mult expression too long (max 200 chars).")
        argv.extend(["--mult", mult])
    argv.extend(tokens)
    return argv


def _worker(flavor: str, argv: list[str], q) -> None:
    """Subprocess entrypoint: run the script's main() and ship back its stdio."""
    out, err = io.StringIO(), io.StringIO()
    try:
        if flavor == "py":
            from schubmult._scripts import schubmult_py as mod
        else:
            from schubmult._scripts import schubmult_double as mod
        with redirect_stdout(out), redirect_stderr(err):
            try:
                mod.main(argv)
            except SystemExit:
                pass
            except Exception:
                err.write(traceback.format_exc())
    finally:
        q.put((out.getvalue(), err.getvalue()))


def _run_script(flavor: str, argv: list[str]) -> tuple[str, str, bool]:
    """Run the script in a child process with a hard timeout.

    Returns (stdout, stderr, timed_out).
    """
    ctx = mp.get_context("spawn")
    q = ctx.Queue()
    p = ctx.Process(target=_worker, args=(flavor, argv, q), daemon=True)
    p.start()
    p.join(COMPUTE_TIMEOUT)
    if p.is_alive():
        p.terminate()
        p.join(1.0)
        if p.is_alive():
            p.kill()
        return ("", f"Computation exceeded {COMPUTE_TIMEOUT}s timeout.\n", True)
    if q.empty():
        return ("", "Worker exited without producing output.\n", False)
    out, err = q.get()
    return (out, err, False)


# ---------- routes ----------

def _csp_header() -> str:
    if ALLOWED_ORIGINS:
        ancestors = "'self' " + " ".join(ALLOWED_ORIGINS)
    else:
        ancestors = "'self'"
    return (
        "default-src 'self'; "
        "script-src 'self' 'unsafe-inline'; "
        "style-src 'self' 'unsafe-inline'; "
        f"frame-ancestors {ancestors};"
    )


@app.after_request
def _add_security_headers(resp):
    resp.headers["Content-Security-Policy"] = _csp_header()
    resp.headers["X-Content-Type-Options"] = "nosniff"
    resp.headers["Referrer-Policy"] = "no-referrer"
    return resp


@app.route("/")
def index():
    return render_template("index.html", embed=False, enable_mult=ENABLE_MULT)


@app.route("/embed")
def embed():
    return render_template("index.html", embed=True, enable_mult=ENABLE_MULT)


@app.route("/api/compute", methods=["POST"])
def compute():
    data = request.get_json(force=True, silent=True) or {}
    flavor = data.get("flavor", "py")
    perms_raw = data.get("perms", "")
    ascode = bool(data.get("ascode", False))
    coprod = bool(data.get("coprod", False))
    display_positive = bool(data.get("display_positive", False))
    mult = (data.get("mult") or "").strip() or None

    if flavor == "py":
        prog = "schubmult_py"
        display_positive = False
    elif flavor == "double":
        prog = "schubmult_double"
    else:
        return jsonify({"ok": False, "error": f"Unknown flavor {flavor!r}"}), 400

    try:
        argv = _build_argv(prog, perms_raw, ascode=ascode, coprod=coprod,
                           display_positive=display_positive, mult=mult)
    except ValueError as e:
        return jsonify({"ok": False, "error": str(e)}), 400

    stdout, stderr, timed_out = _run_script(flavor, argv)
    return jsonify({
        "ok": not timed_out,
        "argv": argv,
        "stdout": stdout,
        "stderr": stderr,
        "timed_out": timed_out,
    })


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False)
