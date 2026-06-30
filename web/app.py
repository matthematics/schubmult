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
import logging
import logging.handlers
import multiprocessing as mp
import os
import re
import shlex
import traceback
from contextlib import redirect_stderr, redirect_stdout

from flask import Flask, jsonify, render_template, request

app = Flask(__name__)

# ---------- access logging ----------
# Logs each /api/compute request with the exact command line that would be
# invoked (argv), but NEVER the output. Path configurable via env var
# SCHUBMULT_ACCESS_LOG; defaults to web/access.log next to this file.
# Set SCHUBMULT_ACCESS_LOG=- to log to stderr instead.
_DEFAULT_LOG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "access.log")
_ACCESS_LOG_PATH = os.environ.get("SCHUBMULT_ACCESS_LOG", _DEFAULT_LOG_PATH)

access_logger = logging.getLogger("schubmult.web.access")
access_logger.setLevel(logging.INFO)
access_logger.propagate = False
if not access_logger.handlers:
    _fmt = logging.Formatter("%(asctime)s %(message)s",
                             datefmt="%Y-%m-%dT%H:%M:%S%z")
    if _ACCESS_LOG_PATH == "-":
        _h: logging.Handler = logging.StreamHandler()
    else:
        try:
            _h = logging.handlers.RotatingFileHandler(
                _ACCESS_LOG_PATH, maxBytes=2_000_000, backupCount=5,
                encoding="utf-8")
        except OSError:
            _h = logging.StreamHandler()
    _h.setFormatter(_fmt)
    access_logger.addHandler(_h)


def _client_ip() -> str:
    xff = request.headers.get("X-Forwarded-For", "")
    if xff:
        return xff.split(",")[0].strip()
    return request.remote_addr or "-"


def _log_access(flavor: str, argv: list[str] | None, *, status: str,
                error: str | None = None) -> None:
    cmdline = shlex.join(argv) if argv else "-"
    ua = request.headers.get("User-Agent", "-").replace("\n", " ")[:200]
    msg = (f"ip={_client_ip()} flavor={flavor} status={status} "
           f"cmd={cmdline!r} ua={ua!r}")
    if error:
        msg += f" error={error!r}"
    access_logger.info(msg)

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
                display_positive: bool, mixed_var: bool, parabolic: str | None,
                mult: str | None) -> list[str]:
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
        if v < 0:
            raise ValueError(f"Negative entry {v} not allowed.")
        if v > MAX_INT_VALUE:
            raise ValueError(f"Permutation entry {v} out of allowed range "
                             f"[0, {MAX_INT_VALUE}].")
        if not ascode and v < 1:
            raise ValueError(f"Permutation entry {v} not allowed: "
                             f"one-line entries must be >= 1. "
                             f"(Use --code if you meant a Lehmer code.)")
        int_count += 1

    if int_count > MAX_PERM_LENGTH:
        raise ValueError(f"Too many integers ({int_count}); limit is "
                         f"{MAX_PERM_LENGTH} per request.")

    if "-" not in tokens and not coprod:
        raise ValueError("Need at least two permutations separated by '-'. "
                         "Example:  3 1 4 2 - 2 5 1 4 3")

    # When tokens are interpreted as one-line permutations (not Lehmer codes),
    # each '-'-separated chunk must be a permutation of {1,...,n}: distinct
    # positive integers forming an initial interval.
    if not ascode:
        chunks: list[list[str]] = [[]]
        for tok in tokens:
            if tok == "-":
                chunks.append([])
            else:
                chunks[-1].append(tok)
        for idx, chunk in enumerate(chunks, start=1):
            if not chunk:
                raise ValueError(f"Permutation #{idx} is empty.")
            try:
                vals = [int(t) for t in chunk]
            except ValueError as e:  # pragma: no cover -- already validated
                raise ValueError(str(e)) from None
            n = len(vals)
            if sorted(vals) != list(range(1, n + 1)):
                raise ValueError(
                    f"Permutation #{idx} ({' '.join(chunk)}) is not a permutation "
                    f"of 1..{n}: entries must be distinct positive integers forming "
                    f"the initial interval 1..{n}. "
                    f"(Use --code if you meant a Lehmer code.)"
                )

    if coprod:
        argv.append("--coprod")
    if ascode:
        argv.append("--code")
    if display_positive:
        argv.append("--display-positive")
    if mixed_var:
        argv.append("--mixed-var")
    # Validate parabolic but do NOT append it yet; --parabolic uses nargs="+" and
    # would greedily swallow the perm tokens. We append it *after* the perms.
    ptoks: list[str] = []
    if parabolic:
        ptoks = parabolic.replace(",", " ").split()
        for t in ptoks:
            if not t.isdigit():
                raise ValueError(f"Invalid parabolic generator {t!r}: expected positive integers.")
            if int(t) < 1 or int(t) > MAX_INT_VALUE:
                raise ValueError(f"Parabolic generator {t} out of range [1, {MAX_INT_VALUE}].")
    mult_tail: list[str] = []
    if mult:
        if not ENABLE_MULT:
            raise ValueError("The --mult option is disabled on this server.")
        if len(mult) > 200:
            raise ValueError("--mult expression too long (max 200 chars).")
        mult_tail = ["--mult", mult]
    argv.extend(tokens)
    if ptoks:
        argv.append("--parabolic")
        argv.extend(ptoks)
    argv.extend(mult_tail)
    return argv


def _worker(flavor: str, argv: list[str], q) -> None:
    """Subprocess entrypoint: run the script's main() and ship back its stdio."""
    out, err = io.StringIO(), io.StringIO()
    try:
        if flavor == "py":
            from schubmult._scripts import schubmult_py as mod
        elif flavor == "groth":
            from schubmult._scripts import grothmult_py as mod
        elif flavor == "double":
            from schubmult._scripts import schubmult_double as mod
        elif flavor == "q":
            from schubmult._scripts import schubmult_q as mod
        elif flavor == "q_double":
            from schubmult._scripts import schubmult_q_double as mod
        else:
            raise ValueError(f"Unknown flavor {flavor!r}")
        with redirect_stdout(out), redirect_stderr(err):
            try:
                mod.main(argv)
            except SystemExit:
                pass
            except Exception:
                err.write(traceback.format_exc())
    finally:
        q.put((out.getvalue(), err.getvalue()))


def _run_inline(flavor: str, argv: list[str]) -> tuple[str, str, bool]:
    """Fallback: run in-process (no timeout). Used when multiprocessing fails."""
    out, err = io.StringIO(), io.StringIO()
    try:
        if flavor == "py":
            from schubmult._scripts import schubmult_py as mod
        elif flavor == "groth":
            from schubmult._scripts import grothmult_py as mod
        elif flavor == "double":
            from schubmult._scripts import schubmult_double as mod
        elif flavor == "q":
            from schubmult._scripts import schubmult_q as mod
        elif flavor == "q_double":
            from schubmult._scripts import schubmult_q_double as mod
        else:
            return ("", f"Unknown flavor {flavor!r}\n", False)
        with redirect_stdout(out), redirect_stderr(err):
            try:
                mod.main(argv)
            except SystemExit:
                pass
            except Exception:
                err.write(traceback.format_exc())
    except Exception:
        err.write(traceback.format_exc())
    return (out.getvalue(), err.getvalue(), False)


def _run_script(flavor: str, argv: list[str]) -> tuple[str, str, bool]:
    """Run the script in a child process with a hard timeout, falling back
    to in-process execution if multiprocessing isn't available.

    Returns (stdout, stderr, timed_out).
    """
    if os.environ.get("SCHUBMULT_DISABLE_SUBPROCESS") == "1":
        return _run_inline(flavor, argv)
    try:
        ctx = mp.get_context("spawn")
        q = ctx.Queue()
        p = ctx.Process(target=_worker, args=(flavor, argv, q), daemon=True)
        p.start()
    except Exception:
        # Fall back to inline if subprocess machinery fails (e.g. some
        # restricted hosts disallow exec/fork).
        return _run_inline(flavor, argv)
    p.join(COMPUTE_TIMEOUT)
    if p.is_alive():
        p.terminate()
        p.join(1.0)
        if p.is_alive():
            p.kill()
        return ("", f"Computation exceeded {COMPUTE_TIMEOUT}s timeout.\n", True)
    if q.empty():
        return ("", "Worker exited without producing output. "
                    "Set SCHUBMULT_DISABLE_SUBPROCESS=1 to run inline.\n", False)
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
    mixed_var = bool(data.get("mixed_var", False))
    parabolic = (data.get("parabolic") or "").strip() or None
    mult = (data.get("mult") or "").strip() or None

    if flavor == "py":
        prog = "schubmult_py"
        display_positive = False
        mixed_var = False
        parabolic = None
    elif flavor == "groth":
        prog = "grothmult_py"
        coprod = False
        display_positive = False
        mixed_var = False
        parabolic = None
    elif flavor == "double":
        prog = "schubmult_double"
        parabolic = None
    elif flavor == "q":
        prog = "schubmult_q"
        display_positive = False
        mixed_var = False
        coprod = False
    elif flavor == "q_double":
        prog = "schubmult_q_double"
        coprod = False
    else:
        return jsonify({"ok": False, "error": f"Unknown flavor {flavor!r}"}), 400

    try:
        argv = _build_argv(prog, perms_raw, ascode=ascode, coprod=coprod,
                           display_positive=display_positive,
                           mixed_var=mixed_var, parabolic=parabolic, mult=mult)
    except ValueError as e:
        _log_access(flavor, None, status="reject", error=str(e))
        return jsonify({"ok": False, "error": str(e)}), 400

    _log_access(flavor, argv, status="run")
    stdout, stderr, timed_out = _run_script(flavor, argv)
    if timed_out:
        _log_access(flavor, argv, status="timeout")
    return jsonify({
        "ok": not timed_out,
        "argv": argv,
        "stdout": stdout,
        "stderr": stderr,
        "timed_out": timed_out,
    })


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False)
