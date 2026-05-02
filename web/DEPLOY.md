# Deploying schubmult web wrapper to PythonAnywhere

This is the simplest path that works with a GoDaddy-hosted WordPress site.
PythonAnywhere has a free tier with full Python and a persistent process; the
WordPress side just embeds an `<iframe>` via the bundled plugin.

## 1. Sign up

Create a free account at <https://www.pythonanywhere.com/>. The free plan is
sufficient for small permutations.

## 2. Upload the code

In a PythonAnywhere **Bash console**:

```bash
git clone --depth=1 https://github.com/matthematics/schubmult.git ~/schubmult
cd ~/schubmult
python3.12 -m venv ~/.venvs/schubmult
source ~/.venvs/schubmult/bin/activate
pip install --upgrade pip
pip install -e . flask
pip cache purge          # IMPORTANT: pip caches ~400MB of wheels
```

Expected disk usage after this:

| What | Size |
|---|---|
| Virtualenv (sympy + symengine + numpy + flask + schubmult) | ~280 MB |
| Source clone (`--depth=1`) | ~10 MB |
| `pip cache purge` saves | ~400 MB |
| **Total** | **~290 MB** |

**Disk-saving variants:**

- *Minimal* (~85 MB total) — skip the heavy compiled deps; only `schubmult_py`
  works, and slower than with symengine:
  ```bash
  pip install sympy cachetools sortedcontainers flask
  pip install -e . --no-deps
  ```
- Drop `pulp` (~36 MB) — only needed for `--display-positive`. It's not in the
  default install anyway unless you add the `[milp]` extra.
- Drop `matplotlib`, `Pillow`, `cryptography` — not used by the web wrapper.

If `symengine` fails to build on the free tier, fall back to:

```bash
SCHUBMULT_NO_SYMENGINE=1 pip install -e .
```

(Sympy alone is enough for `schubmult_py`; double Schubert polynomials may be
slower without symengine.)

## 3. Create the web app

In the PythonAnywhere dashboard:

1. **Web → Add a new web app**.
2. Choose **Manual configuration**, Python 3.12.
3. After it's created, scroll to **Code** and set:
   - **Source code:** `/home/YOURUSER/schubmult`
   - **Working directory:** `/home/YOURUSER/schubmult`
4. Under **Virtualenv**, enter `/home/YOURUSER/.venvs/schubmult`.
5. Click the **WSGI configuration file** link and replace its contents with:

```python
import os, sys

# Your WordPress origin(s):
os.environ['SCHUBMULT_ALLOWED_ORIGINS'] = 'https://your-wordpress-site.com,https://www.your-wordpress-site.com'
os.environ['SCHUBMULT_COMPUTE_TIMEOUT'] = '8'
os.environ['SCHUBMULT_MAX_PERM_LENGTH'] = '24'
# Leave SCHUBMULT_ENABLE_MULT unset (i.e. disabled) for public deployments.

project_home = '/home/YOURUSER/schubmult/web'
if project_home not in sys.path:
    sys.path.insert(0, project_home)

from app import app as application  # noqa: E402
```

6. **Reload** the web app (green button at the top of the Web tab).

Visit `https://YOURUSER.pythonanywhere.com/embed` to confirm the widget loads.

## 4. Connect it to WordPress

1. On your WordPress site, **Plugins → Add New → Upload Plugin**.
2. Upload `wordpress-plugin/schubmult-embed.php` (zip it first if WP requires
   a directory archive).
3. Activate. Go to **Settings → Schubmult Embed**.
4. Paste `https://YOURUSER.pythonanywhere.com/embed` and save.
5. In any page/post add the shortcode:

```
[schubmult height="640"]
```

Done — visitors can now compute Schubert products on your WordPress site.

## Optional: using Render / Fly / your own VPS

The same `web/app.py` works with any WSGI host. Generic recipe:

```bash
pip install gunicorn
gunicorn -w 2 -b 0.0.0.0:8000 web.app:app
```

Set the same environment variables (`SCHUBMULT_ALLOWED_ORIGINS`, etc.) and
point the WordPress plugin at `https://your-host/embed`.

## Hardening notes

The bundled `web/app.py` already:

- Validates permutation tokens (integers only, bounded range).
- Rejects requests above `SCHUBMULT_MAX_PERM_LENGTH` integers (default 24).
- Runs each computation in a child process with a hard
  `SCHUBMULT_COMPUTE_TIMEOUT` (default 8s). The child is `terminate()`d /
  `kill()`ed if it overruns.
- Disables the `--mult` polynomial input by default (it goes through
  `sympify`). Only enable it on trusted, authenticated deployments by setting
  `SCHUBMULT_ENABLE_MULT=1`.
- Sets `Content-Security-Policy: frame-ancestors` to your allowed origins.
- Sets `X-Content-Type-Options: nosniff`, `Referrer-Policy: no-referrer`.

For higher-traffic deployments, put a rate limiter (e.g.
[`flask-limiter`](https://flask-limiter.readthedocs.io/)) in front of
`/api/compute`.
