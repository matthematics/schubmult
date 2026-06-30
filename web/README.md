# schubmult web wrapper

A small Flask app that exposes `schubmult_py`, `grothmult_py`, and `schubmult_double` through a
web form. Designed to be embeddable in another page via `<iframe>`.

## Run

```bash
conda activate schubmult_312
pip install flask
python web/app.py
```

Then visit:

- `http://localhost:5000/` — full page
- `http://localhost:5000/embed` — bare widget (no `<html>` chrome)

## Embed

```html
<iframe src="https://your-host/embed" width="800" height="600"
        style="border: 1px solid #ccc;"></iframe>
```

## API

`POST /api/compute` with JSON body:

```json
{
  "flavor": "py" | "groth" | "double" | "q" | "q_double",
  "perms": "3 1 2 - 2 1 3",
  "ascode": false,
  "coprod": false,
  "display_positive": false,
  "mult": ""
}
```

Returns `{"ok": true, "argv": [...], "stdout": "...", "stderr": "..."}`.

## Notes

- The app calls each script's `main(argv)` directly (no subprocess), so it
  shares the importing process's Python environment.
- For production, wrap in a real WSGI server (gunicorn / waitress) and put a
  rate limiter in front; permutation arguments are validated but the
  `--mult` polynomial expression is passed to `sympify` and should be treated
  as untrusted input on a public deployment.
