# schubmult web wrapper

A small Flask app that exposes the `schubmult_*` and `grothmult_*` CLI scripts (ordinary, double,
quantum and quantum double Schubert and Grothendieck products) through a web form. Designed to be
embeddable in another page via `<iframe>`.

## Run

```bash
conda activate schubmult_312
pip install flask
python web/app.py
```

Then visit:

- `http://localhost:5000/` — full page
- `http://localhost:5000/embed` — bare widget (no `<html>` chrome)

Select **Download output as a text file instead of displaying it** before
computing to save the result as a UTF-8 `.txt` file. When unchecked, results
appear in the result box as before. The command and any errors or warnings
remain visible in either mode. Downloads contain the returned standard output
and have a separate **100 MB (100,000,000 UTF-8 bytes)** limit by default,
configurable with `SCHUBMULT_MAX_DOWNLOAD_BYTES`. Displayed results retain the
existing `SCHUBMULT_MAX_OUTPUT_BYTES` limit (default 1,048,576 characters).
Output exceeding its limit includes a truncation notice, which adds a few bytes
beyond the limit. Errors retain the display limit, and computation timeouts still
apply in either mode.

No result file is stored on the host: output is buffered in server memory and
sent as JSON, then the browser creates the download. Large downloads can consume
several times their size in server and browser memory; this is not streaming.

For **grothmult_double** and **grothmult_q_double**, select **Multiplicative
variables** to calculate at `beta=-1` with coefficient coordinates `1-1/y_i`
(also `1-1/z_i` when mixed variables are selected), using `ChangeOfVariables`.
The checkbox is disabled and cleared for other flavors. It is off by default;
leave it unchecked for the existing symbolic-beta convention. **Simplify
coefficients** cancels the resulting rational coefficients. The coordinate
conversion leaves quantum parameters unchanged and also works with parabolic
projection and output downloads.

## Embed

```html
<iframe src="https://your-host/embed" width="800" height="600"
        style="border: 1px solid #ccc;"></iframe>
```

If you sandbox the iframe, include `allow-downloads` along with the permissions
needed to run the widget. The bundled WordPress plugin includes this permission.

## API

`POST /api/compute` with JSON body:

```json
{
  "flavor": "py" | "groth" | "double" | "groth_double" | "q" | "groth_q" | "q_double" | "groth_q_double",
  "perms": "3 1 2 - 2 1 3",
  "ascode": false,
  "coprod": false,
  "display_positive": false,
  "multiplicative_variables": false,
  "download": false,
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
