# schubmult web API

This document describes the HTTP API exposed by the schubmult web wrapper. The
same endpoint powers the `[schubmult]` WordPress widget; it is also usable
directly from `curl`, JavaScript, or any HTTP client.

**Base URL:** `https://schubmult.pythonanywhere.com`
(replace with your own host if you self-deploy)

All endpoints are CORS-friendly; the only cross-origin restriction is on
`<iframe>` embedding (controlled by the `Content-Security-Policy:
frame-ancestors` header on the server).

---

## Endpoints

| Method | Path           | Purpose                                    |
|--------|----------------|--------------------------------------------|
| `GET`  | `/`            | Full HTML page with the widget.            |
| `GET`  | `/embed`       | Chrome-less version of the widget for `<iframe>`. |
| `POST` | `/api/compute` | Run a Schubert polynomial computation.     |

---

## `POST /api/compute`

### Request

`Content-Type: application/json`. Body fields:

| Field              | Type    | Default | Description                                                             |
|--------------------|---------|---------|-------------------------------------------------------------------------|
| `flavor`           | string  | `"py"`  | Which kernel: `py`, `groth`, `double`, `q`, or `q_double`. See below.   |
| `perms`            | string  | `""`    | One-line permutations separated by `-`. E.g. `"3 1 2 - 2 1 3"`. With `ascode=true`, these are Lehmer codes instead. |
| `ascode`           | boolean | `false` | Treat each token sequence as a Lehmer code rather than a permutation.   |
| `coprod`           | boolean | `false` | Compute the coproduct of a single permutation along a split index list. Only honoured for `py`/`double`. With `coprod=true`, `perms` is one permutation followed by `-` followed by the split indices, e.g. `"1 4 2 3 - 1 2"`. |
| `display_positive` | boolean | `false` | (`double`/`q_double` only) Solve a MILP to display coefficients in the root-positive form. Slower. |
| `mixed_var`        | boolean | `false` | (`double`/`q_double` only) Use two variable sets `y`, `z` instead of just `y`. |
| `parabolic`        | string  | `""`    | (`q`/`q_double` only) Space-separated positive integers giving the simple-reflection generators of the parabolic subgroup, e.g. `"1 3"`. Empty = non-parabolic. |
| `mult`             | string  | `""`    | (Disabled by default on this host for security.) Polynomial factor parsed by SymPy. |

### Permutation input rules (when `ascode=false`)

Each `-`-separated chunk must be a **bona fide permutation of `{1,…,n}`**:
distinct positive integers forming the initial interval. The server will
reject input like `1 1 2` or `3 1 4` with an explanatory `400`. If you wanted
a Lehmer code (which has no such constraint), set `ascode=true`.

### Flavor / option compatibility matrix

| Option             | `py` | `groth` | `double` | `q` | `q_double` |
|--------------------|:----:|:-------:|:--------:|:---:|:----------:|
| `coprod`           |  ✓   |    —    |    ✓     |  —  |     —      |
| `display_positive` |  —   |    —    |    ✓     |  —  |     ✓      |
| `mixed_var`        |  —   |    —    |    ✓     |  —  |     ✓      |
| `parabolic`        |  —   |    —    |    —     |  ✓  |     ✓      |

Unsupported options are silently ignored for that flavor.

### Limits

- **`SCHUBMULT_MAX_PERM_LENGTH`** integers per request (default 24, may be
  lower on this host). Each integer counts; `-` separators don't.
- **Per-entry range:** integers must be in `[-64, 64]`.
- **Compute timeout:** `SCHUBMULT_COMPUTE_TIMEOUT` seconds (default 8). On
  timeout, the response has `ok: false` and `timed_out: true`.

### Response

`Content-Type: application/json`:

```json
{
  "ok": true,
  "argv": ["schubmult_py", "3", "1", "2", "-", "2", "1", "3"],
  "stdout": "1  (4, 1, 2, 3)\n",
  "stderr": "",
  "timed_out": false
}
```

| Field       | Description                                                                |
|-------------|----------------------------------------------------------------------------|
| `ok`        | `true` if the script ran to completion, `false` on validation error or timeout. |
| `argv`      | The full CLI invocation that was executed (useful for reproducing locally). |
| `stdout`    | Captured standard output of the script. This is the actual answer.         |
| `stderr`    | Captured standard error (warnings, debug messages).                        |
| `timed_out` | `true` if the request was killed by the server-side timeout.               |

On a validation error the response is HTTP `400` with:

```json
{ "ok": false, "error": "Permutation #1 (1 1 2) is not a permutation of 1..3 ..." }
```

### Stdout format

The script prints one term per line:

```
<coefficient>  <permutation as Python tuple>
```

For example `schubmult_py 3 1 2 - 2 1 3` prints:

```
1  (4, 1, 2, 3)
```

For `double` and `q_double` the coefficient is a polynomial in `y` (or in `y`,
`z` if `mixed_var=true`); for `q` and `q_double` the permutation may be
preceded by quantum-deformation factors `q_i`.

---

## Examples

### Ordinary Schubert product
```bash
curl -X POST https://schubmult.pythonanywhere.com/api/compute \
  -H 'Content-Type: application/json' \
  -d '{"flavor":"py","perms":"3 1 2 - 2 1 3"}'
```

### Grothendieck product
```bash
curl -X POST https://schubmult.pythonanywhere.com/api/compute \
  -H 'Content-Type: application/json' \
  -d '{"flavor":"groth","perms":"3 1 2 - 2 1 3"}'
```

### Double Schubert product, mixed variables, root-positive display
```bash
curl -X POST https://schubmult.pythonanywhere.com/api/compute \
  -H 'Content-Type: application/json' \
  -d '{"flavor":"double","perms":"3 1 2 - 2 1 3",
       "mixed_var":true,"display_positive":true}'
```

### Quantum Schubert product with parabolic subgroup
```bash
curl -X POST https://schubmult.pythonanywhere.com/api/compute \
  -H 'Content-Type: application/json' \
  -d '{"flavor":"q","perms":"3 1 2 - 2 1 3","parabolic":"1"}'
```

### Lehmer-code input
```bash
curl -X POST https://schubmult.pythonanywhere.com/api/compute \
  -H 'Content-Type: application/json' \
  -d '{"flavor":"py","perms":"2 0 - 1 0","ascode":true}'
```

### Coproduct (split a permutation along an index list)
```bash
curl -X POST https://schubmult.pythonanywhere.com/api/compute \
  -H 'Content-Type: application/json' \
  -d '{"flavor":"py","perms":"1 4 2 3 - 1 2","coprod":true}'
```

### Browser / JavaScript
```js
const r = await fetch('https://schubmult.pythonanywhere.com/api/compute', {
  method: 'POST',
  headers: {'Content-Type': 'application/json'},
  body: JSON.stringify({flavor: 'py', perms: '3 1 2 - 2 1 3'}),
});
const j = await r.json();
console.log(j.ok ? j.stdout : j.error);
```

### Python
```python
import requests
r = requests.post(
    "https://schubmult.pythonanywhere.com/api/compute",
    json={"flavor": "double", "perms": "3 1 2 - 2 1 3", "mixed_var": True},
)
data = r.json()
print(data["stdout"])
```

---

## Error reference

| HTTP | Cause                                                | Example `error` text                                 |
|------|------------------------------------------------------|------------------------------------------------------|
| 400  | Unknown flavor                                       | `Unknown flavor 'foo'`                               |
| 400  | Empty input                                          | `No permutation supplied.`                           |
| 400  | Non-integer token                                    | `Invalid token 'x': expected integers and '-' separators.` |
| 400  | Out-of-range integer                                 | `Permutation entry 99 out of allowed range [-64, 64].` |
| 400  | Too many integers                                    | `Too many integers (50); limit is 24 per request.`   |
| 400  | Missing `-` separator (and not `coprod`)             | `Need at least two permutations separated by '-'.`   |
| 400  | Tokens are not a real permutation (and not `ascode`) | `Permutation #1 (3 1 4) is not a permutation of 1..3 ...` |
| 400  | Bad parabolic generator                              | `Parabolic generator 0 out of range [1, 64].`        |
| 200, `timed_out=true` | Computation exceeded the timeout       | `stderr` may contain partial diagnostics.            |

---

## Embedding

Use the bundled WordPress plugin (`schubmult-embed.php`) and the shortcode

```
[schubmult]
```

The iframe auto-resizes to fit the widget. See `web/DEPLOY.md` for the full
deployment guide.

---

## Source

The web wrapper lives in [`web/app.py`](https://github.com/matthematics/schubmult/blob/main/web/app.py)
of the [schubmult repository](https://github.com/matthematics/schubmult).
The underlying CLI tools `schubmult_py`, `grothmult_py`, `schubmult_double`,
`schubmult_q`, `schubmult_q_double` are documented in the project README.
