# schubmult API docs

Generated from source docstrings with [pydoc-markdown](https://github.com/NiklasRosenstein/pydoc-markdown)
and published as a static site with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/).

- [index.md](index.md): the site homepage.
- [API.md](API.md): the whole package as one consolidated reference (not part of the built site; see `mkdocs.yml`'s `exclude_docs`).
- [modules/](modules/README.md): one page per module, mirroring `src/schubmult`. This is what the site's navigation is generated from.

## Hosting

Pushes to `main` that touch `src/`, `docs/`, or `mkdocs.yml` trigger
[.github/workflows/docs.yml](../.github/workflows/docs.yml), which regenerates the API
reference, builds the MkDocs site, and deploys it to GitHub Pages
(`https://matthematics.github.io/schubmult/`). Requires the repo's
**Settings -> Pages -> Source** set to "GitHub Actions" (one-time setup).

## Regenerating locally

```bash
pip install -e .[docs]
python docs/generate_docs.py   # refresh docs/modules/ and docs/API.md from docstrings
mkdocs serve                   # live preview at http://127.0.0.1:8000
mkdocs build --site-dir site   # static build, matches what CI deploys
```

This reads the shared loader/processor/renderer settings from
[pydoc-markdown.yml](../pydoc-markdown.yml) at the repo root. `schubmult._scripts.unlinted`
(research scratch space) is excluded; everything else under `src/schubmult` is included,
and undocumented members are dropped (`documented_only: true`).

To render a single module manually:

```bash
pydoc-markdown pydoc-markdown.yml -m schubmult.mult.single
```
