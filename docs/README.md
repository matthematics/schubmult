# schubmult API docs

Generated from source docstrings with [pydoc-markdown](https://github.com/NiklasRosenstein/pydoc-markdown).

- [API.md](API.md): the whole package as one consolidated reference.
- [modules/](modules/README.md): one page per module, mirroring `src/schubmult`.

## Regenerating

```bash
python docs/generate_docs.py
```

This reads the shared loader/processor/renderer settings from
[pydoc-markdown.yml](../pydoc-markdown.yml) at the repo root. `schubmult._scripts.unlinted`
(research scratch space) is excluded; everything else under `src/schubmult` is included,
and undocumented members are dropped (`documented_only: true`).

To render a single module manually:

```bash
pydoc-markdown pydoc-markdown.yml -m schubmult.mult.single
```
