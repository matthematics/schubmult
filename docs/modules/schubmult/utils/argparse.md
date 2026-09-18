<a id="schubmult.utils.argparse"></a>

# schubmult.utils.argparse

Shared command-line parsing for the ``schubmult_*`` scripts (`schub_argparse`).

<a id="schubmult.utils.argparse.schub_argparse"></a>

#### schub\_argparse

```python
def schub_argparse(prog_name,
                   description,
                   argv,
                   quantum=False,
                   yz=False,
                   coprod=True)
```

Parse the common CLI of the ``schubmult_*`` scripts and return ``(args, formatter)``.

Permutations are given as space-separated integers separated by ``-`` (or Lehmer codes with
``--code``); ``quantum`` adds the ``--parabolic`` options, ``yz`` the double-variable and
``--display-positive`` options, ``coprod`` the ``--coprod`` mode. ``--display-mode`` selects the
``formatter`` (a callable rendering expressions as LaTeX, pretty, basic, sympy, or ``None`` for
raw). Hidden ``-g`` dumps the parsed arguments to a JSON file for the script test suite and
exits. Also initializes SymPy printing and logging.

