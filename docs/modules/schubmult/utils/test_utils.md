<a id="schubmult.utils.test_utils"></a>

# schubmult.utils.test\_utils

Helpers for the test suite: locating JSON test data and inspecting SymPy/SymEngine expression trees.

<a id="schubmult.utils.test_utils.generate_all"></a>

#### generate\_all

```python
def generate_all(module, filename)
```

Print an import block and ``__all__`` list for the public names defined in ``filename`` (dev helper).

<a id="schubmult.utils.test_utils.get_json"></a>

#### get\_json

```python
def get_json(file: str)
```

Load ``<file>.json`` from the test data directory.

<a id="schubmult.utils.test_utils.load_json_test_names"></a>

#### load\_json\_test\_names

```python
def load_json_test_names(this_dir)
```

Names (without ``.json``) of all test case files in the data subdirectory ``this_dir``.

<a id="schubmult.utils.test_utils.print_args"></a>

#### print\_args

```python
def print_args(poly)
```

Nested string of the argument types of an expression tree (for debugging printing issues).

<a id="schubmult.utils.test_utils.sympify_args"></a>

#### sympify\_args

```python
def sympify_args(poly)
```

Convert a SymPy expression to SymEngine, recursing into ``Mul``/``Pow``/``Add`` when direct
conversion fails.

