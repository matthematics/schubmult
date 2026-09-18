<a id="schubmult.utils.parsing"></a>

# schubmult.utils.parsing

Parsing of coefficient expressions from the command line.

<a id="schubmult.utils.parsing.parse_coeff"></a>

#### parse\_coeff

```python
def parse_coeff(coeff_str, latex=False)
```

Sympify ``coeff_str`` and map any symbol ``name_i`` to ``GeneratingSet(name)[i]`` so the result
uses the package's interned variables. LaTeX input is not yet supported (returns ``None``).

