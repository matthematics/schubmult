<a id="schubmult.utils.logging"></a>

# schubmult.utils.logging

Thin wrappers around the standard ``logging`` module.

<a id="schubmult.utils.logging.init_logging"></a>

#### init\_logging

```python
def init_logging(debug=False)
```

Configure root logging at DEBUG (if ``debug``) or ERROR with a timestamped ``file:line`` format.

<a id="schubmult.utils.logging.get_logger"></a>

#### get\_logger

```python
def get_logger(name)
```

``logging.getLogger(name)``.

