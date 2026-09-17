<a id="schubmult"></a>

# schubmult

<a id="schubmult.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name: str)
```

Lazy import exported names. This allows `import schubmult` to succeed
even if some optional dependencies are missing.

<a id="schubmult.__dir__"></a>

#### \_\_dir\_\_

```python
def __dir__()
```

Include lazily-exported names in dir()

