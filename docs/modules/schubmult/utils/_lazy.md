<a id="schubmult.utils._lazy"></a>

# schubmult.utils.\_lazy

`LazyAttr`: a stand-in for a module attribute that imports the module on first real use.

<a id="schubmult.utils._lazy.LazyAttr"></a>

## LazyAttr Objects

```python
class LazyAttr()
```

Proxy for ``getattr(import_module(modname), attr)``, resolved on first call, attribute access,
or use as a base class. ``isinstance``/``issubclass`` checks resolve it only once ``modname`` has
been imported, since no instance of a class can exist before its module is loaded.

<a id="schubmult.utils._lazy.lazy_from"></a>

#### lazy\_from

```python
def lazy_from(modname, *names)
```

``LazyAttr`` proxies for ``from modname import *names``.

