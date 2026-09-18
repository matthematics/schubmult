<a id="schubmult.symbolic.symmetric_polynomials.qelem_sym"></a>

# schubmult.symbolic.symmetric\_polynomials.qelem\_sym

Quantum factorial elementary symmetric polynomials ``E_q(p, k, xvars, yvars)`` as SymPy atoms.

The quantum deformation adds ``q_i`` terms for adjacent pairs ``x_i, x_{i+1}`` of the *positional*
generators (positions taken from ``x_var``), so the variable set need not be an initial segment.
``expand_func`` evaluates via `elem_sym_positional_poly_q`. Alias: `QFactorialElemSym`.

<a id="schubmult.symbolic.symmetric_polynomials.qelem_sym.E_q"></a>

## E\_q Objects

```python
class E_q(ElemSym_base)
```

Quantum factorial elementary symmetric atom; same canonicalization as `E`, plus ``x_var``/``q_var``
generating sets fixing the positions of the generators and the ``q`` parameters.

<a id="schubmult.symbolic.symmetric_polynomials.qelem_sym.elem_sym_positional_poly_q"></a>

#### elem\_sym\_positional\_poly\_q

```python
def elem_sym_positional_poly_q(p,
                               k,
                               varl1,
                               varl2,
                               x_var=GeneratingSet("x"),
                               q_var=GeneratingSet("q"))
```

Quantum factorial elementary symmetric polynomial of degree ``p`` in the generators ``varl1[:k]``.

Recursion on the last generator ``x_l``: the classical two terms plus, when ``x_{l-1}`` (resp.
``x_{l+1}``) is also among the generators, ``q_{l-1}`` (resp. ``q_l``) times the degree ``p - 2``
polynomial with that pair removed.

