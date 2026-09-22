<a id="schubmult.mult.groth_quantum"></a>

# schubmult.mult.groth\_quantum

Multiplication kernel for (single) quantum beta-Grothendieck polynomials.

``grothmult_q`` computes ``G^q_u(x) * G^q_v(x)`` in the ``G^q`` basis: the ``y = z = 0``
specialization of ``groth_quantum_double.grothmult_q_double``, by the same route as
``groth.grothmult_py`` -- expand ``G_v`` into Schubert polynomials, push each ``S_{v'}``
through the strict-theta v-path layers, and use the quantum K-Pieri support
``quantum_pieri_chains`` with the binomial closed form ``groth_elem_sym_coeff`` evaluated at
the quantum chain length and weighted by ``q^D``.

``G^q_v = Q(G_v)`` is the Lenart--Maeno quantization (``groth_quantum_double.lm_quantize``);
at ``beta = -1`` these are the quantum Grothendieck polynomials of Lenart--Maeno with
``Q_j = q_j``, and the ``e_p`` Pieri rule ``grothmult_q_pieri`` is then equivalent to the
Naito--Sagaki quantum K Pieri theorem (arXiv:2211.01578): the sign ``(-1)^{len - p}`` and the
marking count ```Mark``` there are the ``beta^{len - m} (-beta)^{L - p'} binom(L, p')`` here,
summed over the ``e_p <-> G_{c[k,p]}`` triangular change of basis.  Conjectural in general;
see ``groth_quantum_double`` for the evidence.

<a id="schubmult.mult.groth_quantum.grothmult_q_pieri"></a>

#### grothmult\_q\_pieri

```python
def grothmult_q_pieri(coeff_dict, p, k, beta=None, q_var=None)
```

Multiply ``sum_u coeff_u G^q_u(x)`` by ``Q(e_p(x_1..x_k))``.

Coefficient of ``G^q_w``: ``q^D * groth_elem_sym_coeff(k, u, w, k - p, beta, length)``
over the quantum K-Pieri support; ``p = k`` is the non-equivariant quantum top block.

<a id="schubmult.mult.groth_quantum.grothmult_q"></a>

#### grothmult\_q

```python
def grothmult_q(perm_dict, v, beta=None, q_var=None)
```

Multiply (single) quantum Grothendieck polynomials, mirroring ``schubmult_q``.

Returns the expansion of ``sum_u coeff_u G^q_u(x) * G^q_v(x)`` in the basis ``{G^q_w(x)}``
as ``{w: coeff_w}``, polynomial in ``beta`` and ``q``.  ``G_v`` is expanded into Schubert
polynomials by ``groth_elem_as_schub_dict`` (quantization is linear over ``beta``) and
each ``S_{v'}`` pushed through ``_qgroth_schub_vpath_mul``.  ``q = 0`` is ``grothmult_py``;
``beta = 0`` is ``schubmult_q``.

<a id="schubmult.mult.groth_quantum.grothmult_q_dict"></a>

#### grothmult\_q\_dict

```python
def grothmult_q_dict(perm_dict1, perm_dict2, beta=None, q_var=None)
```

Product of two coefficient dicts: ``sum_v coeff2_v grothmult_q(perm_dict1, v, ...)``.

