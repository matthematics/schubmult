<a id="schubmult.rings.schubert.chevalley"></a>

# schubmult.rings.schubert.chevalley

Lenart--Postnikov :math:`K_T`-Chevalley formula in type :math:`A`.

Implements Theorem 6.1 / Proposition 14.5 of Lenart--Postnikov, *Affine Weyl
groups in K-theory and representation theory* (arXiv:math/0309207)::

    e^lambda . [O_u] = sum_{w, mu} c^{lambda, mu}_{u, w} x^mu [O_w]

    c^{lambda, mu}_{u, w} = sum_J (-1)^{n(J)}

summed over subsets ``J = {j_1 < ... < j_s}`` of a fixed lambda-chain such that

(a) ``u > u r_{j_1} > ... > u r_{j_1} ... r_{j_s} = w`` is a saturated
    decreasing chain in Bruhat order, and
(b) ``-mu = u r_{j_1} ... r_{j_s} (-lambda)``,

with ``n(J)`` the number of negative roots among ``beta_{j_1}, ..., beta_{j_s}``.

Weights are integer vectors in the ``epsilon`` basis; a root ``(a, b)`` denotes
``eps_a - eps_b`` and is positive iff ``a < b``.

<a id="schubmult.rings.schubert.chevalley.fundamental_weight"></a>

#### fundamental\_weight

```python
def fundamental_weight(k, n)
```

``omega_k = eps_1 + ... + eps_k`` as a length-``n`` vector.

<a id="schubmult.rings.schubert.chevalley.lambda_chain"></a>

#### lambda\_chain

```python
@cache
def lambda_chain(weight, n)
```

Reduced ``lambda``-chain for ``weight`` in ``A_{n-1}``, via Prop. 6.7.

Returns a tuple of ``(beta, root, k)``: ``root`` is the positive root
``alpha`` of the affine reflection ``r_j = s_{alpha, k}``, and ``beta`` is
the signed root ``b(r_j)`` whose sign determines ``(-1)^{n(J)}``.

<a id="schubmult.rings.schubert.chevalley.kt_chevalley_coefficients"></a>

#### kt\_chevalley\_coefficients

```python
def kt_chevalley_coefficients(u, weight, n=None)
```

Chevalley coefficients for ``e^weight . [O_u]``.

Returns ``{(w, mu): coefficient}`` with ``mu`` an integer weight vector.

