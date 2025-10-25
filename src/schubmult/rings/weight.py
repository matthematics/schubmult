from typing import Sequence, Tuple


class Weight(tuple):
    """Immutable integral weight: tuple of ints."""
    @classmethod
    def from_seq(cls, seq: Sequence[int]) -> "Weight":
        return Weight(tuple(int(x) for x in seq))

    def as_tuple(self) -> Tuple[int, ...]:
        return tuple(self)


def permute_weight(perm, weight: Sequence[int], dot_action: bool = False) -> Weight:
    """
    Return the weight obtained by the Weyl-group action of `perm` on `weight`.

    Uses the inverse permutation (~perm) on the fly to read preimages without
    allocating an intermediate images list.

    Convention (type A character / variable permutation):
      - perm[i] returns the image of i+1 as a 1-based int.
      - (w·λ)_i = λ_{w^{-1}(i)} for the plain action.
      - If dot_action is True: w·λ = w(λ+ρ) − ρ with ρ = (n-1,...,0).
    """
    lam = tuple(int(x) for x in weight)
    n = len(lam)

    # use inverse permutation object; (~perm)[i] should give preimage of i+1 as 1-based int
    inv_perm = ~perm

    if dot_action:
        rho = tuple(n - 1 - i for i in range(n))
        lam_rho = tuple(lam[i] + rho[i] for i in range(n))
        permuted = []
        for i in range(n):
            try:
                pre = int(inv_perm[i])  # pre is 1-based index of preimage of i+1
                if 1 <= pre <= n:
                    permuted.append(lam_rho[pre - 1])
                else:
                    permuted.append(lam_rho[i])
            except Exception:
                permuted.append(lam_rho[i])
        result = tuple(permuted[i] - rho[i] for i in range(n))
        return Weight.from_seq(result)

    # plain action: use inverse permutation to pick coordinates
    result = []
    for i in range(n):
        try:
            pre = int(inv_perm[i])
            if 1 <= pre <= n:
                result.append(lam[pre - 1])
            else:
                result.append(lam[i])
        except Exception:
            result.append(lam[i])
    return Weight.from_seq(result)


def fixes_weight(perm, weight: Sequence[int], dot_action: bool = False) -> bool:
    """Return True iff perm fixes the integral weight under chosen action."""
    return permute_weight(perm, weight, dot_action) == Weight.from_seq(weight)

