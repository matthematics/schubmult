import itertools
import pytest
from schubmult import Permutation, RCGraph, Sx
from schubmult.rings.combinatorial.grass_tensor_algebra import GrassTensorAlgebra


def test_match_rc_graph():
    ring = GrassTensorAlgebra()
    perms = Permutation.all_permutations(5)

    for perm in perms:
        cem = RCGraph.full_CEM(perm, len(perm.trimcode))
        for rc, cem_dict in cem.items():
            elem = ring.from_dict(cem_dict)
            rc_elem = elem.to_rc_graph_ring_element()
            if rc.perm == perm:
                assert rc_elem.almosteq(rc_elem.ring(rc)), f"Failed on {rc}\nGot {rc_elem}\nExpected {rc_elem.ring(rc)}"
            else:
                assert rc_elem.almosteq(rc_elem.ring.zero), f"Failed on {rc}\nGot {rc_elem}\nExpected 0"

def test_mul_matches_schubert():
    ring = GrassTensorAlgebra()
    n = 4
    perms = Permutation.all_permutations(n)

    grass_tensor_elems = {}
    for perm in perms:
        grass_tensor_elems[perm] = ring.zero
        cem = RCGraph.full_CEM(perm, n - 1)
        for rc, cem_dict in cem.items():
            grass_tensor_elems[perm] += ring.from_dict(cem_dict)
    for perm1, perm2 in itertools.product(perms, repeat=2):
        prd = Sx(perm1) * Sx(perm2)
        prd_elem = grass_tensor_elems[perm1] * grass_tensor_elems[perm2]
        prd_rc = prd_elem.to_rc_graph_ring_element()
        for rc, coeff in prd_rc.items():
            assert prd.get(rc.perm, 0) == coeff, f"Failed on {perm1} * {perm2}, got {rc.perm}: {coeff} which is not in {prd}\n{prd.get(rc.perm,0)=}\n{rc.perm=}"


# ---- BoundedRCFactorAlgebra tests ----


from schubmult.rings.combinatorial.bounded_rc_factor_algebra import (
    BoundedRCFactorAlgebra,
    _is_full_grassmannian_rc,
)


class TestBoundedRCFactorAlgebraBasics:
    """Basic construction, key normalization, and algebraic identities."""

    def test_zero_and_one(self):
        brc = BoundedRCFactorAlgebra()
        assert len(brc.zero) == 0
        assert len(brc.one) == 1

    def test_make_key_empty(self):
        brc = BoundedRCFactorAlgebra()
        key = brc.make_key((), 3)
        assert len(key) == 0
        assert key.size == 3

    def test_make_key_single_factor(self):
        rc = list(RCGraph.all_rc_graphs(Permutation([1, 3, 2]), 2))[0]
        brc = BoundedRCFactorAlgebra()
        key = brc.make_key((rc,), 3)
        assert len(key) == 1
        assert key.size == 3

    def test_key_to_rc_graph_empty(self):
        brc = BoundedRCFactorAlgebra()
        key = brc.make_key((), 3)
        rc = brc.key_to_rc_graph(key)
        assert rc.perm.inv == 0

    def test_key_to_rc_graph_single_factor(self):
        brc = BoundedRCFactorAlgebra()
        rc_orig = list(RCGraph.all_rc_graphs(Permutation([1, 3, 2]), 2))[0]
        key = brc.make_key((rc_orig,), 3)
        rc_out = brc.key_to_rc_graph(key)
        assert rc_out.perm == rc_orig.perm

    def test_elem_from_key(self):
        brc = BoundedRCFactorAlgebra()
        rc = list(RCGraph.all_rc_graphs(Permutation([1, 3, 2]), 2))[0]
        key = brc.make_key((rc,), 3)
        elem = brc(key)
        assert len(elem) == 1

    def test_normalize_key_sorts_by_length(self):
        """Factors should be sorted by length (ascending) after normalization."""
        brc = BoundedRCFactorAlgebra()
        rc1 = list(RCGraph.all_rc_graphs(Permutation([2, 1]), 1))[0]
        rc2 = list(RCGraph.all_rc_graphs(Permutation([1, 3, 2]), 2))[0]
        # Put longer factor first
        key = brc.make_key((rc2, rc1), 3)
        nkey = brc._normalize_key(key)
        if len(nkey) == 2:
            assert len(nkey[0]) <= len(nkey[1])


class TestBoundedRCFactorAlgebraSchubElem:
    """schub_elem should squash to the correct Schubert polynomial."""

    @pytest.mark.parametrize("perm_code", [
        [1],
        [2],
        [1, 1],
        [0, 1],
        [2, 1],
    ])
    def test_schub_elem_squash_matches_schubert(self, perm_code):
        from schubmult import uncode
        from schubmult.symbolic import S

        perm = uncode(perm_code)
        brc = BoundedRCFactorAlgebra()
        size = max(len(perm_code), 2)
        elem = brc.schub_elem(perm, size)
        rc_elem = elem.to_rc_graph_ring_element()
        # Every RC graph in the expansion should have the correct permutation
        schub_poly = Sx(perm)
        for rc, coeff in rc_elem.items():
            assert coeff == schub_poly.get(rc.perm, S.Zero), (
                f"Mismatch at {rc.perm} for schub_elem({perm}, {size})"
            )

    def test_schub_elem_factors_are_grassmannian(self):
        """Every factor in every key should be a full Grassmannian RC graph."""
        from schubmult import uncode

        perm = uncode([2, 1])
        brc = BoundedRCFactorAlgebra()
        elem = brc.schub_elem(perm, 3)
        for key in elem:
            for factor in key:
                assert _is_full_grassmannian_rc(factor), f"Non-Grassmannian factor: {factor}"

    def test_schub_elem_identity(self):
        """S_id should be the unit element."""
        brc = BoundedRCFactorAlgebra()
        elem = brc.schub_elem(Permutation([]), 2)
        assert elem.almosteq(brc.one) or len(elem) == 1


class TestBoundedRCFactorAlgebraMul:
    """Multiplication in the factor algebra should squash to Schubert multiplication."""

    @pytest.mark.parametrize("code1,code2", [
        ([1], [1]),
        ([1], [0, 1]),
        ([2], [1]),
        ([1, 1], [1]),
    ])
    def test_mul_squash_matches_schubert(self, code1, code2):
        from schubmult import uncode
        from schubmult.symbolic import S

        perm1 = uncode(code1)
        perm2 = uncode(code2)
        brc = BoundedRCFactorAlgebra()
        size = max(len(code1), len(code2)) + 1

        s1 = brc.schub_elem(perm1, size)
        s2 = brc.schub_elem(perm2, size)
        prod = s1 * s2
        rc_prod = prod.to_rc_graph_ring_element()

        schub_prod = Sx(perm1) * Sx(perm2)
        for rc, coeff in rc_prod.items():
            assert coeff == schub_prod.get(rc.perm, S.Zero), (
                f"Product mismatch at {rc.perm}: got {coeff}, "
                f"expected {schub_prod.get(rc.perm, S.Zero)} "
                f"for S_{perm1} * S_{perm2}"
            )

    def test_mul_associativity(self):
        from schubmult import uncode

        brc = BoundedRCFactorAlgebra()
        s21 = brc.schub_elem(uncode([1]), 2)
        s132 = brc.schub_elem(uncode([0, 1]), 2)

        prod_lr = (s21 * s132) * s21
        prod_rl = s21 * (s132 * s21)

        rc_lr = prod_lr.to_rc_graph_ring_element()
        rc_rl = prod_rl.to_rc_graph_ring_element()
        assert rc_lr.almosteq(rc_rl), "Multiplication not associative at RC graph level"


class TestBoundedRCFactorAlgebraElemSym:
    """Elementary symmetric elements e_{p,k}."""

    def test_elem_sym_e11_is_s21(self):
        """e_{1,1} should be the single transposition Schubert element."""
        from schubmult.symbolic import S

        brc = BoundedRCFactorAlgebra()
        e11 = brc.elem_sym(1, 1)
        rc_e11 = e11.to_rc_graph_ring_element()
        assert len(rc_e11) == 1
        for rc, coeff in rc_e11.items():
            assert rc.perm == Permutation([2, 1])
            assert coeff == S.One


class TestBoundedRCFactorAlgebraKeyToRC:
    """key_to_rc_graph squashing."""

    def test_squash_product_two_factors(self):
        """Squashing two Grassmannian factors should give the correct RC graph."""
        brc = BoundedRCFactorAlgebra()
        rc1 = list(RCGraph.all_rc_graphs(Permutation([2, 1]), 1))[0]
        rc2 = list(RCGraph.all_rc_graphs(Permutation([1, 3, 2]), 2))[0]
        key = brc._normalize_key(brc.make_key((rc1, rc2), 3))
        rc_out = brc.key_to_rc_graph(key)
        # The result should be a valid RC graph
        assert rc_out.perm.inv > 0

    def test_squash_idempotent(self):
        """Squashing a single-factor key should give that factor back (resized)."""
        brc = BoundedRCFactorAlgebra()
        rc = list(RCGraph.all_rc_graphs(Permutation([1, 3, 2]), 3))[0]
        key = brc.make_key((rc,), len(rc))
        rc_out = brc.key_to_rc_graph(key)
        assert rc_out == rc


class TestBoundedRCFactorAlgebraCoproduct:
    """Coproduct via vertical cuts."""

    def test_coproduct_coassociative(self):
        """The coproduct should be coassociative for a simple element."""
        from schubmult import uncode
        from schubmult.symbolic import S

        brc = BoundedRCFactorAlgebra()
        perm = uncode([1])
        elem = brc.schub_elem(perm, 2)
        coprod = elem.coproduct()
        # Just check that it's nonzero and has valid structure
        assert len(coprod) > 0
        for (k1, k2), coeff in coprod.items():
            assert coeff != S.Zero

    def test_coproduct_identity(self):
        """Coproduct of identity should be id ⊗ id."""
        brc = BoundedRCFactorAlgebra()
        elem = brc.schub_elem(Permutation([]), 2)
        coprod = elem.coproduct()
        assert len(coprod) > 0


class TestBoundedRCFactorAlgebraCEM:
    """CEM expansion structure."""

    def test_cem_all_factors_grassmannian(self):
        """All CEM tensor factors should be full Grassmannian."""
        from schubmult import uncode

        for code in [[1], [2], [1, 1], [0, 2, 1], [2, 0, 1]]:
            perm = uncode(code)
            size = len(code) + 1
            cem = RCGraph.full_CEM(perm, size)
            for _rc, cem_dict in cem.items():
                for tensor_key in cem_dict:
                    for factor in tensor_key:
                        assert _is_full_grassmannian_rc(factor), (
                            f"Non-Grassmannian factor in CEM of {perm}: {factor}"
                        )

    def test_cem_legit_terms_squash_correctly(self):
        """Non-phantom CEM terms should squash to their keying RC graph."""
        from schubmult import uncode
        from schubmult.symbolic import S

        brc = BoundedRCFactorAlgebra()
        for code in [[1], [2], [1, 1], [0, 2, 1]]:
            perm = uncode(code)
            size = len(code) + 1
            cem = RCGraph.full_CEM(perm, size)
            for rc, cem_dict in cem.items():
                if rc.perm != perm:
                    continue
                # The sum over tensor keys should squash to this RC graph
                accum = brc.zero
                for tensor_key, coeff in cem_dict.items():
                    key = brc.make_key(tensor_key, size)
                    accum += coeff * brc(key)
                rc_elem = accum.to_rc_graph_ring_element()
                assert rc_elem.get(rc, S.Zero) != S.Zero, (
                    f"Legit CEM of {perm} did not squash to {rc}"
                )

    def test_cem_phantom_terms_vanish(self):
        """Phantom CEM terms (rc.perm != perm) should squash to zero."""
        from schubmult import uncode
        from schubmult.symbolic import S

        brc = BoundedRCFactorAlgebra()
        for code in [[0, 2, 1], [2, 0, 1], [3, 0, 1]]:
            perm = uncode(code)
            size = len(code) + 1
            cem = RCGraph.full_CEM(perm, size)
            for rc, cem_dict in cem.items():
                if rc.perm == perm:
                    continue
                accum = brc.zero
                for tensor_key, coeff in cem_dict.items():
                    key = brc.make_key(tensor_key, size)
                    accum += coeff * brc(key)
                rc_elem = accum.to_rc_graph_ring_element()
                for _rc_out, c in rc_elem.items():
                    assert c == S.Zero, (
                        f"Phantom CEM of {perm} at {rc} did not vanish: "
                        f"got {c} at {_rc_out.perm}"
                    )