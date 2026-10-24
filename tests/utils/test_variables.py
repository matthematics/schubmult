import pytest

from schubmult import ChangeOfVariables, GeneratingSet
from schubmult.symbolic import S, expand, symbols, sympify_sympy
from schubmult.symbolic.poly.variables import DEF_GENSET_SIZE, CustomGeneratingSet, MaskedGeneratingSet, NotEnoughGeneratorsError


def test_generating_set_grows_geometrically_and_preserves_symbols():
    gens = GeneratingSet("growth")
    assert len(gens) == DEF_GENSET_SIZE == 10
    original = tuple(gens)
    original_hash = hash(gens)
    assert gens[10] == symbols("growth_10")
    assert len(gens) == 20
    assert gens[20] == symbols("growth_20")
    assert len(gens) == 40
    assert gens[1000] == symbols("growth_1000")
    assert len(gens) >= 1001
    assert all(gens[i] is symbol for i, symbol in enumerate(original))
    assert hash(gens) == original_hash
    assert GeneratingSet("growth") is gens
    assert gens(1001) == gens[1000]


def test_generating_set_slices_expand_only_for_explicit_bounds():
    gens = GeneratingSet("slice_growth")
    assert len(tuple(gens)) == len(gens[:]) == 10
    assert gens[-1] == symbols("slice_growth_9")
    with pytest.raises(IndexError):
        gens[-11]
    assert gens[8:15:2] == tuple(symbols(f"slice_growth_{i}") for i in (8, 10, 12, 14))
    assert len(gens) == 20
    assert gens[110:105:-2] == tuple(symbols(f"slice_growth_{i}") for i in (110, 108, 106))
    size = len(gens)
    assert tuple(gens) == gens[:]
    assert len(gens) == size
    with pytest.raises(ValueError):
        gens[::0]
    with pytest.raises(TypeError):
        gens[1.5]


def test_generating_set_recognizes_unmaterialized_symbols():
    gens = GeneratingSet("name_with_underscores")
    assert gens.index(symbols("name_with_underscores_105")) == 105
    assert gens.index("name_with_underscores_205") == 205
    assert gens.index(sympify_sympy("name_with_underscores_305")) == 305
    assert symbols("name_with_underscores_405") in gens
    assert gens.index(symbols("name_with_underscores_01")) == -1
    assert gens.index(gens[1] + 1) == -1
    assert gens.index([]) == -1
    assert gens.index(None) == -1
    size = len(gens)
    assert gens.index(symbols("other_10000")) == -1
    assert len(gens) == size


def test_masked_generating_set_remains_consistent_after_growth():
    gens = GeneratingSet("masked_growth")
    masked = MaskedGeneratingSet(gens, [2, 4, 105])
    assert len(masked) == 8
    assert masked[2] == gens[3]
    assert masked[3] == gens[5]
    assert masked[103] == gens[106]
    assert masked.index(gens[106]) == 103
    assert masked._index_lookup[gens[106]] == 103
    assert masked.index(gens[105]) == -1
    assert masked.index(symbols("masked_growth_205")) == 202
    assert masked[150:155] == [masked[i] for i in range(150, 155)]
    assert masked[-1] == gens[-1]
    nested = MaskedGeneratingSet(masked, [1])
    assert nested[200] == masked[201]
    assert nested.index(nested[200]) == 200
    selected = masked.complement()
    assert list(selected) == [gens[0], gens[2], gens[4], gens[105]]
    assert not selected.is_infinite
    with pytest.raises(IndexError):
        selected[4]


def test_custom_generating_sets_remain_finite():
    gens = CustomGeneratingSet(symbols("finite_0 finite_1 finite_2"))
    masked = MaskedGeneratingSet(gens, [1])
    assert list(MaskedGeneratingSet(CustomGeneratingSet(()), ())) == []
    assert not gens.is_infinite
    assert list(masked) == [gens[0], gens[2]]
    assert masked[1:5] == [gens[2]]
    assert list(masked.complement()) == [gens[0], gens[1]]
    with pytest.raises(IndexError):
        gens[3]
    with pytest.raises(IndexError):
        masked[2]


def test_change_of_variables_indexes_and_substitutes_generators():
    source = GeneratingSet("source")
    target = GeneratingSet("target")
    change = ChangeOfVariables(source, target, lambda index, new_gset: 1 - 1 / new_gset[index])

    assert change.label == target.label
    assert change[0] == 1 - 1 / target[0]
    assert change(2) == change[1]
    assert change.index(source[1]) == 1

    expr = source[0] + source[1] ** 2
    expected = (1 - 1 / target[0]) + (1 - 1 / target[1]) ** 2
    assert expand(change.subs(expr) - expected) == S.Zero


def test_change_of_variables_leaves_other_symbols_untouched():
    source = GeneratingSet("source")
    target = GeneratingSet("target")
    other = GeneratingSet("other")
    change = ChangeOfVariables(source, target, lambda index, new_gset: new_gset[index] ** 2)

    assert change.subs(source[0] + other[0]) == target[0] ** 2 + other[0]


def test_change_of_variables_supports_unbounded_generating_sets():
    source = GeneratingSet("cov_source")
    target = GeneratingSet("cov_target")
    change = ChangeOfVariables(source, target, lambda i, gens: gens[i] ** 2)
    assert change[105] == target[105] ** 2
    assert change[150:153] == tuple(target[i] ** 2 for i in range(150, 153))
    assert change.subs(symbols("cov_source_205") + 1) == target[205] ** 2 + 1


@pytest.mark.parametrize("kind", ["e", "E", "h", "H", "E_q"])
def test_symmetric_polynomials_materialize_required_generators(kind):
    from schubmult.symbolic.symmetric_polynomials.complete_sym import H, h
    from schubmult.symbolic.symmetric_polynomials.elem_sym import E, e
    from schubmult.symbolic.symmetric_polynomials.qelem_sym import E_q

    x = GeneratingSet(f"symmetric_x_{kind}")
    y = GeneratingSet(f"symmetric_y_{kind}")
    cls = {"e": e, "E": E, "h": h, "H": H, "E_q": E_q}[kind]
    polynomial = cls(1, 105, x) if kind in ("e", "h") else cls(1, 105, x, y)
    assert polynomial.numvars == 105
    assert polynomial.genvars == tuple(sorted(x[:105], key=lambda v: sympify_sympy(v).sort_key()))
    assert x[104] in polynomial.free_symbols
    with pytest.raises(NotEnoughGeneratorsError):
        cls(1, 105, CustomGeneratingSet(x[:2])) if kind in ("e", "h") else cls(1, 105, CustomGeneratingSet(x[:2]), y)


def test_schubert_ring_accepts_descents_beyond_cached_prefix():
    from schubmult import Permutation, SingleSchubertRing

    gens = GeneratingSet("large_ring")
    ring = SingleSchubertRing(gens)
    permutation = Permutation([*range(1, 105), 106, 105])
    element = ring(permutation)
    assert ring(list(permutation)) == element
    shifted_ring = SingleSchubertRing(MaskedGeneratingSet(gens, range(1, 105)))
    assert shifted_ring([2, 1]).expand() == gens[105]


@pytest.mark.parametrize("down", [False, True])
def test_double_schubert_ring_distinguishes_cache_size_from_finite_limits(down):
    from schubmult import Permutation
    from schubmult.rings.schubert.double_schubert_ring import DoubleSchubertRing, DoubleSchubertRingDown

    gens = GeneratingSet(f"double_large_ring_{down}")
    coefficients = GeneratingSet("double_large_ring_coeffs")
    cls = DoubleSchubertRingDown if down else DoubleSchubertRing
    ring = cls(gens, coefficients)
    permutation = Permutation([*range(1, 105), 106, 105])
    assert permutation in ring(permutation)
    assert permutation in ring(list(permutation))
    assert len(gens) == 10
    finite_ring = cls(CustomGeneratingSet(gens[:3]), coefficients)
    with pytest.raises(NotEnoughGeneratorsError):
        finite_ring(permutation)
