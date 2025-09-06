import sys

from schubmult import ASx, Permutation
from schubmult.abc import x, y, z
from schubmult.rings import FA, FreeAlgebra, FreeAlgebraBasis, NilHeckeRing, SchubertBasis, SingleSchubertRing, Sx, TensorRing
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule
from schubmult.symbolic import S, expand_seq
from schubmult.utils.perm_utils import artin_sequences


def main():
    n = int(sys.argv[1])

    seqs = artin_sequences(n - 1)
    NH = NilHeckeRing(x)
    yring = SingleSchubertRing(y)
    zring = SingleSchubertRing(z)
    dual_schub_tensor_square_nil_ring = TensorRing(FreeAlgebra(SchubertBasis) @ FreeAlgebra(SchubertBasis), NH)
    schub_dualschub_nil_ring = TensorRing(Sx([]).ring, dual_schub_tensor_square_nil_ring)
    main_ring1 = TensorRing(zring, schub_dualschub_nil_ring)
    main_ring = TensorRing(NH, main_ring1)

    unit_rc_module = RCGraphModule({RCGraph(): 1})
    result = main_ring.zero

    for seq in seqs:
        rc_graph_module_term = FA(*seq) * unit_rc_module
        # test_elem = ASx([]).ring.from_dict({(rc.perm, len(seq)): coeff for rc, coeff in rc_graph_module_term.items()})
        # schubert_elem = FA(*seq).change_basis(SchubertBasis)
        # assert  test_elem == schubert_elem, f"Basis change failed for {seq}: {test_elem} != {schubert_elem}"
        # MODULE ACTION IS THE KEY
        accum_term = main_ring.zero
        rc_graph_module_term = RCGraphModule({rc: v for rc, v in rc_graph_module_term.items() if len(rc.perm) <= n})
        coprod = FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(), SchubertBasis, SchubertBasis)
        coprod = coprod.ring.from_dict({k: v for k, v in coprod.items() if len(k[0][0]) <= n and len(k[1][0]) <= n})
        accum_term = main_ring.zero
        for rc_graph1, coeff1 in rc_graph_module_term.items():
            for rc_graph2, coeff2 in rc_graph_module_term.items():
                accum_term += coeff1 * coeff2 * main_ring.ext_multiply(NH(rc_graph1.perm),main_ring1.ext_multiply(
                    zring(expand_seq(seq, z)),
                    schub_dualschub_nil_ring.ext_multiply(Sx(expand_seq(seq, x)), dual_schub_tensor_square_nil_ring.ext_multiply(coprod, dual_schub_tensor_square_nil_ring.rings[1](rc_graph2.perm))),
                ))

        result += accum_term

    Permutation.print_as_code = True
    result_dict = {}
    failed = False
    for key, value in result.items():
        right_nil_perm = key[1][1][1][1]
        left_nil_perm = key[0]
        middle_schub_perm = key[1][1][0]
        coprod_perm_pair = key[1][1][1][0]
        left_schub_perm = key[1][0]
        if right_nil_perm == middle_schub_perm and left_nil_perm == right_nil_perm:
            result_dict[coprod_perm_pair] = result_dict.get(coprod_perm_pair, yring.zero) + value * yring(left_schub_perm)

    for coprod_key in result_dict:
        ((left_coprod_perm, _), (right_coprod_perm, _)) = coprod_key
        if any(len(permperm) > n for permperm, val in (Sx(left_coprod_perm) * Sx(right_coprod_perm)).items() if val != S.Zero):
            continue
        testval = result_dict[coprod_key] - yring(left_coprod_perm) * yring(right_coprod_perm)
        if testval.expand() != S.Zero:
            print(f"Failures for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: {[(perm2, key) for perm2, key in testval.items() if key != 0]}")
            failed = True
        else:
            print(f"Success for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: Sx({left_coprod_perm.trimcode})*Sx({right_coprod_perm.trimcode})={result_dict[coprod_key]}")

    if not failed:
        print("YEAH!!!")
        print("YEAH!!!", file=sys.stderr)


if __name__ == "__main__":
    main()
