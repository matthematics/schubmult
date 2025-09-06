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
    dual_schub_tensor_square = TensorRing(FreeAlgebra(SchubertBasis) @ FreeAlgebra(SchubertBasis))
    main_ring = (TensorRing(NH, zring)@TensorRing(Sx([]).ring, NH)) @ dual_schub_tensor_square

    unit_rc_module = RCGraphModule({RCGraph(): 1})
    result = main_ring.zero

    # 100% positive!
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

        main_ring1 = main_ring.rings[0].rings[0]
        main_ring2 = main_ring.rings[0].rings[1]
        first_term = main_ring1.zero
        second_term = main_ring2.zero
        zr_elem = zring(expand_seq(seq, z))
        xr_elem = Sx(expand_seq(seq, x))
        for rc_graph, coeff in rc_graph_module_term.items():
            first_term += coeff * main_ring1.ext_multiply(NH(rc_graph.perm), zr_elem)
            second_term += coeff * main_ring2.ext_multiply(xr_elem, NH(rc_graph.perm))

        accum_term += main_ring.ext_multiply(main_ring.rings[0].ext_multiply(first_term,second_term),coprod)
        result += accum_term

    Permutation.print_as_code = True
    result_dict = {}
    failed = False
    for key, value in result.items():
        right_nil_perm = key[0][1][1]
        left_nil_perm = key[0][0][0]
        right_schub_perm = key[0][1][0]
        coprod_perm_pair = key[1]
        left_schub_perm = key[0][0][1]
        #if right_nil_perm == right_schub_perm and left_nil_perm == right_nil_perm:
        if right_nil_perm == right_schub_perm and left_nil_perm == right_schub_perm:
            # print(f"{left_nil_perm.trimcode=},{right_schub_perm.trimcode=},{left_schub_perm.trimcode=},{right_nil_perm.trimcode=},{coprod_perm_pair=} contributes {value}?")
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
