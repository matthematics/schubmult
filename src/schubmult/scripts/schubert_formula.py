import sys

from schubmult import Permutation
from schubmult.abc import x, y, z
from schubmult.rings import FA, FreeAlgebra, FreeAlgebraBasis, NilHeckeRing, SchubertBasis, SingleSchubertRing, Sx, TensorRing
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule
from schubmult.symbolic import S, expand_seq
from schubmult.utils.perm_utils import artin_sequences


def main():
    n = int(sys.argv[1])

    seqs = artin_sequences(n - 1)

    yring = SingleSchubertRing(y)
    zring = SingleSchubertRing(z)
    dual_schub_tensor_square_nil_ring = TensorRing(FreeAlgebra(SchubertBasis) @ FreeAlgebra(SchubertBasis), NilHeckeRing(x))
    schub_dualschub_nil_ring = TensorRing(Sx([]).ring, dual_schub_tensor_square_nil_ring)
    main_ring = TensorRing(zring, schub_dualschub_nil_ring)

    unit_rc_module = RCGraphModule({RCGraph(): 1})
    result = main_ring.zero

    for seq in seqs:
        rc_graph_module_term = FA(*seq) * unit_rc_module
        accum_term = main_ring.zero
        rc_graph_module_term = RCGraphModule({rc: v for rc, v in rc_graph_module_term.items() if len(rc.perm) <= n})
        coprod = FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(), SchubertBasis, SchubertBasis)
        coprod = coprod.ring.from_dict({k: v for k, v in coprod.items() if len(k[0][0]) <= n and len(k[1][0]) <= n})
        for rc_graph, coeff in rc_graph_module_term.items():
            accum_term += coeff * main_ring.ext_multiply(
                zring(expand_seq(seq, z)),
                schub_dualschub_nil_ring.ext_multiply(Sx(expand_seq(seq, x)), dual_schub_tensor_square_nil_ring.ext_multiply(coprod, dual_schub_tensor_square_nil_ring.rings[1](rc_graph.perm))),
            )
        result += accum_term

    Permutation.print_as_code = True
    result_dict = {}
    failed = False
    for key, value in result.items():
        if key[1][1][1] == key[1][0]:
            result_dict[key[1][1][0]] = result_dict.get(key[1][1][0], yring.zero) + value * yring(key[0])

    for perm in result_dict:
        if any(len(permperm) > n for permperm, val in (Sx(perm[0][0]) * Sx(perm[1][0])).items() if val != S.Zero):
            continue
        testval = result_dict[perm] - yring(perm[0][0]) * yring(perm[1][0])
        if testval.expand() != S.Zero:
            print(f"Failures for {(perm[0][0].trimcode, perm[1][0].trimcode)}: {[(perm2, key) for perm2, key in testval.items() if key != 0]}")
            failed = True
        else:
            print(f"Success for {(perm[0][0].trimcode, perm[1][0].trimcode)}: Sx({perm[0][0].trimcode})*Sx({perm[1][0].trimcode})={result_dict[perm]}")

    if not failed:
        print("YEAH!!!")
        print("YEAH!!!", file=sys.stderr)


if __name__ == "__main__":
    main()
