import sys
from functools import cache

from schubmult import ASx, Permutation, uncode
from schubmult.abc import x, y, z
from schubmult.rings import FA, FreeAlgebra, FreeAlgebraBasis, MonomialBasis, NilHeckeRing, PolynomialAlgebra, SchubertBasis, SingleSchubertRing, Sx, TensorRing
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, TensorModule
from schubmult.symbolic import S, expand_seq
from schubmult.utils.perm_utils import artin_sequences


def main():
    n = int(sys.argv[1])

    seqs = artin_sequences(n - 1)
    xring = PolynomialAlgebra(MonomialBasis(x, n - 1))

    yring = PolynomialAlgebra(MonomialBasis(y,n-1))
    zring = PolynomialAlgebra(MonomialBasis(z,n-1))

    unit_rc_module = RCGraphModule({RCGraph(): 1})
    result = TensorModule()

    # 100% positive!
    for seq in seqs:
        rc_graph_module_term = FA(*seq) * unit_rc_module
        # MODULE ACTION IS THE KEY
        rc_graph_module_term = RCGraphModule({rc: v for rc, v in rc_graph_module_term.items() if len(rc.perm) <= n})
        coprod = FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(), SchubertBasis, SchubertBasis)
        coprod = coprod.ring.from_dict({k: v for k, v in coprod.items() if len(k[0][0]) <= n and len(k[1][0]) <= n})

        zr_elem = zring(seq)
        xr_elem = xring(seq)
        module1 = TensorModule.ext_multiply(rc_graph_module_term, rc_graph_module_term)
        module2 = TensorModule.ext_multiply(rc_graph_module_term, rc_graph_module_term)
        result += TensorModule.ext_multiply(TensorModule.ext_multiply(module1, module2),coprod)

    result2 = result

    result = TensorModule()

    @cache
    def cached_fa(seq):
        return Sx([]).ring.from_dict({perm1[0]: coeff for perm1, coeff in FA(*seq).change_basis(SchubertBasis).items()})

    for key, value in result2.items():
        right_graph = key[0][1][1]
        left_graph = key[0][0][0]
        right_graph2 = key[0][1][0]
        coprod_perm_pair = key[1]
        left_graph2 = key[0][0][1]
        if right_graph == right_graph2 and left_graph == left_graph2:
            module1 = TensorModule.ext_multiply(value * Sx(left_graph.perm),cached_fa(left_graph.length_vector()))
            module2 = TensorModule.ext_multiply(Sx(right_graph.perm),cached_fa(right_graph.length_vector()))
            result += TensorModule.ext_multiply(TensorModule.ext_multiply(module1, module2),(ASx([]).ring@ASx([]).ring)(coprod_perm_pair))


    Permutation.print_as_code = True
    result_dict = {}
    failed = False
    for key, value in result.items():
        #right_graph = key[0][1][1]
        right_nil_perm = key[0][1][1]
        #$right_nil_perm = key[0][1][1].perm
        #left_graph = key[0][0][0]
        left_nil_perm = key[0][0][0]
        right_schub_perm = key[0][1][0]
        coprod_perm_pair = key[1]
        left_schub_perm = key[0][0][1]
        if right_nil_perm == right_schub_perm and left_nil_perm == right_nil_perm:
            result_dict[coprod_perm_pair] = result_dict.get(coprod_perm_pair, 0) + value * Sx(left_schub_perm)

    for coprod_key in result_dict:
        ((left_coprod_perm, _), (right_coprod_perm, _)) = coprod_key
        if any(len(permperm) > n for permperm, val in (Sx(left_coprod_perm) * Sx(right_coprod_perm)).items() if val != S.Zero):
            continue
        testval = result_dict[coprod_key] - Sx(left_coprod_perm) * Sx(right_coprod_perm)
        if testval.expand() != S.Zero:
            print(f"Failures for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: {testval=} {result_dict[coprod_key]=}")
            failed = True
        else:
            print(f"Success for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: Sx({left_coprod_perm.trimcode})*Sx({right_coprod_perm.trimcode})={result_dict[coprod_key]}")

    if not failed:
        print("YEAH!!!")
        print("YEAH!!!", file=sys.stderr)


if __name__ == "__main__":
    main()
