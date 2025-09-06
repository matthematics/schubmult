import sys

from schubmult import Permutation
from schubmult.abc import x, y, z
from schubmult.rings import FA, FreeAlgebra, FreeAlgebraBasis, MonomialBasis, NilHeckeRing, PolynomialAlgebra, SchubertBasis, SingleSchubertRing, Sx, TensorRing, WordBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule
from schubmult.symbolic import S, expand_seq
from schubmult.utils.perm_utils import artin_sequences


def main():
    n = int(sys.argv[1])

    # These functions must be defined/imported elsewhere in your codebase
    seqs = artin_sequences(n - 1)

    yring = SingleSchubertRing(y)
    zring = SingleSchubertRing(z)
    ring = TensorRing(FreeAlgebra(SchubertBasis) @ FreeAlgebra(SchubertBasis), NilHeckeRing(x))
    ring2 = TensorRing(Sx([]).ring, ring)
    ring3333 = TensorRing(ring.rings[1], ring2)
    ringbob = TensorRing(zring, ring3333)

    
    mod1 = RCGraphModule({RCGraph(): 1})
    result = ringbob.zero

    for seq in seqs:
        modmod = (FA(*seq) * mod1).as_nil_hecke(x)
        # print(f"{seq=}")
        # print(f"{modmod=}")
        ding = ringbob.zero
        modmod = modmod.ring.from_dict({k: v for k, v in modmod.items() if len(k) <= n})
        coprod = FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(), SchubertBasis, SchubertBasis)
        coprod = coprod.ring.from_dict({k: v for k, v in coprod.items() if len(k[0][0]) <= n and len(k[1][0]) <= n})
        for nil_key in modmod.keys():
            ding += ringbob.ext_multiply(
                zring(expand_seq(seq, z)),
                ring3333.ext_multiply(ring.rings[1](nil_key), ring2.ext_multiply(Sx(expand_seq(seq, x)), ring.ext_multiply(coprod, ring.rings[1](nil_key)))),
            )
        result += ding

    # print(result)
    Permutation.print_as_code = True
    separate = {}
    failed = False
    for key, value in result.items():
        if key[1][1][1][1] == key[1][1][0]:
            # print(f"{key[1][1][0]=} {key[0]=}", file=sys.stderr)
            separate[key[1][1][1][0]] = separate.get(key[1][1][1][0], yring.zero) + value * yring(key[0])

    for perm in separate:
        if any(len(permperm) > n for permperm, val in (Sx(perm[0][0]) * Sx(perm[1][0])).items() if val != S.Zero):
            continue
        # print(f"Test {(perm[0][0].trimcode, perm[1][0].trimcode)}: {separate[perm]}")
        testval = separate[perm] - yring(perm[0][0]) * yring(perm[1][0])
        if testval.expand() != S.Zero:
            print(f"Failures for {(perm[0][0].trimcode, perm[1][0].trimcode)}: {[(perm2, key) for perm2, key in testval.items() if key != 0]}")
            failed = True
        else:
            print(f"Success for {(perm[0][0].trimcode, perm[1][0].trimcode)}: Sx({perm[0][0].trimcode})*Sx({perm[1][0].trimcode})={separate[perm]}")

    if not failed:
        print("YEAH!!!")
        print("YEAH!!!", file=sys.stderr)


if __name__ == "__main__":
    main()
