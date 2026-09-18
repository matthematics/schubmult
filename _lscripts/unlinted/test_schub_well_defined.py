from schubmult import *

if __name__ == "__main__":
    import sys
    from sympy import pretty_print

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    r = BoundedRCFactorAlgebra()
    doms = [perm for perm in perms if perm.is_dominant]
    for perm in perms:
        schub_elem = None
        for dom in doms:
            if len(dom.trimcode) != len(perm.trimcode):
                continue
            if (perm * (~dom)).inv == dom.inv - perm.inv:
                partition = tuple((~dom).trimcode)
                try:
                    test_schub_elem = r.schub_elem(perm, len(perm.trimcode), partition=partition)
                except ValueError as e:
                    #print(f"ERROR: perm={perm}, dom={dom}, partition={partition}, error={e}")
                    continue
                if schub_elem is None:
                    schub_elem = test_schub_elem
                else:
                    if not schub_elem.almosteq(test_schub_elem):
                        print(f"FAIL: perm={perm}, dom={dom}, partition={partition}")
                        print(f"  schub_elem: {schub_elem}")
                        print(f"  test_schub_elem: {test_schub_elem}")
                        sys.exit(1)
                    print(f"PASS: perm={perm} {partition=}")