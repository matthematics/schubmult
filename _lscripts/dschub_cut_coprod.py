from schubmult import *
from schubmult.rings.free_algebra import *

r = RCGraphRing()
def _canonical_rc(rc):
    if rc.is_quasi_yamanouchi:
        return rc
    #row = 1
    for i in range(1, len(rc)):
        if max(rc[i], default=0) < min(rc[i - 1], default=0) and min(rc[i - 1], default=0) >= i + 1:
            new_rc = [*rc]
            new_rc[i] = tuple([*rc[i - 1], *rc[i]])
            new_rc[i - 1] = ()
            return _canonical_rc(RCGraph(new_rc))
    raise ValueError(f"Should never reach here {rc=} {rc.is_quasi_yamanouchi=}")
gcache = {}
acache = {}

def fast_coproduct(w, length, rc=None):
    global gcache
    if (w, length) in gcache:
        return gcache[(w, length)]
    if rc is None:
        rc = RCGraph.principal_rc(w, length)
    #length = len(rc)

    if length <= 1:
        return ASx(w, length).coproduct()
    mid = length // 2
    bottom_rc, top_rc = rc.vertical_cut(mid)
    bperm, tperm = bottom_rc.perm, top_rc.perm
    left_coprod = fast_coproduct(bperm, len(bottom_rc), bottom_rc)
    right_coprod = fast_coproduct(tperm, len(top_rc), top_rc)
    full_coprod = left_coprod * right_coprod
    mixed_elem = ASx(bperm, len(bottom_rc)) * ASx(tperm, len(top_rc))
    for (perm, _), coeff in mixed_elem.items():
        if perm != w:
            full_coprod -= coeff * fast_coproduct(perm, length)
    gcache[(w, length)] = full_coprod
    return full_coprod

def fast_schub(w, length, rc=None):
    global acache
    if (w, length) in acache:
        return acache[(w, length)]
    
    if length <= 1:
        return FA(*w.pad_code(length))
    if w.inv == 0:
        return FA(*w.pad_code(length))
    if rc is None:
        rc = RCGraph.principal_rc(w, length)
    bottom_rc, top_rc = rc.vertical_cut(length // 2)
    bperm, tperm = bottom_rc.perm, top_rc.perm
    left_schub = fast_schub(bperm, len(bottom_rc), bottom_rc)
    right_schub = fast_schub(tperm, len(top_rc), top_rc)
    base = left_schub * right_schub
    rc_prod = r(bottom_rc) * r(top_rc)
    #base = 
    for rc2, coeff in rc_prod.items():
        if rc2.perm != w:
            base -= coeff * fast_schub(rc2.perm, length)
    acache[(w, length)] = base
    return base

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    if n < 2:
        print("Please provide n >= 2")
        sys.exit(1)
    perms = Permutation.all_permutations(n)
    FSlideDual = FreeAlgebra(FundamentalSlideBasis)
    for w in perms:
        print("Testing result for", w)
        
        actual = ASx(w, n - 1).change_basis(FundamentalSlideBasis)
        test_base = ASx(w, n - 1).change_basis(WordBasis)
        test = 0
        seen_by_weight = {}
        for word, coeff in test_base.items():
            #fslideish = #r.from_dict({rc: coeff2 for rc, coeff2 in r.monomial(*word).items() if rc.is_quasi_yamanouchi})
            fslideish = sum([r(rc0) for rc0 in {_canonical_rc(rc) for rc in r.monomial(*word)}])
            
            for rc, coeff2 in fslideish.items():
                if rc.length_vector not in seen_by_weight:
                    seen_by_weight[rc.length_vector] = rc.perm_word
                elif seen_by_weight[rc.length_vector] != rc.perm_word:
                    continue
                test += coeff * coeff2 * FSlideDual(*rc.length_vector)
        assert test.almosteq(actual), f"Failed for {w}, got {test} but expected {actual}\n{test-actual=}"
        print(f"Passed for {w}")
    