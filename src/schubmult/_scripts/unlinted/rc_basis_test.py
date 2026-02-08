from schubmult import *
from sympy import pretty_print

r = RCGraphRing()

def _word_to_rc(word):
    r = RCGraphRing()

    ret = r.one

    for a in word:
        ret = ret * r(RCGraph.one_row(a))
    return ret

def rc_schubert_element(w, n):
    words = ASx(w, n).change_basis(WordBasis)
    
    ret = r.zero
    for word, coeff in words.items():
        ret += coeff * _word_to_rc(word)
    return ret

def rc_fa_elem(seq):
    schubs = FA(*seq).change_basis(SchubertBasis)
    moib = r.zero
    for (perm, length), coeff in schubs.items():
        moib += coeff * rc_schubert_element(perm, length)
    return moib

def prin_calc(w, n):
    if n == 0 and w.inv == 0:
        return r.one
    if n == 1 and len(w.trimcode) == 1:
        return r(RCGraph.one_row(w.trimcode[0]))
    if n < len(w.trimcode):
        return r.zero
    if w.inv == 0:
        return r(RCGraph.one_row(0))**n
    
    lowcode = w.trimcode[:-1]
    low_perm = uncode(lowcode)
    lower_prin = prin_calc(low_perm, n-1)
    bigger_sum = lower_prin * r(RCGraph.one_row(w.trimcode[-1]))
    new_sum = r.from_dict(bigger_sum)
    for rc, coeff in bigger_sum.items():
        if rc.perm != w:
            new_sum -= coeff * rc_fa_elem(rc.length_vector)
    return new_sum
        

if __name__ == "__main__":
    # Test the basis of the representation of the symmetric group on the cohomology of the regular semisimple Hessenberg variety.
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    for w in perms:
        
        #for ln in range(len(w.trimcode), n):
        ln = len(w.trimcode)
        rw = r.zero
        pretty_print(rc_schubert_element(w, ln))
        # assert len([v for v in rw if rw[v] != 0]) == 1
        # print(f"Yay {w,n}")
        # spoing = ASx(w, n).change_basis(WordhBasis)
        # for word, coeff in spoing.items():