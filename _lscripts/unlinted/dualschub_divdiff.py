from schubmult import *
from schubmult.rings.free_algebra import FreeAlgebraElement, FA

def dual_schub_divdiff(free_algebra_element, index):
    elem = FA.from_dict(free_algebra_element.change_basis(WordBasis))
    # [b, a]
    ret_elem = 0
    i = index - 1
    for word, coeff in elem.items():
        new_elem = FA(*word[:i])
        if index + 1 > len(word):
            continue
        
        a, b = word[i], word[i + 1]
        if a >= b:
            new_elem = new_elem * sum([(FA(a - 1 - j, b + j)  - FA(b + j, a - 1 - j)) for j in range(a)])
        else:
            new_elem = new_elem * sum([(FA(a + j, b - 1 - j)  - FA(b - 1 - j, a + j)) for j in range(b)])
        new_elem *= FA(*word[i+2:])
        if a < b - 1:
            new_elem = -new_elem
        ret_elem += coeff * new_elem
    if ret_elem == 0:
        return 0
    return ret_elem.change_basis(free_algebra_element.ring._basis)

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    for perm in perms:
        if perm.inv == 0:
            continue
        dual_schub = ASx(perm, n - 1)
        for index in range(1, n):
            divdiffed = dual_schub_divdiff(dual_schub, index)
            print(f"{perm.trimcode}:d{index} => {dual_schub} => {divdiffed}")
