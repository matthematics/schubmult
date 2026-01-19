from schubmult import *


def is_reduced(word):
    return Permutation.ref_product(*word).inv == len(word)

def w_i(word, i):
    return word[:i] + word[i + 1:]

def w_i_minus(word, i):
    letter = word[i]
    if letter > 1:
        return word[:i] + [letter - 1] + word[i + 1:]
    return [a + 1 for a in word[:i]] + [letter] + [b + 1 for b in word[i + 1:]]

def find_reduced_fail(word, inserted):
    perm = Permutation.ref_product(*word)
    a_start, b_start = perm.right_root_at(inserted, word=word)
    positive = False
    if a_start > b_start:
        positive = True
    for i in range(len(word)):
        if i == inserted:
            continue
        a, b = perm.right_root_at(i, word=word)
        if not positive and a > b:
            return i
        if positive and a == b_start and b == a_start:
            return i
    return None

def little_bump(word, i):
    assert is_reduced(w_i(word, i)) 
    new_word = w_i_minus(word, i)
    while not is_reduced(new_word):
        i = find_reduced_fail(new_word, i)
        new_word = w_i_minus(new_word, i)
    return new_word

if __name__ == "__main__":
    pass