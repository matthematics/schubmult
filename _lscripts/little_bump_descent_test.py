from schubmult import *

from schubmult.utils.perm_utils import little_bump_pos, find_reduced_fail, is_reduced

def test_descent_bump(word):
    perm = Permutation.ref_product(*word)
    word = list(word)
    
    i = max(perm.descents())
    if len([a for a in word if a == 1])>1:
        return
    if i == 0:
        return
    root = (i + 1, i + 2)
    index = next(iter([i for i in range(len(word)) if Permutation._right_root_at(i, word) == root]), None)
    if index is None:
        raise ValueError(f"Permutation {perm} does not have a root at {root}")
    if word[index] == 1:
        return
    modified_word = [*word[:index], word[index] - 1, *word[index + 1:]]
    index2 = find_reduced_fail(modified_word, index)
    if index2 is None:
        return
    if index2 >= index:
        raise ValueError(f"find_reduced_fail returned index {index2} which is greater than the original index {index} {word=} {modified_word=}")
    if word[index2] == 1 or Permutation._right_root_at(index2, word) == (1,2):
        return
    bumped_word = little_bump_pos(word, index)
    bumped_word2 = [*little_bump_pos(word[:index], index2), word[index] - 1, *word[index + 1:]]
    if bumped_word != tuple(bumped_word2):
        raise ValueError(f"Little bump at descent {i} failed: {bumped_word} != {tuple(bumped_word2)}")

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        if perm.inv == 0:
            continue
        if all(a > 0 for a in perm.trimcode):
            continue
        print(f"PPin {perm}")
        for word in perm.all_reduced_words():
            print(f"Word {word}")
            test_descent_bump(word)
            
        print(f"Potato {perm}")