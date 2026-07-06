from schubmult import *
from schubmult.symbolic import prod
from schubmult.utils.tuple_utils import pad_tuple

def reduced_set_valued_compatible_sequences(perm: Permutation):
    ret = []
    for wc in WCGraph.all_wc_graphs(perm):
        increasing, recording = IncreasingTableau.hecke_column_insert_rsk(wc.compatible_sequence, wc.perm_word)
        assert tuple(pad_tuple(recording.weight[:len(wc.length_vector)], len(wc))) == wc.length_vector, f"Length vector mismatch for {wc.perm.trimcode}: {tuple(recording.weight)} != {wc.length_vector}"
        ret.append((increasing, recording))
    return ret
    # stack = [(rc.perm_word, tuple([frozenset({a}) for a in rc.compatible_sequence])) for rc in RCGraph.all_rc_graphs(perm)]
    # ret = set()
    # #seen = set()
    # while stack:
    #     perm_word, comp_seq = stack.pop()
    #     if (perm_word, comp_seq) not in ret:
    #         ret.add((perm_word, comp_seq))
    #         #seen.add((perm_word, comp_seq))
    #     for i in range(len(perm_word)):
    #         max_val = perm_word[i]
    #         min_val = 1
    #         if i < len(perm_word) - 1:
    #             max_val = min(max_val, min(comp_seq[i + 1]))
    #             if perm_word[i] < perm_word[i + 1]:
    #                 if max_val == min(comp_seq[i + 1]):
    #                     max_val = min(comp_seq[i + 1]) - 1
    #         if i > 0:
    #             min_val = max(min_val, max(comp_seq[i - 1]))
    #             if perm_word[i] > perm_word[i - 1]:
    #                 if min_val == max(comp_seq[i - 1]):
    #                     min_val = max(comp_seq[i - 1]) + 1
    #         if min_val > max_val or min_val <= 0 or max_val <= 0:
    #             continue
    #         for new_val in range(min_val, max_val + 1):
    #             if new_val not in comp_seq[i]:
    #                 new_comp_seq = list(comp_seq)
    #                 st = set(comp_seq[i])
    #                 st.add(new_val)
    #                 new_comp_seq[i] = frozenset(st)
    #                 #ret.add(tuple(new_comp_seq))
    #                 stack.append((perm_word, tuple(new_comp_seq)))
    # return ret

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        svcs = reduced_set_valued_compatible_sequences(perm)
        poly = 0
        for seq11 in svcs:
            toadd = 1
            word = seq11[0].row_word
            if len(word) > perm.inv:
                print(f"{word=}")
            seq = seq11[1].weight
            for i, s in enumerate(seq):
                if s == 0:
                    continue
                toadd *= prod([Gx.genset[i + 1] ** s])
            poly += Gx._beta**(sum(seq) - perm.inv) * toadd
        diff = (poly - Gx(perm).expand()).expand() 
        assert diff == 0, f"Grothendieck polynomial for {perm} does not match sum over reduced set-valued compatible sequences, \n{diff=}"    
        #print(f"{perm.trimcode}: {len(svcs)} reduced set-valued compatible sequences")
        print("Pataad")