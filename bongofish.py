from schubmult import *
from itertools import permutations
from schubmult.symbolic import *
from schubmult.schub_lib.schub_lib import pull_out_var
from schubmult.abc import *
import sys

n = int(sys.argv[1])

perms = [Permutation(perm) for perm in permutations([i for i in range(1,n+1)])]

# print("Testing CEM pattern criterion")

# for perm in perms:
#     cem_expansion = expand(Sx(perm).in_CEM_basis())
#     avoids_patterns = not (perm.has_pattern([1,4,3,2]) or perm.has_pattern([1,4,2,3]) or perm.has_pattern([4,1,3,2]) or perm.has_pattern([3,1,4,2]))
#     has_one_term = not (str(cem_expansion).find(" + ") != -1 or str(cem_expansion).find(" - ") != -1)
#     assert avoids_patterns == has_one_term

# print("Success")

# perm_index = {}

# for i, perm in enumerate(perms):
#     perm_index[perm] = i

# perms.sort(key=lambda p: p._unique_key)

# bongset = set()
# for pingus1 in perms:
#     if pingus1.inv == 0:
#         continue
#     for pingus2 in perms:
#         if pingus2.inv == 0:
#             continue
#         bd = (DSx(pingus1) * DSx(pingus2, "z")).expand(deep=False)
#         bcode1 = (~((~pingus1).minimal_dominant_above())).code
#         bcode2 = (~((~pingus2).minimal_dominant_above())).code

#         bcode = [*bcode1]
#         for i, b2 in enumerate(bcode2):
#             if i >= len(bcode):
#                 bcode += [b2]
#             else:
#                 bcode[i] += b2
#         luv = uncode(bcode)
#         for w in bd:
#             if len(w)<=n:
#                 bongset.add((w * (~luv))*(~(w * (~luv))).minimal_dominant_above())

# reduced?


#A = [2, 3, 4]
A = [3, 2, 1]
for perm in perms:
# for afs in range(1):
#     perm=Permutation((1, 4, 2, 5, 3, 6, 8, 7))
    bongset = set()
    if perm.inv == 0:
        continue
    if max(perm.descents()) < A[-1]:
        continue
    if any((d>=min(A) - 1 and d<=max(A)-1 for d in perm.descents())):
        continue
    refs = []
    seq = []
    perm_stack = [perm]
    ref_stack = [[]]
    seq_stack = [[]]
    for i, a in enumerate(A):
        new_perm_stack = []
        new_ref_stack = []
        new_seq_stack = []
        while len(perm_stack) != 0:
            perm0 = perm_stack.pop()    
            seq0 = seq_stack.pop()
            ref0 = ref_stack.pop()
            for index_list, new_perm in pull_out_var(a - sum([1 for bi in A[:i] if bi < a]), perm0):
                #fesspool = perm0.pattern_at(*list(range(a - i - 1, len(perm0))))
                #print(f"{fesspool=}")
                new_perm_stack.append(new_perm)
                #index_list_sort = [fesspool[(~perm0)[j-1] - a + i] for j in index_list]
                index_list_sort = [*index_list]
                #print(f"{index_list=}")
                #print(f"{index_list_sort=}")
                index_list_sort.sort(reverse=True)
                ref_add = [i+i0 for i0 in index_list_sort]
                new_ref_stack.append(ref0 + ref_add)
                new_seq_stack.append(seq0 + [i+1]*(len(index_list)))
        perm_stack = new_perm_stack
        ref_stack = new_ref_stack
        seq_stack = new_seq_stack
    print(f"{perm.code=}")
    for ref, seq, new_perm in zip(ref_stack,seq_stack,perm_stack):
        print(f"{new_perm=}")
        print(seq)
        print(ref)
        flob = Permutation([])
        boing = NilPlactic.from_word([])
        for r in ref:
            if flob[r-1] > flob[r]:
                print("Fail")
            flob = flob.swap(r-1, r)
            boing = boing.ed_insert(r)
        print(f"{boing=}")
        shape = [len(b) for b in boing._word]
        from schubmult.utils.perm_utils import p_trans
        shape = p_trans(shape)
        #shape.reverse()
        bongset.add((tuple(shape),tuple(new_perm.code)))
    print(bongset)
        

# perm=(1, 4, 2, 5, 3, 6, 8, 7)
# new_perm=()
# [1, 1, 1, 3]
# [7, 3, 2, 4]

# new_perm=(2, 1)
# [1, 1, 1]
# [7, 3, 2]

# new_perm=()
# [1, 1, 1, 2]
# [7, 3, 2, 4]

# new_perm=(2, 1)
# [1, 1, 3]
# [7, 3, 3]

# Fail     
# semistandard coxeter braid group       
        


# for pingus in perms:
#     bongset.add(pingus*(~pingus).minimal_dominant_above())
exit(0)
for bong in bongset:
    print(bong.code)

# into the product poset of the maximal parabs

# pull out gen generates RC graphs?
# Permutation.print_as_code=True
# perm0 = uncode([3,2,0,4])

# def rc_graph_set(perm, reorg_perm):
#     if perm.inv == 0:
#         #print(f"{perm=}")
#         return [tuple([(),(),(),()])]
#     # ret = [tuple([(),()])]
#     index = (~reorg_perm)[0]
#     # indexes = ~reorg_perm
#     ret = []
#     L = pull_out_var(index, perm)
#     print(f"{index=}")
#     for index_list, new_perm in L:
#         print(f"{perm=}, {new_perm=}")
#         print(f"{index_list=}")
#         to_add_first = [*perm]
#         print(f"{to_add_first=}")
#         to_add_first=["+" if a in index_list else "-" for a in to_add_first]
#         print(f"{to_add_first=}")
#         # print(f"{perm.code=}, {new_perm.code=}")
#         # print(f"{((~new_perm)*perm).inv=} {perm.inv-new_perm.inv=}")
#         #if ((~new_perm)*perm).inv != perm.inv-new_perm.inv:
#         rc_set = rc_graph_set(new_perm, Permutation([r - 1 for r in reorg_perm if r != 1]))
#         lsort = sorted(index_list, reverse=True)
#         # lsort = sorted(index_list)
        

#         # to_add_first[index] = "*"
#         for labels, word, perml, grid in rc_set:
#             #new_labels = tuple([1]*len(index_list)+[reorg_perm[label + 1 for label in labels])
#             new_labels = tuple([index]*len(index_list)+[(label if label < index else label + 1) for label in labels])
#             new_word = tuple(lsort+list(word))
#             new_perml = tuple([perm,*perml])
#             new_grid = [to_add_first]
#             for i, perm_old in enumerate(grid):
#                 to_add = [*perm_old[:index - 1],"*",*perm_old[index - 1:]]
#                 to_add += ["-"]*(-(len(to_add)-len(perm)))
#                 new_grid.append(to_add)
#             ret.append((new_labels, new_word, new_perml,new_grid))
#     return ret

# rp = Permutation([3, 1, 2])  # reorganization permutation
# # rp = Permutation([])
# rs_set = rc_graph_set(perm0,~rp)

# for labels, word, perml,grid in rs_set:
#     # print(f"{labels=}, {word=}")
#     # print(f"{perml=}")
#     print(f"{labels=}")
#     for i in range(len(grid)):
#         print(f"{rp[i]}: "+" ".join([str(x) for x in grid[i]]))
#     print("")
#     perm2 = Permutation([])
#     flob=[]
#     for i, word_s in enumerate(word):
#         #perm2 = perm2.swap(word[rp[i]-1]+labels[rp[i]-1]-2, word[rp[i]-1]+labels[rp[i]-1]-1)
#         flob += [word[i] + rp[labels[i]-1]-1]
#         perm2 = perm2.swap(flob[-1]-1, flob[-1])
#     # print(f"{perm2=}")
#     # print(f"{flob=}")
