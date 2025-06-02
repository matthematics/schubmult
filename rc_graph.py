from schubmult import *
from itertools import permutations
from schubmult.symbolic import *
from schubmult.schub_lib.schub_lib import pull_out_var
from schubmult.abc import *
import sys

# n = int(sys.argv[1])

# perms = [Permutation(perm) for perm in permutations([i for i in range(1,n+1)])]

# print("Testing CEM pattern criterion")

# for perm in perms:
#     cem_expansion = expand(Sx(perm).in_CEM_basis())
#     avoids_patterns = not (perm.has_pattern([1,4,3,2]) or perm.has_pattern([1,4,2,3]) or perm.has_pattern([4,1,3,2]) or perm.has_pattern([3,1,4,2]))
#     has_one_term = not (str(cem_expansion).find(" + ") != -1 or str(cem_expansion).find(" - ") != -1)
#     assert avoids_patterns == has_one_term

# print("Success")

# pull out gen generates RC graphs?
Permutation.print_as_code=True
perm0 = uncode([3,2,0,4])

def rc_graph_set(perm, reorg_perm):
    if perm.inv == 0:
        #print(f"{perm=}")
        return [tuple([(),(),(),()])]
    # ret = [tuple([(),()])]
    index = (~reorg_perm)[0]
    # indexes = ~reorg_perm
    ret = []
    L = pull_out_var(index, perm)
    print(f"{index=}")
    for index_list, new_perm in L:
        print(f"{perm=}, {new_perm=}")
        print(f"{index_list=}")
        to_add_first = [*perm]
        print(f"{to_add_first=}")
        to_add_first=["+" if a in index_list else "-" for a in to_add_first]
        print(f"{to_add_first=}")
        # print(f"{perm.code=}, {new_perm.code=}")
        # print(f"{((~new_perm)*perm).inv=} {perm.inv-new_perm.inv=}")
        #if ((~new_perm)*perm).inv != perm.inv-new_perm.inv:
        rc_set = rc_graph_set(new_perm, Permutation([r - 1 for r in reorg_perm if r != 1]))
        lsort = sorted(index_list, reverse=True)
        # lsort = sorted(index_list)
        

        # to_add_first[index] = "*"
        for labels, word, perml, grid in rc_set:
            #new_labels = tuple([1]*len(index_list)+[reorg_perm[label + 1 for label in labels])
            new_labels = tuple([index]*len(index_list)+[(label if label < index else label + 1) for label in labels])
            new_word = tuple(lsort+list(word))
            new_perml = tuple([perm,*perml])
            new_grid = [to_add_first]
            for i, perm_old in enumerate(grid):
                to_add = [*perm_old[:index - 1],"*",*perm_old[index - 1:]]
                to_add += ["-"]*(-(len(to_add)-len(perm)))
                new_grid.append(to_add)
            ret.append((new_labels, new_word, new_perml,new_grid))
    return ret

rp = Permutation([3, 1, 2])  # reorganization permutation
# rp = Permutation([])
rs_set = rc_graph_set(perm0,~rp)

for labels, word, perml,grid in rs_set:
    # print(f"{labels=}, {word=}")
    # print(f"{perml=}")
    print(f"{labels=}")
    for i in range(len(grid)):
        grid2 = []
        for j in range(len(grid[i])):
            if grid[i][j] == "*":
                grid2 = ["*"] + grid2
            else:
                grid2.append(grid[i][j])
        print(f"{rp[i]}: "+" ".join([str(x) for x in grid2]))
    print("")
    perm2 = Permutation([])
    flob=[]
    for i, word_s in enumerate(word):
        #perm2 = perm2.swap(word[rp[i]-1]+labels[rp[i]-1]-2, word[rp[i]-1]+labels[rp[i]-1]-1)
        flob += [word[i] + rp[labels[i]-1]-1]
        perm2 = perm2.swap(flob[-1]-1, flob[-1])
    # print(f"{perm2=}")
    # print(f"{flob=}")
