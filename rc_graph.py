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
# perm0 = uncode([2, 3])
perm0 = uncode([2,0,4])

def rc_graph_set(perm, reorg_perm, floin):
    #if perm.inv == 0:
    if floin == 0:
        #print(f"{perm=}")
        return [tuple([(),(),(),()])]
    # ret = [tuple([(),()])]
    index = (~reorg_perm)[0]
    # indexes = ~reorg_perm
    ret = []
    L = pull_out_var(index, perm)
    #print(f"{index=}")
    for index_list, new_perm in L:
        #print(f"{perm=}, {new_perm=}")
        #print(f"{index_list=}")
        to_add_first = [*perm]
        #print(f"{to_add_first=}")
        # pos
        # to_add_first=["+" if a in index_list else "-" for a in to_add_first]
        #to_add_first=["+" if a + 1 in index_list else "-" for a in range(len(to_add_first))] # TEMP
        to_add_first=["+" if a + 1 in index_list else "-" for a in range(len(to_add_first))] # TEMP
        #print(f"{to_add_first=}")
        # print(f"{perm.code=}, {new_perm.code=}")
        # print(f"{((~new_perm)*perm).inv=} {perm.inv-new_perm.inv=}")
        #if ((~new_perm)*perm).inv != perm.inv-new_perm.inv:
        rc_set = rc_graph_set(new_perm, Permutation([r - 1 for r in reorg_perm if r != 1]), floin - 1)
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
            # new_term = term * prod([x[index] - y[il] for il in index_list])
            ret.append((new_labels, new_word, new_perml,new_grid))
    return ret
# KEY POLYNOMIAL OPERTAORS POSITIVE SCHUBERT
# Cauchy kernel
# WE CAN BACK OUT THE CAUCHY KERNEL FROM THE PULLOUTVAR
# SCHUB BASIS PULL OUT FACTOR
# FULL FACTORIZATION OF THE CAUCHY KERNEL
# nilPlactic key Cauchy kernel
# rc graph pull out var nilplactic
# identify Bruhat word

# rp = Permutation([3, 1, 2])  # reorganization permutation
from schubmult.rings import *

ring = NilHeckeRing(x)

def key_polynomial(comp):
    perm = Permutation.sorting_perm(comp, reverse=True)
    print(perm)
    part = [comp[perm[i] - 1 ] for i in range(len(comp))]
    return ring.isobaric(perm).apply(Sx(uncode(part))).as_polynomial().expand()

def perm_to_key(perm):
    if len(perm) == 0:
        return {NilPlactic(()): S.One}
    
    ret = {}
    stack = [(perm,(),1, S.One, [])]

    while len(stack) > 0:
        current_perm, word, index, poly, word2 = stack.pop()
        if current_perm.inv == 0:
            np_elem = word
            print(f"{(word, word2)}")
            # for i in range(2):
            #     print(f"{i+1}:{NilPlactic.reverse_insert_rsk(word, word2, '*', i + 1)}")
            # assert tuple(Plactic(word2)._word) == tuple(word2)
            ret[np_elem] = ret.get(np_elem, S.Zero) + poly
            continue        
        L = pull_out_var(1, current_perm)
        for index_list, new_perm in L:
            index_list.sort(reverse=True)
            new_word = [*word]
            new_word2 = [*word2]
            for index2 in index_list:
                new_word, new_word2 = NilPlactic.ed_insert_rsk(new_word, new_word2, index+index2 - 1, index)            
            stack.append((new_perm, new_word, index + 1, poly * prod([x[index] - y[a] for a in index_list]), new_word2))
    return ret

def rc_graph_eg(perm):
    # if len(perm) == 0:
    #     return (perm, {NilPlactic(()): S.One}
    
    ret = set()
    stack = [(perm,(),1, S.One, [])]

    while len(stack) > 0:
        current_perm, word, index, poly, word2 = stack.pop()
        if current_perm.inv == 0:
            print(f"{(word, word2)}")
            # for i in range(2):
            #     print(f"{i+1}:{NilPlactic.reverse_insert_rsk(word, word2, '*', i + 1)}")
            # assert tuple(Plactic(word2)._word) == tuple(word2)
            # ret[np_elem] = ret.get(np_elem, S.Zero) + poly
            ret.add((word, word2))
            continue        
        L = pull_out_var(1, current_perm)
        for index_list, new_perm in L:
            index_list.sort(reverse=True)
            new_word = [*word]
            new_word2 = [*word2]
            for index2 in index_list:
                new_word, new_word2 = NilPlactic.ed_insert_rsk(new_word, new_word2, index+index2 - 1, index)            
            stack.append((new_perm, new_word, index + 1, poly * prod([x[index] - y[a] for a in index_list]), new_word2))
    return ret

if __name__ == "__main__":
    # print(key_polynomial([3,0,5,1]))
    # free module over sym
    bob = Sx(uncode([0,2,1]))
    bob2 = bob * Sx(uncode([0,1]))
    print(perm_to_key(next(iter(bob.keys()))))
    print("prong")
    for key in bob2:
        print(perm_to_key(key))
    #Sx(uncode([0,2,1])).in_CEM_basis().expand()
    panko = uncode([0,2])
    graphs = rc_graph_eg(panko)
    # for i in range(1, 5):
    for graph in graphs:
        word1, word2 = graph
        word1 = [word1[i] - word2[i] + 1 for i in range(len(word1))]
        print(f"{(word1, word2)}")
    print("fat bacon")
    graphs = rc_graph_eg(~panko)
    # for i in range(1, 5):
    for graph in graphs:
        word1, word2 = graph
        word1 = [word1[i] - word2[i] + 1 for i in range(len(word1))]
        print(f"{(word1, word2)}")
    exit(0)
    # factor Cauchy kernel
    n = int(sys.argv[1])
    index = int(sys.argv[2])
    perms = [Permutation(perm) for perm in permutations([i for i in range(1, n+1)])]
    frog_dict = {}
    for perm in perms:
        L = pull_out_var(index, perm)
        for index_list, new_perm in L:
            # L2 = pull_out_var(index+1, new_perm)
            # for index_list2, new_perm2 in L2:
                #frog_dict[new_perm2] = frog_dict.get(new_perm2,S.Zero) + prod([x[index+1] - y[a] for a in index_list2])*prod([x[index] - y[a] for a in index_list])
            frog_dict[tuple(sorted(index_list))] = frog_dict.get(tuple(sorted(index_list)), []) + [new_perm]
    for pain in frog_dict:
        print(f"{pain}: {frog_dict[pain]}")    
    #print(frog_dict)
    exit(0)
    rp1 = Permutation([1,3,2])
    rp2 = Permutation([3,1,2])
    rs_set1 = rc_graph_set(perm0,~rp1,len(perm0))
    rs_set2 = rc_graph_set(perm0,~rp2,len(perm0))

    # SAME SPOT
    def print_bagel(rs_set,rp):
        for labels, word, perml,grid in rs_set:
            # print(f"{labels=}, {word=}")
            # print(f"{perml=}")
            print(f"{labels=}")
            term = prod([x[l] - y[ww] for l, ww in zip(labels, word)])
            print(f"{term=}")
            #top_line = set()
            bacon = []
            farple = []
            for i in range(len(grid)):
                grid2 = []
                #bump = 0
                spain = []
                boing = 0
                for j in range(len(grid[i])):
                    if grid[i][j] == "*":
                        grid2 = ["*"] + grid2
                        spain.append("*")
                        boing += 1
                    else:
                        #top_line.add(j + 1)
                        grid2.append(grid[i][j])
                        spain.append(perml[i][j - boing])
                farple.append(spain)
                bacon.append(grid2)
            bacon = grid
            #top_line = sorted(top_line)
            print("   "+" ".join([str(x) for x in perm0]))
            for i, grid2 in enumerate(bacon):
                print(f"{rp[i]}: "+" ".join([str(x) for x in grid2])+"  "+" ".join([str(x) for x in farple[i]]))
            print("")
            perm2 = Permutation([])
            flob=[]
            for i, word_s in enumerate(word):
                #perm2 = perm2.swap(word[rp[i]-1]+labels[rp[i]-1]-2, word[rp[i]-1]+labels[rp[i]-1]-1)
                flob += [word[i] + rp[labels[i]-1]-1]
                perm2 = perm2.swap(flob[-1]-1, flob[-1])
            # print(f"{perm2=}")
            # print(f"{flob=}")

    print_bagel(rs_set2,rp2)
    print("FROF")
    print_bagel(rs_set1,rp1)