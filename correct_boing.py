from schubmult import *
from itertools import permutations
from schubmult.symbolic import *
from schubmult.schub_lib.schub_lib import pull_out_var
from schubmult.abc import *
from schubmult.utils.perm_utils import has_bruhat_ascent
import sys

def shiftup(pp, i):
    if i == 0:
        return pp
    return shiftup(Permutation([1] + [poip+1 for poip in pp]), i - 1)

def reduced_word(perm):
    if perm.inv == 0:
        return []
    to_ret = []
    for iasd in range(len(perm)-1):
        if perm[iasd] > perm[iasd + 1]:
            perm = perm.swap(iasd, iasd + 1)
            to_ret = [iasd+1] + to_ret
    return [*reduced_word(perm), *to_ret]

def root_at(word, index):
    aa = word[index]
    bb = word[index] + 1
    for i in range(index-1, -1, -1):
        aa = Permutation.ref_product(word[i])[aa - 1]
        bb = Permutation.ref_product(word[i])[bb - 1]
    return (aa, bb)

def root_at_right(word, index):
    aa = word[index]
    bb = word[index] + 1
    for i in range(index + 1, len(word)):
        aa = Permutation.ref_product(word[i])[aa - 1]
        bb = Permutation.ref_product(word[i])[bb - 1]
    return (aa, bb)



def pull_out_var_reflections(vnum, v):
    v = Permutation(v)
    vup = v
    if vnum >= len(v):
        return [[[], v, (), []]]
    vpm_list = [(vup, 0, [])]
    ret_list = []
    for p in range(len(v) + 1 - vnum):
        vpm_list2 = []
        for vpm, b, refs in vpm_list:
            if vpm[vnum - 1] == len(v) + 1:
                vpm2 = [*vpm]
                vpm2.pop(vnum - 1)
                vp = permtrim(vpm2)
                ret_list += [
                    [
                        [v[i] for i in range(vnum, len(v)) if ((i > len(vp) and v[i] == i) or (i <= len(vp) and v[i] == vp[i - 1]))],
                        vp, vpm, refs,
                    ],
                ]
            for j in range(vnum, len(vup) + 2):
                if vpm[j] <= b:
                    continue
                for i in range(vnum):
                    if has_bruhat_ascent(vpm, i, j):
                        vpm_list2 += [(vpm.swap(i, j), vpm[j], refs + [(i + 1, j + 1)] )]
        vpm_list = vpm_list2
    for vpm, b, refs in vpm_list:
        if vpm[vnum - 1] == len(v) + 1:
            vpm2 = [*vpm]
            vpm2.pop(vnum - 1)
            vp = permtrim(vpm2)
            ret_list += [
                [
                    [v[i] for i in range(vnum, len(v)) if ((i > len(vp) and v[i] == i) or (i <= len(vp) and v[i] == vp[i - 1]))] , vp, vpm, refs,
                    
                ],
            ]
    return ret_list

def code_word(perm):
    cd = perm.code
    word = []
    for i, code_value in enumerate(cd):
        word += list(range(code_value + i,i,-1))
    return word

def normal_word(perm, A):
    # cd = perm.code
    # word = []
    # for i, code_value in enumerate(cd):
    #     word += list(range(code_value + i,i,-1))
    # return word
    if len(A) == 0:
        return perm
    word = []
    found = True
    A2 = [*A[1:]]
    index = A[0]
    perm2 = perm
    while found:
        found = False
        for i in range(len(perm2) - 1):
            if i == index - 1:
                continue
            if perm2[i] > perm2[i + 1]:
                found = True
                perm2 = perm2.swap(i, i+1)
    return [*code_word(perm2), *normal_word((~perm2)*perm, A2)]



def delete_ref_from_word(word, ref):
    # print(f"{word=}")
    # print(f"{ref=}")
    # print(f"{[root_at(word, i) for i in range(len(word))]}")
    for i in range(len(word)):
        if root_at_right(word, i) == ref:
            return [*word[:i],*word[i+1:]]
    
    raise ValueError("can't do")


if __name__ == "__main__":

    n = int(sys.argv[1])

    perms = [Permutation(perm) for perm in permutations([i for i in range(1,n+1)])]

    
    for perm in perms:
        if perm.inv == 0:
            continue
        # if max(perm.descents()) < n-2:
        #     continue

        # if max(perm.descents()) < A[-1]:
        #     continue
        dd = max(perm.descents()) + 1
        #A = list(range(1,dd-1))
        if dd != 5:
            continue
        #A = [2, 3, 5, 1, 4]
        # A = [5,4,3,2,1]
        A = [5,4,1,2,3]
        n = len(perm) + 1
        Aspots = ~Permutation(A)
        print(f"{perm=}")
        end_perm = Permutation([n + 1 - a for a in Aspots] + list(range(1,n-4)))
        print(f"{end_perm=}")
        word = code_word(end_perm)
        assert Permutation.ref_product(*word) == end_perm
        results = [[]]
        fatuous_results = [[]]
        perm_stack = [perm]
        weight_stack = [[]]
        end_stack = [perm]
        index_stack = [[]]
        for a_index, varspot in enumerate(A):            
            print(f"{a_index=} {varspot=}")

            new_results = []
            new_perm_stack = []
            new_weight_stack = []
            new_fatuous_results = []
            new_end_stack = []
            new_index_stack = []
            while len(perm_stack) > 0:
                new_refs = results.pop()
                new_perm = perm_stack.pop()
                new_weight = weight_stack.pop()
                new_indexes = index_stack.pop()
                #the_index = varspot - len([bi for bi in A[:a_index] if bi < varspot])
                # refpulllist = pull_out_var_reflections(the_index, new_perm)
                # for index_list, new_new_perm, fat_perm, refs in refpulllist:
                vup = new_perm # this is a fat perm
                # if vvarspot >= len(new_perm):
                #     return [[[], v, (), []]]
                # target = max([arb for arb in new_perm.rslice(varspot - 1, len(new_perm))]) + 1
                # print(f"{target=}")
                target = n - a_index #max(max([arb for arb in new_perm.rslice(varspot - 1, len(new_perm))]),len(perm)) + 1

                # if (~new_perm)[target] <= varspot:
                #     while (~new_perm)[target] <= varspot:
                #         target += 1

                vpm_list = [(vup, 0, [])]
                ret_list = []
                #for p in (range(target - varspot)):
                while len(vpm_list) > 0:
                    vpm_list2 = []
                    for vpm, b, refs in vpm_list:
                        if vpm[varspot - 1] > target:
                            continue
                        if vpm[varspot - 1] == target:
                            
                            #print(f"Found {vpm}")
                            # vpm2 = [*vpm]
                            # vpm2.pop(vnum - 1)
                            vp = permtrim(vpm)
                            # ret_list += [
                            #     [
                            #         [v[i] for i in range(vnum, len(v)) if ((i > len(vp) and v[i] == i) or (i <= len(vp) and v[i] == vp[i - 1]))],
                            #         vp, vpm, refs,
                            #     ],
                            # ]
                            # the_target = n - a_index
                            the_target = target
                            vvp = Permutation([arb for arb in vp if arb < the_target])
                            print(f"Found {vp=} {vvp=}")
                            new_perm0 = Permutation([arb for arb in new_perm if arb <= the_target])
                            index_list = [new_perm0[i] for i in range(varspot - len([aa for aa in A[:a_index] if aa < varspot]), len(new_perm0)) if ((i > len(vvp) and new_perm0[i] == i and new_perm0[i]<target) or (i <= len(vvp) and new_perm0[i] == vvp[i - 1] and new_perm0[i] < target))]
                            new_perm_stack.append(vp)
                            new_results.append([*new_refs,*refs])
                            new_new_weight = new_weight + ([varspot] * len(index_list))
                            new_weight_stack.append(new_new_weight)
                            #new_index_stack.append([*new_indexes,*list(sorted(index_list,reverse=True))])
                            new_index_stack.append([*new_indexes,*list(sorted(index_list))])
                 
                        for j in range(varspot, len(vup) + 2):
                            if vpm[j] <= b or vpm[j]>target or j + 1 in A[:a_index]:
                                continue
                            for i in range(varspot):
                                if i + 1 in A[:a_index] or vpm[i] > target:
                                    continue
                                if has_bruhat_ascent(vpm, i, j):
                                    vpm_list2 += [(vpm.swap(i, j), vpm[j], refs + [(i + 1, j + 1)] )]
                    vpm_list = vpm_list2
                for vpm, b, refs in vpm_list:
                    if vpm[varspot - 1] == target:
                        vp = permtrim(vpm)
                        # ret_list += [
                        #     [
                        the_target = n - a_index
                        vvp = Permutation([arb for arb in vp if arb < the_target])
                        new_perm0 = Permutation([arb for arb in new_perm if arb <= the_target])
                        index_list = [new_perm0[i] for i in range(varspot - len([aa for aa in A[:a_index] if aa < varspot]), len(new_perm0)) if ((i > len(vvp) and new_perm0[i] == i and new_perm0[i]<target) or (i <= len(vvp) and new_perm0[i] == vvp[i - 1] and new_perm0[i] < target))]
                                
                        #     ],
                        new_results.append([*new_refs,*refs])
                 
                        new_new_weight = new_weight + ([varspot] * len(index_list))
                        new_weight_stack.append(new_new_weight)
                        #new_index_stack.append([*new_indexes,*list(sorted(index_list,reverse=True))])
                        new_index_stack.append([*new_indexes,*list(sorted(index_list))])
                    

                    
            fatuous_results = new_fatuous_results
            results = new_results
            perm_stack = new_perm_stack
            weight_stack = new_weight_stack
            end_stack = new_end_stack
            index_stack = new_index_stack

        print(f"printing {perm=}")
        print(f"{results=}")

        if len(results) == 0:
            print("Fail")
        
        for the_word, the_weight, the_perm, index_weight in zip(results, weight_stack, perm_stack, index_stack):
            end_perm = the_perm
            
            # print(f"{the_perm=}")
            # print(f"{the_word=}")
            #print(f"{the_fat=}")
            # for bongle in the_word:
            #     frof = end_perm.inv
            #     print(f"{bongle=}")
            #     print(f"{end_perm=}")
            #     end_perm = end_perm.swap(bongle[0] - 1, bongle[1] - 1)
            #     assert frof == end_perm.inv - 1
            
            word = normal_word(end_perm, A)

            #word = code_word(end_perm)
            
            start_word = [*word]
            # print(f"{end_perm=}")
            # print(f"{end_perm.inv=} {len(the_word)=} {perm.inv=}")
            new_end_perm = Permutation.ref_product(*start_word)
            # print(f"{new_end_perm=}")
            assert new_end_perm == end_perm
            # fat_perm = end_perm
            # print(f"{the_word=}")
            # print(f"{fat_perm=}")
            # for aa, bb in reversed(the_word):
            #     print(f"{aa,bb=}")
            #     elv = fat_perm.inv
            #     fat_perm = fat_perm.swap(aa-1, bb-1)
            #     assert elv == fat_perm.inv + 1
            #     print(f"{fat_perm=}")
            # print(f"{fat_perm=} {perm=}")
            # assert fat_perm == perm
            
            for r in reversed(the_word):
                start_word = delete_ref_from_word(start_word, r)
            print("--------")
            the_weight.sort()
            print(the_weight)
            index_weight.sort(reverse=True)
            print(index_weight)
            print(start_word)

    exit(0)
