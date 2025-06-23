from schubmult import *
from itertools import permutations
from schubmult.symbolic import *
from schubmult.schub_lib.schub_lib import pull_out_var
from schubmult.abc import *
import sys

n = int(sys.argv[1])

perms = [Permutation(perm) for perm in permutations([i for i in range(1,n+1)])]
from schubmult.rings import SingleSchubertRing

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


# def test_new_pull_out_var(index, perm):
#     print(f"{index=}")
#     tpoly = Sx(perm).as_polynomial()
#     L = pull_out_var(index, perm)
#     porn_set = set()
#     gset = MaskedGeneratingSet(x, [index])
#     SRing = SingleSchubertRing(gset)
#     test_poly = S.Zero
#     print(f"{perm=}")
#     for index_list, new_perm in L:
#         print(type(new_perm))
#         lnp = list(new_perm)
#         if len(lnp) == 0:
#             lnp = [1]
#         print(f"{lnp=} {new_perm[:index-1]=} {new_perm[index-1:max(lnp)]=} {new_perm=}")
#         reset_perm = Permutation(list(new_perm.rslice(0,index-1))+[max(list(perm))]+list(new_perm.rslice(index-1,(~new_perm)[max(list(perm))]-2)))
#         print(f"{reset_perm=}")
#         pperm = (~perm) * reset_perm
        
#         cycs = pperm.get_cycles()
#         print(f"{pperm=} first {cycs=}")
#         cycs = [c for c in pperm.get_cycles() if index not in c]
#         cycs = [tuple([cc - 1 for cc in c]) for c in cycs]
#         print(f"{cycs=}")
#         reset_perm = reset_perm * ~Permutation.from_cycles(cycs)
#         tuff_perm = Permutation(list(reset_perm[:index-1])+list(reset_perm[index:]))
#         porn_set.add(tuff_perm)
#         print(f"{new_perm=} {tuff_perm=} {perm=}")
#     for permo in porn_set:
#         print(f"{permo=}")
#         amt = perm.inv - permo.inv
#         sporkle = SRing.from_dict({permo: 1}).as_polynomial()*expand_func(h(amt, index - 1,x[1:]))
#         print(f"{sporkle=}")
#         test_poly += sporkle
#     print(f"{perm=} {expand(test_poly-tpoly)=}")
#     print(f"{test_poly=} {tpoly=}")
    
# for i in range(2, n):
#     for perm in perms:
#         if perm.inv == 0:
#             continue
#         test_new_pull_out_var(i, perm)

# simple refs up, H refs down!
# vt_1 .. t_k(2 through n-1 cycle) = v'
# v t_{1,b_1} .. t_{i, b_i} = (v'(n-1) through v'(i) cycle)

# v (<i reflections) = (s_n through s_i, spots deleted)v' patternmap_removei(v')
# v  = (s_n through s_i, spots deleted)v' patternmap_removei(v') (remove h reflections)
# this allows RC graph building, unbiject reflections
# weight preserving bijections for fixed v between (power of x_i) x RC graphs for larger pull out perms -> (simple refs) x (excess h weights for <=i -1, length is strip length + actual pattern - length - perm length) x (RC graph for reflections in the remainder pattern)
# n down to i simple ref strip piece x patternmap(v'), length additive in pattern, number of reflections that go beneath is excess weight < i, remainder is i weight

# Schubert poly with n variables: Permutation \sigma of [n] x (\sigma-patternized reflection strips, length additive multiply to v) x (weights <= sigma i)

# print("Testing CEM pattern criterion")

# v (<i reflections) (i reflections) = v' (n to i simple refs)

# v (<i reflections) = (n to i simple refs) conj(v') 
# 
# (n to i simple refs)patternmap(v') = l(v') + n - i

# v (<i reflections) = (n to i simple refs with dels)patternmap(v') l(v) + excess = l(v') + l(ref strip)
# l(v) - l(v') number of is
# v = (ref strip)patternmap(v'') l(v'')

# maximal reduced word pulling backward
# peel off the smaller elements?
# for perm in perms:
#     cem_expansion = expand(Sx(perm).in_CEM_basis())
#     avoids_patterns = not (perm.has_pattern([1,4,3,2]) or perm.has_pattern([1,4,2,3]) or perm.has_pattern([4,1,3,2]) or perm.has_pattern([3,1,4,2]))
#     has_one_term = not (str(cem_expansion).find(" + ") != -1 or str(cem_expansion).find(" - ") != -1)
#     assert avoids_patterns == has_one_term

# reduced word -> i's -> excess weights which are less

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



#A = [n - 3, n - 2, n - 1]
#Aperm = Permutation([1, 2, 3, 4, 5])
for perm in perms:
# for afs in range(1):
#     perm=Permutation((1, 4, 2, 5, 3, 6, 8, 7))
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
    A = [3, 4, 1, 2, 5]
    # if len(A) == 0:
    #     continue
    # if dd<3 or any((d>=min(A) - 1 and d<max(A)-1 for d in perm.descents())):
    #     continue
    poly = Sx(perm)
    print(f"{dd=} {A=}")
    polcop = poly.coproduct(*A)
    songset = set(polcop.keys())
    bongset = set()
    bongdict = {}
    subs_dict = {x[i]: S.One for i in range(10)}
    songdict = {(k1, k2): (Sx(k1).as_polynomial().subs(subs_dict))*v for (k1, k2), v in polcop.items()}
    refs = []
    seq = []
    perm_stack = [perm]
    ref_stack = [Permutation([])]
    bongo_stack = [[]]
    seq_stack = [[]]
    #A = Aperm[:max(len(Aperm),max(perm.descents()) + 1)]
    #A = [a for a in A if a <= max(perm.descents()) + 1]
    #print(f"{A=}")
    #TEST
    #A
    used = []
    #flipdunk = len(A) - 1
    for i, a in enumerate(A):
        #i = flipdunk - j
        used.append(a)
        new_perm_stack = []
        new_ref_stack = []
        new_seq_stack = []
        new_bongo_stack = []
        while len(perm_stack) != 0:
            perm0 = perm_stack.pop()    
            seq0 = seq_stack.pop()
            ref0 = ref_stack.pop()
            bongo0 = bongo_stack.pop()
            the_index = a - len([bi for bi in used if bi < a])
            for index_list, new_perm in pull_out_var(a - len([bi for bi in used if bi < a]), perm0):#pull_out_var(a - sum([1 for bi in A[:i] if bi < a]), perm0):
                #fesspool = perm0.pattern_at(*list(range(a - i - 1, len(perm0))))
                #print(f"{fesspool=}")
                # if new_perm.rslice(0, the_index - 1) != perm0.rslice(0, the_index - 1):
                #     continue
                #new_perm_stack.append(new_perm)
                #index_list_sort = [fesspool[(~perm0)[j-1] - a + i] for j in index_list]
                index_list_sort = [*index_list]
                #print(f"{index_list=}")
                #print(f"{index_list_sort=}")
                index_list_sort.sort(reverse=True)
                #ref_add = [i+i0 for i0 in index_list_sort]
                #TEST
                # big_n = len(perm0)
                # perm_cycle = Permutation([big_n + 1, *list(range(1,big_n + 1))])
                # # print(f"{big_n=}")
                # # new_new_perm = new_perm
                # # for ii in range(big_n - 2, -1, -1):
                # #     print(f"{ii=} {new_new_perm=}")
                # #     new_new_perm = new_new_perm.swap(ii, ii+1)
                # #     print(f"{ii=} {new_new_perm=}")
                # # cycle =((~new_new_perm) * perm0) #perm_cycle * (~((~perm0) * (new_new_perm)))
                # #cycle = ((~new_perm) * perm0)
                # #cycle = (~)
                # new_new_perm = [big_n + 1, *list(new_perm)]
                # mx = max(new_perm) if new_perm.inv != 0 else 0
                # while len(new_new_perm) < big_n + 1:
                #     new_new_perm += [mx + 1]
                #     mx += 1
                # new_new_perm = Permutation(new_new_perm)
                # print(f"{new_new_perm=}")
                # cycle = ((~perm0) * new_new_perm)
                # #((~new_new_perm) * perm0) * perm_cycle
                # if cycle.inv != len(index_list):
                #     print(f"{cycle=}")
                #     # print(f"{perm_cycle=}")
                #     # print(f"{new_perm*perm_cycle=}")
                #     print(f"{new_new_perm=}")
                #     print(f"{perm0=}")
                #     print(f"{new_perm=}")
                #     print(f"{index_list=}")
                #     exit(0)

                # ref_add = []
                # while cycle.inv != 0:
                #     md = max(cycle.descents())
                #     ref_add += [md + 1]
                #     cycle = cycle.swap(md, md+1)
                #ref_add = list([(~perm0)[i0 - 1] + i - 1 for i0 in index_list_sort])
                #ref_add = list(reversed([j + 1 for j in range(len(perm0)) if j > the_index - 1 and perm0[j] == new_perm[j-1]]))
                # (14)1 3 5 4 6 2
                
                #ref_add = list([i0 + i for i0 in index_list_sort])
                # ref_add = list([i0 + i for i0 in index_list_sort])
                # p_perm_inv = perm0.inv - (perm0.inv - new_perm.inv)
                # cyc = perm0 * (~Permutation([1] + [j+1 for j in new_perm]))#(~Permutation.ref_product(*list(range(p_perm_inv + 1 + i,1 + i,-1))))
                # ref_add = []
                # for joi in range(len(cyc) - 1):
                #     if cyc[joi] > cyc[joi + 1]:
                #         ref_add =  [joi + 1 + i] + ref_add
                #         cyc = cyc.swap(joi, joi + 1)
                #cyc = Permutation.ref_product(*index_list_sort)
                full_ref_prod = Permutation.ref_product(*list(range(max(the_index+2,len(perm0)+1), 0, -1)))
                shifted_new_perm = Permutation([1] + [jj+1 for jj in new_perm])
                lower_ref_prod = Permutation.ref_product(*list(range(1, the_index)))
                up_perm = new_perm * (full_ref_prod*(lower_ref_prod))
                # up_perm = new_perm * (full_ref_prod)
                build_perm = (up_perm * (~((shifted_new_perm))))
                
                print(f"{up_perm=}")
                print(f"{shifted_new_perm=} {the_index=}")
                print(f"{reduced_word(build_perm)=}")
                
                if build_perm.inv != up_perm.inv - new_perm.inv:
                    #raise Exception
                    print("Big guns")
                    # up_perm_fat = full_ref_prod
                    build_perm = full_ref_prod
                    # lower ref, (1, the_index), .. (1,2)
                    for bp in range(the_index, 1, -1):
                        if build_perm[shifted_new_perm[0] - 1] < build_perm[shifted_new_perm[bp - 1] - 1]:
                            print(f"{shifted_new_perm[0]=} {shifted_new_perm[bp - 1]=} {bp=}")
                            print(f"{build_perm=}")
                            shifted_new_perm = shifted_new_perm.swap(0, bp - 1)
                            wordo = reduced_word(build_perm)
                            print(f"{wordo=}")
                            wordo.reverse()
                            for ii in range(len(wordo)):
                                aa, bb = root_at(wordo, ii)
                                isnp = ~shifted_new_perm
                                if isnp[aa - 1] > isnp[bb - 1]:
                                    print(f"Found {aa=} {bb=}")
                                    isnp = isnp.swap(aa-1, bb-1)
                                    shifted_new_perm = ~isnp
                                    #wordo.reverse()
                                    build_perm = ~Permutation.ref_product(*(wordo[:ii]+wordo[ii+1:]))
                        else:
                            build_perm = build_perm.swap(shifted_new_perm[0] - 1, shifted_new_perm[bp - 1] - 1)
                print(f"new {shifted_new_perm=}")
                
                # if build_perm.inv != up_perm.inv - new_perm.inv:
                #     print(f"{build_perm=}")
                #     raise Exception        
                cyc = ((~perm0) * up_perm)

                # the_cycles = cyc.get_cycles()
                # print(f"{the_cycles}")
                # cyc2 = ~cyc
                # reflecs = []
                # from schubmult.utils.perm_utils import cyclic_sort, has_bruhat_ascent
                # perm1 = perm0
                # last = 0
                # while perm1 != up_perm:
                #     print(f"{perm1=} {up_perm} {cyc2=}")
                #     for aa in range(the_index - 1, -1, -1):
                #         for bb in range(the_index,len(up_perm)):
                #             if has_bruhat_ascent(perm1, aa, bb):
                #                 last = perm1[bb]
                #                 cyc2 = cyc2.swap(aa, bb)
                #                 perm1 = perm1.swap(aa, bb)
                # # for cycle in the_cycles:
                # #     cycle = cyclic_sort([*cycle])
                # #     print(f"{cycle=}")
                # #     spot = min(cycle)
                # #     cycle.remove(spot)
                # #     cycle = [cycle[-1]] + cycle[:-1]
                # #     #cycle.reverse()
                # #     for c in cycle:
                # #         reflecs += [(spot, c)]
                # print(reflecs)
                # # while up_perm != perm0:
                # #     for a0 in range(1, the_index)
                # #reflecs.reverse()
                # spindle = Permutation([])
                # for a0, b0 in reflecs:
                #     spindle = spindle.swap(a0-1, b0-1)
                #     if shifted_new_perm[a0-1] < shifted_new_perm[b0-1]:
                #         if build_perm[shifted_new_perm[a0-1]-1] < build_perm[shifted_new_perm[b0 - 1]-1]:
                #             print("oh no")
                #             print(f"{build_perm=} {a0=} {b0=}")
                #         build_perm = build_perm.swap(shifted_new_perm[a0 - 1]-1, shifted_new_perm[b0 - 1] - 1)
                #     else:
                #         shifted_new_perm = shifted_new_perm.swap(a0-1, b0-1)
                # print(f"{spindle=} {~cyc=}")
                #assert spindle == ~cyc



                # build_perm = full_ref_prod * shifted_new_perm
                # build_perm *= Permutation.ref_product(*list(range(1,the_index)))
                # build_perm = build_perm * (~shifted_new_perm)
                #build_perm = build_perm * ((~perm0) * build_perm)
                # print(f"piff {cyc=} {index_list_sort=} {perm0=} {new_perm=}")
                # if len(index_list) == ((~cyc) * perm0).inv:
                #     print(f"piff {cyc=} {index_list_sort=} {perm0=}")
                #     ref_add = list([i0 for i0 in index_list_sort])
                #     new_ref_stack.append(ref_add + ref0)#[*list(range(len(index_list) + i,0 + i,-1))])
                # else:
                # if i == 0 and perm0[0] != new_perm[0]:
                #     ref_add = list([i0 for i0 in index_list_sort])
                #     new_ref_stack.append(ref0 + ref_add)
                    
                # do real

                # bongperm = perm0 * (~Permutation([1] + [jj+1 for jj in new_perm]))

                #[*list(range(len(index_list) + i,0 + i,-1))])
                
                
                puggle = reduced_word(~cyc)
                pingus = cyc
                shiftedp = shifted_new_perm
                up_perm1 = up_perm
                from schubmult.utils.perm_utils import has_bruhat_descent
                last = 0
                
                while pingus.inv > 0:
                    cycles = pingus.get_cycles()
                    for c in cycles:
                        aa = c[-2]
                        bb = c[-1]
                        print(f"Bingo {the_index=} {aa=} {bb=} {up_perm1=} {pingus=} {perm0=}")
                        up_perm1 = up_perm1.swap(aa-1, bb-1)
                        pingus = pingus.swap(aa-1, bb-1)
                        if shiftedp[aa - 1] > shiftedp[bb - 1]:
                            shiftedp = shiftedp.swap(aa-1, bb-1)
                        elif build_perm[shiftedp[aa-1] - 1] > build_perm[shiftedp[bb-1] - 1]:
                            build_perm = build_perm.swap(shiftedp[aa-1] - 1, shiftedp[bb - 1] - 1)
                            
                            print(f"{reduced_word(build_perm)=}")
                            break
                        else:
                            print("Binglehoffy")
                            raise Exception
                        
                        
                print(f"{up_perm1=} {perm0=}")
                assert up_perm1 == perm0
                assert build_perm * shiftedp == perm0

                    
                    # print(f"{puggle=}")
                    # new_puggle = []
                    # for ai in range(len(puggle) - 1, -1, -1):
                    #     aa, bb = root_at(puggle, ai)
                    #     print(f"{aa,bb}")
                    #     if aa <= the_index and bb>the_index and has_bruhat_descent(up_perm1, aa - 1, bb - 1):# and up_perm1[bb - 1] > last:
                    #         last = up_perm1[bb - 1]
                    #         new_puggle = puggle[:ai] + new_puggle
                    #         break
                    #     else:
                    #         new_puggle = [puggle[ai]] + new_puggle
                    # if len(puggle) == len(new_puggle):
                    #     raise Exception
                    # puggle = reduced_word(Permutation.ref_product(*new_puggle))
                            

                # for po in puggle:
                #     if shiftedp[po-1] > shiftedp[po]:
                #         shiftedp = shiftedp.swap(po-1, po)
                #     else:
                #         #created_word += [po]
                #         build_perm = build_perm.swap(shiftedp[po-1] -1,shiftedp[po] -1)

                #build_perm *= (shifted_new_perm*(~cyc)*(~shifted_new_perm))
                #build_perm *= (shifted_new_perm*(Permutation.ref_product(*created_word))*(~shifted_new_perm))
                # if perm0.inv == 0:
                #     build_perm = Permutation([])
                if ((~build_perm)*perm0).inv != perm0.inv - build_perm.inv:
                    raise Exception
                    
                word0 = reduced_word(build_perm)
                assert len(word0) == build_perm.inv
                assert Permutation.ref_product(*word0) == build_perm
                
                print(f"{word0=}")
                outperm = (~build_perm)*perm0
                outperm = Permutation([ap - 1 for ap in outperm if ap > 1])
                new_perm_stack.append(outperm)
                new_ref_stack.append(ref0 * shiftup(build_perm,i))
                new_seq_stack.append([a]*(len(index_list))+ seq0)
                new_bongo_stack.append(bongo0 + [ww + i for ww in word0])
        perm_stack = new_perm_stack
        ref_stack = new_ref_stack
        seq_stack = new_seq_stack
        bongo_stack = new_bongo_stack
    print(f"{perm.code=}")
    for ref, seq, new_perm, bongo in zip(ref_stack,seq_stack,perm_stack,bongo_stack):
        print(f"{new_perm=}")
        print(seq)
        print(bongo)
        print(ref)
        
        flob = Permutation([])
        # boing = NilPlactic.from_word([])
        for r in bongo:
            if flob[r-1] > flob[r]:
                print("Fail")
            flob = flob.swap(r-1, r)
        #     boing = boing.ed_insert(r)
        print(f"{flob=}")
        if flob != ref:
            print("Failed")
        # print(f"{boing=}")
        # #shape = [len(b) for b in boing._word]
        # # pony_perm = Permutation([])
        # # spuggle = list(boing._word)
        # # spuggle.reverse()
        # # for ind in range(len(spuggle[-1])):
        # #     for ind2 in range(len(spuggle)):
        # #         if ind < len(spuggle[ind2]):
        # #             pony_perm = pony_perm.swap(spuggle[ind2][ind] - 1, spuggle[ind2][ind])            
        bongset.add((ref,new_perm))
        bongdict[(ref,new_perm)] = bongdict.get((ref,new_perm),0) + 1
    print(f"{bongset=}")
    print(f"{songset=}")
    #print(f"{songset==bongset}")
    print(f"{bongset.issuperset(songset)=}")
    #print(f"{bongset.issuperset(songset)=}")
    print(f"{songset.difference(bongset)=}")
    print(f"{bongdict=}")
    print(f"{songdict=}")
        

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
