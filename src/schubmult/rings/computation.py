def actual_code_to_compute_result(perm, length):
    from schubmult import schubmult_py
    from schubmult.rings.free_algebra import FA
    from schubmult.rings.rc_graph_module import RCGraph

    if perm.inv == 0:
        if length == 0:
            mod = [(RCGraph(), RCGraph())]
        else:
            mod = [(RCGraph(() * length), RCGraph(() * length))]
   
    rc_set = set((FA(*perm.trimcode, *((0,) * (length - len(perm.trimcode)))) * RCGraph()).value_dict.keys())
    consideration_set = {(k[0], k[1]): v for k, v in (FA(*perm.trimcode, *((0,) * (length - len(perm.trimcode)))).coproduct() * (RCGraph() @ RCGraph())).value_dict.items()}

    consider_dict = {}
    for (rc1, rc2), v in consideration_set.items():
        consider_dict[(rc1.perm, rc2.perm)] = consider_dict.get((rc1.perm, rc2.perm), set())
        consider_dict[(rc1.perm, rc2.perm)].add((rc1, rc2))

    ret_elem = None

    for perm1, perm2 in consider_dict:
        for rc_graph in sorted(rc_set):
            if rc_graph.perm != perm:
                val = int(schubmult_py({perm1: S.One}, perm2).get(rc_graph.perm, 0))
                lst = sorted(consider_dict[(perm1, perm2)])
                for i in range(val):
                    consider_dict[(perm1, perm2)].remove(lst[i])

    for k, v in consider_dict.items():
        if ret_elem is None:
            ret_elem = list(v)
        else:
            ret_elem += list(v)
    with lock:
        shared_cache_dict[perm] = ret_elem
    return ret_elem