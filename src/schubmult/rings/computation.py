def actual_code_to_compute_result(perm, length):
    from schubmult import schubmult_py
    from schubmult.rings.free_algebra import FA
    from schubmult.rings.rc_graph_module import RCGraph

    # special handling of trivial initial case
    if perm.inv == 0:
        if length == 0:
            return [(RCGraph(), RCGraph())]
        return [(RCGraph(() * length), RCGraph(() * length))]

    # generate the objects to check. This should not be included in the time complexity
    rc_set = set((FA(*perm.trimcode, *((0,) * (length - len(perm.trimcode)))) * RCGraph()).value_dict.keys())
    consideration_set = {(k[0], k[1]): v for k, v in (FA(*perm.trimcode, *((0,) * (length - len(perm.trimcode)))).coproduct() * (RCGraph() @ RCGraph())).value_dict.items()}

    consider_dict = {}
    # this are the pairs a_1, a_2. Generating the ordered set. Also shouldn't be included in the time complexity
    for (rc1, rc2), v in consideration_set.items():
        consider_dict[(rc1.perm, rc2.perm)] = consider_dict.get((rc1.perm, rc2.perm), set())
        consider_dict[(rc1.perm, rc2.perm)].add((rc1, rc2))

    # index of secret element
    secret_val = {}
    # actual secret element
    secret_element = {}

    for perm1, perm2 in consider_dict:
        for rc_graph in sorted(rc_set):
            if rc_graph.perm != perm:
                val = int(schubmult_py({perm1: S.One}, perm2).get(rc_graph.perm, 0))
                # sum up the index
                secret_val[(perm1, perm2)] = secret_val.get((perm1, perm2), 0) + val

    # grab the actual element at the index
    for k, st in consider_dict.items():
        lst = sorted(st)
        v = secret_val.get(k, 0)
        try:
            secret_element[k] = lst[v]
        except IndexError:
            pass

    ret_elem = []
    for k, v in consider_dict.items():
        if k in secret_element:
            for v2 in v:
                # polynomial time verification step
                if secret_element[k] <= v2:
                    ret_elem.append(v2)

    return ret_elem
