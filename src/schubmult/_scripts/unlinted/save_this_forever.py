for seq in aseqs:
        poly = Sx(expand_seq(seq, x))
        for perm0, coeff0 in poly.items():
            if len(perm0) > n:
                continue
            #solution_module3 += coeff0 * TensorModule.ext_multiply(ASx(perm0, n-1) * unit_rc_module,FA(*seq).coproduct() * unit_tensor_rc_module)
            for rc in rc_graphs[perm0]:
                rc_module = RCGraphModule(dict.fromkeys(rc_graphs_by_weight[perm0].get(rc.length_vector(), set()),1))
                solution_module3 += coeff0 * TensorModule.ext_multiply(rc_module,FA(*seq).coproduct() * unit_tensor_rc_module)
