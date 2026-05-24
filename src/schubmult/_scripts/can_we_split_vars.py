from schubmult import *
from schubmult.utils.schub_lib import hw_elementary_tensors
from schubmult.rings.combinatorial.bounded_rc_factor_algebra import BoundedRCFactorAlgebra
from sympy import sstr
import itertools

def _tensor_to_rcs(weight_tensor, descents):
    rcs = []
    for i, desc in enumerate(descents):
        wt = len(weight_tensor[i])
        weight = [0] * desc
        for w in weight_tensor[i]:
            weight[w - 1] = 1
        elem_rc = next(iter(RCGraph.all_rc_graphs(uncode([0]*(desc - wt) + [1] * wt), desc, weight=tuple(weight))))
        rcs.append(elem_rc)
    return CrystalGraphTensor(*rcs)

def all_tensors(weights, descents):
    hw_tensors = hw_elementary_tensors(weights, descents)
    for hw_tensor_weight in hw_tensors:
        hw_tensor = _tensor_to_rcs(hw_tensor_weight, descents)
        for tensor in hw_tensor.full_crystal:
            yield tensor


def split_on_indices(perm, indices):
    from schubmult.rings.polynomial_algebra import Schub, MonomialBasis
    from schubmult.utils.schub_lib import mu_A
    n = len(perm)
    indices = sorted(indices)
    r = BoundedRCFactorAlgebra()
    if indices[-1] >= n:
        raise ValueError("Indices must be less than the length of the permutation")
    descents_in = tuple(reversed([n - i for i in indices]))
    descents_out = tuple([i for i in range(1, n) if i not in descents_in])
    perm_coefficient = perm
    w0 = Permutation.w0(n)
    poly = Schub(perm, n - 1).change_basis(MonomialBasis)
    coeff_dict = {}
    rcs_seen = set()
    mu1 = uncode(mu_A(Permutation.w0(n).trimcode, [j - 1 for j in indices]))
    mu2 = uncode(mu_A(Permutation.w0(n).trimcode, [j - 1  for j in range(1, n) if j not in indices]))
    # mu1 = uncode(mu_A(Permutation.w0(n).trimcode, [j - 1 for j in descents_in]))
    # mu2 = uncode(mu_A(Permutation.w0(n).trimcode, [j - 1 for j in descents_out]))
    for perm_rc in RCGraph.all_rc_graphs(perm, n - 1):
        
        weights1 = perm_rc.length_vector
        #weights = tuple(reversed([n - 1 - j - w for j, w in enumerate(weights1)]))
        comp_weight = [0] * (n - 1)
        for row in range(n - 1):
            for col in range(n - 1 - row):
                if not perm_rc.has_element(row + 1, col + 1):
                    comp_weight[col] += 1
        weights = tuple(reversed(comp_weight))
        # weights = tuple(reversed([n - 1 - j - w for j, w in enumerate((perm_coefficient).pad_code(n - 1))]))
        #weights = perm_coefficient.pad_code(n - 1)
        # print(descents_in, descents_out)
        # print(n)
        # mu1 = ~uncode(tuple(reversed(descents_in)))
        # mu2 = ~uncode(tuple(reversed(descents_out)))
        # weights_in = tuple([weights[i - 1] if i in descents_in else 0 for i in range(1, n)])
        # weights_out = tuple([weights[i - 1] if i in descents_out else 0 for i in range(1, n)])
        # coeff_dict = {}
        # prin_rc = RCGraph.principal_rc(perm * Permutation.w0(n), n - 1)
        # prin_rc_hw, seq = prin_rc.to_highest_weight()
        # descents = sorted(descents_in + descents_out)
        #rcs_seen = set()    
        tlist1 = list(all_tensors([weights[d - 1] for d in descents_in], descents_in))
        tlist2 = list(all_tensors([weights[d - 1] for d in descents_out], descents_out))
        
        for t1, t2 in itertools.product(tlist1, tlist2):
            interleaved = [None] * (n - 1)
            for i, d in enumerate(descents_in):
                interleaved[d - 1] = t1.factors[i]
            for i, d in enumerate(descents_out):
                interleaved[d - 1] = t2.factors[i]
            interleaved_tensor = CrystalGraphTensor(*interleaved)
            rc = r.key_to_rc_graph(r.make_key(tuple(interleaved_tensor), n))
            if rc.perm == perm*w0 and rc.is_principal:
                # if not rcs_seen:
                #     rcs_seen.add(rc)
                # if rc not in rcs_seen:
                #     continue
                # if rc in rcs_seen:
                #     continue
                # rcs_seen.add(rc)
                # rc0 = rc
                # #interleaved_tensor = interleaved_tensor.reverse_raise_seq(seq)
                # for ttt in interleaved_tensor.full_crystal:
                #     #ttt = interleaved_tensor
                #     t1 = tuple([ttt.factors[i - 1] for i in descents_in])
                #     t2 = tuple([ttt.factors[i - 1] for i in descents_out])
                    
                #     # t1 = _tensor_to_rcs(tensor_in, descents_in)
                #     # t2 = _tensor_to_rcs(tensor_out, descents_out)

                
                
                rc1 = r.key_to_rc_graph(r.make_key(t1, n - 1))
                rc2 = r.key_to_rc_graph(r.make_key(t2, n - 1))
                if (rc1, rc2) in rcs_seen:
                    continue
                rcs_seen.add((rc1, rc2))
                # tensor_weight = tuple([rcc.perm.inv for rcc in interleaved_tensor])
                firstperm, secondperm = rc1.perm * ~mu1, rc2.perm * ~mu2
                if firstperm.inv != mu1.inv - rc1.perm.inv or secondperm.inv != mu2.inv - rc2.perm.inv:
                    continue
                coeff_dict[(rc1, rc2, firstperm, secondperm)] = coeff_dict.get((rc1, rc2, firstperm, secondperm), 0) + 1
    return coeff_dict

if __name__ == "__main__":
    from schubmult.utils.argparse import schub_argparse
    # n = int(sys.argv[1])
    # perms = Permutation.all_permutations(n)
    # for perm in perms:
    #     for indices in itertools.combinations(range(1, n), n // 2):
    #         print(f"Permutation: {perm}, Indices: {indices}")
    #         split_on_indices(perm, indices)
    import sys
    args, formatter = schub_argparse(
            "schubmult_py",
            "Compute products of ordinary Schubert polynomials",
            argv=sys.argv[1:],
        )
    print("Arguments:", args)
    mult = args.mult
    mulstring = args.mulstring

    perms = args.perms

    for perm in perms:
        try:
            for i in range(len(perm)):
                perm[i] = int(perm[i])
        except Exception as e:
            print("Permutations must have integer values")
            raise e

    ascode = args.ascode
    Permutation.print_as_code = ascode
    pr = args.pr
    coprod = args.coprod
    raw_result_dict = {}
    if ascode:
        perms[0] = uncode(perms[0])
    else:
        perms[0] = Permutation(perms[0])
    pos = [*perms[1]]
    pos.sort()
    coeff_dict = split_on_indices(perms[0], pos)

    for (rc1, rc2, firstperm, secondperm), val in coeff_dict.items():
        #val = coeff_dict[(firstperm, secondperm)]
        if val != 0:
            print(f"{val} {sstr(firstperm)} {sstr(secondperm)}")