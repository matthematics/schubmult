from ._funcs import (
    mult_poly,
    schubmult,
)
from symengine import sympify
from schubmult._base_argparse import schub_argparse
from schubmult.perm_lib import (
    inverse,
    theta,
    permtrim,
    inv,
    mulperm,
    code,
    uncode,
    trimcode,
)


def main():
    try:
        args, formatter = schub_argparse("schubmult_py", "Compute products of ordinary Schubert polynomials")

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
        pr = args.pr
        coprod = args.coprod

        if coprod:
            if ascode:
                perms[0] = tuple(permtrim(uncode(perms[0])))
            pos = [*perms[1]]
            pos.sort()
            mperm = perms[0]

            cd = code(mperm)
            perms[0] = mperm

            while cd[-1] == 0:
                cd.pop()
            k = len(pos)
            n = len(perms[0])
            kcd = [pos[i] - i - 1 for i in range(len(pos))] + [n + 1 - k for i in range(k, n)]
            N = len(kcd)
            kperm = inverse(uncode(kcd))
            coeff_dict = {tuple(permtrim(kperm)): 1}
            coeff_dict = schubmult(coeff_dict, tuple(permtrim([*perms[0]])))

            inv_kperm = inv(kperm)
            inverse_kperm = inverse(kperm)
            if pr:
                for perm, val in coeff_dict.items():
                    downperm = mulperm(list(perm), inverse_kperm)
                    if inv(downperm) == inv(perm) - inv_kperm:
                        flag = True
                        for i in range(N):
                            if downperm[i] > N:
                                flag = False
                                break
                        if not flag:
                            continue
                        firstperm = downperm[0:N]
                        secondperm = [downperm[i] - N for i in range(N, len(downperm))]
                        if val != 0:
                            # firstcode = code(firstperm)
                            # while len(firstcode)>0 and firstcode[-1] == 0:
                            # firstcode.pop()
                            # secondcode = code(secondperm)
                            # while len(secondcode)>0 and secondcode[-1] == 0:
                            # secondcode.pop()
                            if ascode:
                                print(f"{val} {trimcode(firstperm)} {trimcode(secondperm)}")
                            else:
                                print(
                                    f"{val} {tuple(permtrim(firstperm))} {tuple(permtrim(secondperm))}"
                                )
        else:
            if ascode:
                for i in range(len(perms)):
                    perms[i] = tuple(permtrim(uncode(perms[i])))

            perms.sort(reverse=True, key=lambda x: sum(theta(inverse(x))) - inv(x))

            coeff_dict = {tuple(permtrim([*perms[0]])): 1}

            for perm in perms[1:]:
                coeff_dict = schubmult(coeff_dict, tuple(permtrim([*perm])))
            if mult:
                mul_exp = sympify(mulstring)
                coeff_dict = mult_poly(coeff_dict, mul_exp)

            if pr:
                for perm, val in coeff_dict.items():
                    if val != 0:
                        if ascode:
                            print(f"{val}  {trimcode(perm)}")
                        else:
                            print(f"{val}  {perm}")
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    main()
