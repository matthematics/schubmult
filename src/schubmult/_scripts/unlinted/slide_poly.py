from schubmult import *
from schubmult.symbolic import S, expand_seq


def left_shifted_weak_comps(comp):
    while len(comp) > 0 and comp[-1] == 0:
        comp = comp[:-1]

    length = len(comp)
    if length == 0:
        yield ()
        return

    nonzero_positions = [idx for idx, value in enumerate(comp) if value != 0]
    nonzero_values = [comp[idx] for idx in nonzero_positions]

    if not nonzero_positions:
        yield tuple(0 for _ in range(length))
        return

    chosen_positions = [0] * len(nonzero_positions)

    def backtrack(which, min_pos):
        if which == len(nonzero_positions):
            out = [0] * length
            for pos, value in zip(chosen_positions, nonzero_values):
                out[pos] = value
            yield tuple(out)
            return

        max_pos = nonzero_positions[which]
        for pos in range(min_pos, max_pos + 1):
            chosen_positions[which] = pos
            yield from backtrack(which + 1, pos + 1)

    yield from backtrack(0, 0)


def monomial_slide(comp, genset):
    ret = S.Zero
    for weak_comp in left_shifted_weak_comps(comp):
        ret += expand_seq(weak_comp, genset)
    return ret

def all_weak_comps(length, max_degree):
    if length == 0:
        yield ()
    elif length == 1:
        for i in range(max_degree + 1):
            yield (i,)
    else:
        for i in range(max_degree + 1):
            for tail in all_weak_comps(length - 1, max_degree - i):
                yield (i,) + tail

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    deg = int(sys.argv[2])
    for comp in all_weak_comps(n, deg):
        print(f"Comp: {comp}, Monomial slide: {monomial_slide(comp, Sx.genset)}")