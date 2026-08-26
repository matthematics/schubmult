from schubmult import *

# Let $a = a_1\cdots a_n$ be a reduced word for $w\in S_\infty$ such that $\maxd(w)=m$. Let $i_1>\cdots >i_k$ be the indexes such that, defining $w^0 = w$, we have that
# $$w^j = w^{j-1}s_{m - 1 + j}$$
# and
# $$a_1\cdots \widehat{a_{i_j}}\cdots\widehat{a_{i_1}}\cdots a_n$$
# is a reduced word for $w^j$. If $k>1$, let $t_{pq}$ be the right reflection at index $i_k$ in $a$ and let 
# $$\widehat{a} = a_1\cdots \widehat{a_{i_k}}\cdots a_n.$$
# That is, only index $i_k$ is deleted. Then under the hypothesis that every step in Little's bumping algorithm moves left in the word,
# and
# $$\mathrm{perm}\left(\ltb{a}{i_1}\right)=w',$$
# then we have that
# $$\mathrm{perm}\left(\ltb{\widehat{a}}{i_1}\right)=w't_{pq}$$
# and $\ell(w't_{pq}) = \ell(w') - 1$.

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        if perm.inv == 0:
            continue
        if perm.max_descent <= 1:
            continue
        for rc in [rcc for rcc in RCGraph.all_rc_graphs(perm) if len(rcc[-1]) == 0]:
            working_perm = rc.perm
            #max_root = [len(rc), len(rc) + 1]
            desc = len(rc)
            while working_perm[desc - 1] > working_perm[desc]:
                working_perm = working_perm.swap(desc - 1, desc)
                desc += 1
            desc -= 1
            if desc == len(rc):
                continue
            reflection = (len(rc), desc + 1)
            rc2 = rc.zero_out_last_row()
            coords = None
            index = None
            for i in range(rc.perm.inv):
                coords2 = rc.left_to_right_inversion_coords(i)
                if rc.right_root_at(*coords2) == reflection:
                    coords = coords2
                    index = i
                    break
            if coords is None:
                raise ValueError(f"Could not find reflection {reflection} in {rc.perm}")
            rc_down1 = rc.toggle_ref_at(*coords).zero_out_last_row()
            rc_down2 = rc2.toggle_ref_at(*rc2.left_to_right_inversion_coords(index))
            assert rc_down1 == rc_down2, f"RC graphs do not match: {rc_down1} != {rc_down2}\nOriginal:\n{rc}"