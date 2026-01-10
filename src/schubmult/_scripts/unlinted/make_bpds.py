from schubmult import *
import numpy as np


def display_perms_vertical(perms):
    """
    Display a list of permutations vertically in a grid format.
    
    For columns after the first, print a * next to numbers that are equal 
    to the number in the same index in the previous column.
    
    Args:
        perms: List or tuple of Permutation objects
    """
    # header = "  ".join(f"Col{i:2d}".ljust(2) for i in range(len(perms)))
    # print(header)
    # print("-" * len(header))
    
    # Print each row (each index in the permutation)
    if not perms:
        return
    n = len(perms)

    grid = np.zeros((n - 1, n - 1), dtype=int)

    for col_idx in range(1, n):
        for row_idx in range(n - col_idx):
            value = perms[col_idx][row_idx]
            prev_value = perms[col_idx - 1][row_idx]
            
            if value == prev_value:
                #cell = (str(disp_value) + "*").ljust(2)
                grid[row_idx, col_idx - 1] = value + 1000  # Mark with offset for star
            else:
                grid[row_idx, col_idx - 1] = value
                        
    for row in grid:
        row_str = "  ".join(
            (str(val - 1000) + "*").ljust(2) if val > 1000 else "  " if val == 0 else str(val).ljust(2)
            for val in row
        )
        print(row_str)


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    w0 = Permutation.w0(n)
    for perm in perms:
        mx = len(perm)
        bpd_set = {((~perm)*Permutation.w0(mx),)}
        # bpd is a tuple of perms
        
        for index in range(1, mx):
            new_bpd_set = set()
            for bpd in bpd_set:
                last_perm = bpd[-1]
                down_perms = elem_sym_perms_op(last_perm, min(index,last_perm.inv), index)
                for down_perm, _ in down_perms:
                    if index == mx -1 and down_perm.inv != 0:
                        continue
                    new_bpd_set.add(bpd + (down_perm,))
            bpd_set = new_bpd_set
        print(f"BPDs for permutation {perm}:")
        for bpd in bpd_set:
            display_perms_vertical(tuple(reversed(bpd)))
            print()