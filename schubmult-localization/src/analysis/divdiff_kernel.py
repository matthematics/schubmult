from schubmult import RCGraph, RCGraphRing

def analyze_divdiff_kernel(n):
    """
    Analyze the kernel of divided differences for RC-graphs.
    
    This function identifies cases where divdiff(rc, i) = 0 for various RC-graphs
    and provides insights into the structure of the kernel.
    """
    perms = Permutation.all_permutations(n)
    rc_ring = RCGraphRing()
    
    kernel_cases = {
        'no_descent': [],
        'no_root': [],
        'not_hw': [],
    }
    
    for perm in perms:
        rc = RCGraph.principal_rc(perm, n-1)
        for i in range(1, n):
            diff = rc_ring(rc).divdiff(i)
            if len(diff) == 0:
                if i not in rc.perm.descents():
                    kernel_cases['no_descent'].append((rc, i))
                else:
                    desc_result = rc.divdiff_desc(i)
                    if desc_result is None:
                        kernel_cases['no_root'].append((rc, i))
                    else:
                        kernel_cases['not_hw'].append((rc, i))
    
    return kernel_cases


def summarize_kernel_cases(kernel_cases):
    """
    Summarize the findings from the kernel analysis.
    """
    print(f"Kernel Analysis Summary:")
    print(f"1. No descent cases: {len(kernel_cases['no_descent'])}")
    print(f"2. Descent but no root cases: {len(kernel_cases['no_root'])}")
    print(f"3. Cases where deletion does not yield highest weight: {len(kernel_cases['not_hw'])}")
    
    for case_type, cases in kernel_cases.items():
        print(f"\n{case_type.replace('_', ' ').title()} Cases:")
        for rc, i in cases:
            print(f"  RC {rc.perm.trimcode}: i={i}")
            print(f"  Weight: {rc.crystal_weight}")
            pretty_print(rc)


if __name__ == "__main__":
    import sys
    
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    kernel_cases = analyze_divdiff_kernel(n)
    summarize_kernel_cases(kernel_cases)