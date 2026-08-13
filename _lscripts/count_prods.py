from schubmult import *
from schubmult.symbolic import expand
from sympy import prod
import math

def _nonzero_filter(u, v1, perms):
    dom = v1.mul_dominant()
    return {k for k, v in (DSx(u) * DSx(dom, "z")).items() if expand(v) != 0 and k.inv - u.inv <= v1.inv}

def check_with_filter(u, v, perms):
    prd = DSx(u) * DSx(v, "z")
    actual_nonzero =  {k for k, v in prd.items() if expand(v) != 0}
    if actual_nonzero != _nonzero_filter(u, v, perms):
        return False
    return True

class TreeNode:
    def __init__(self, val=0, left=None, right=None):
        self.val = val
        self.left = left
        self.right = right
        
    def __repr__(self):
        return f"TreeNode({self.val},{self.left},{self.right})"
    
    def __eq__(self, other):
        if not isinstance(other, TreeNode):
            return False
        def _valeq(val1, val2):
            return (val1 is None and val2 is None) or val1 == val2
        return _valeq(self.val, other.val) and _valeq(self.left, other.left) and _valeq(self.right, other.right)
    
    def total_size(self):
        if self.left is None and self.right is None:
            return 1
        left_size = self.left.total_size() if self.left else 0
        right_size = self.right.total_size() if self.right else 0
        return 1 + left_size + right_size
    
    def extension_count(self):
        import math
        if self.left is None and self.right is None:
            return 1
        result = 1
        if self.left:
            result *= self.left.extension_count()
        if self.right:
            result *= self.right.extension_count()
        result *= math.comb(self.total_size() - 1, self.left.total_size() if self.left else 0)
        return result

    def __hash__(self):
        return hash((self.val, self.left, self.right))

def construct_tree_from_perm(iperm, n):
    """
    Constructs a binary tree from a 132-avoiding permutation.
    """
    # if len(perm) == 0:
    #     return None  # Base case: empty permutation corresponds to an empty tree
    
    # # Find the maximum element and its index
    # if isinstance(perm, Permutation) and n is None:
    #     n = len(perm)
    # if n is not None and isinstance(perm, Permutation):
    #     max_idx = perm[n - 1]
    # else:
    #     #max_val = max(perm)
    #     max_idx = perm[-1]
    # #max_idx = perm.index(max_val)
    
    # # The maximum element becomes the root
    # root = TreeNode(perm[max_idx])
    
    # # Elements to the left of the max form the left subtree
    # root.left = construct_tree_from_perm(perm[:max_idx])
    
    # # Elements to the right of the max form the right subtree
    # root.right = construct_tree_from_perm(perm[max_idx + 1:])
    
    # return root
    perm = ~iperm
    root = TreeNode(perm[n-1])
    for i in range(n - 2, -1, -1):
        current = root
        while True:
            if perm[i] < current.val:
                if current.left is None:
                    current.left = TreeNode(perm[i])
                    break
                else:
                    current = current.left
            else:
                if current.right is None:
                    current.right = TreeNode(perm[i])
                    break
                else:
                    current = current.right
    return root

def dominant_interval_size(dom_perm):
    if dom_perm.inv == 0:
        return 1
    # if not dom_perm.is_dominant:
    #     raise ValueError(f"{dom_perm.trimcode} not dominant")
    dpn = ~dom_perm
    dpc = dpn.code
    divv = prod([dpn[i] - dpc[i] if i < len(dpc) else dpn[i] for i in range(len(dpn))])
    return math.factorial(len(dom_perm))//divv

def bruh_size(perm):
    from schubmult.utils.perm_utils import has_bruhat_descent
    stack = [perm]
    seen = set()
    while len(stack) > 0:
        new_perm = stack.pop()
        if new_perm in seen:
            continue
        seen.add(new_perm)
        for i in range(len(new_perm) - 1):
            for j in range(i + 1, len(new_perm)):
                if has_bruhat_descent(new_perm, i, j):
                    stack.append(new_perm.swap(i, j))
    return len(seen)


def potato_sum(n):
    perms = [pp for pp in Permutation.all_permutations(n) if not pp.has_pattern([3,1,2])]
    bs = 0
    for perm in perms:
        md = perm.mul_dominant()
        da_perm = perm * (~md)
        bs = max(bs, bruh_size(da_perm))
    print(f"{bs} pants size vs {len(perms)}")

def actual_interval_size(perm1, *perm2):
    #@perm = perm2
    seen = set()
    iset = perm1.inversion_set
    for perm in perm2:
        stack = [perm]
        
        
        while len(stack) > 0:
            new_perm = stack.pop()
            if new_perm in seen:
                continue
            if not iset.issubset(new_perm.inversion_set):
                continue
            seen.add(new_perm)
            for i in (~new_perm).descents():
                stack.append(~((~new_perm).swap(i, i + 1)))
    return len(seen)

def actual_bruhat_interval_size(perm1, *perm2, descents=None):
    #perm = perm2
    
    seen = set()
    #iset = perm1.inversion_set
    for perm in perm2:
        stack = [perm]
        while len(stack) > 0:
            new_perm = stack.pop()
            if new_perm in seen:
                continue
            if not perm1.bruhat_leq(new_perm):#iset.issubset(new_perm.inversion_set):
                continue
            if descents is None or new_perm.descents().issubset(descents):
                seen.add(new_perm)
            for i in range(len(new_perm) - 1):#(~new_perm).descents():
                for j in range(i + 1, len(new_perm)):
                    if new_perm[i] > new_perm[j]:
                        stack.append(new_perm.swap(i, j))
    return len(seen)


def actual_interval(perm):
    stack = [perm]
    seen = set()
    while len(stack) > 0:
        new_perm = stack.pop()
        if new_perm in seen:
            continue
        seen.add(new_perm)
        for i in (~new_perm).descents():
            stack.append(~((~new_perm).swap(i, i + 1)))
    return seen


def count_nonzero(prd):
    return len({k for k,v in prd.items() if expand(v) != 0})

def count_nonzero_interval(perm1, perm2):
    return len({k for k,v in prd.items() if expand(v) != 0})


def elem(p, k, vv):
    return DSx(uncode([0] * (k - p) + [1] * p), vv)

#def is_irreducible(perm):

# count_nonzero <= z(u) + us_k

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    perms2 = [pp for pp in perms if pp.is_dominant]
    perms_by_descent = {k: {pp: pp for pp in perms if max(pp.descents(), default=-1) + 1 == k} for k in range(1, n)}
    mx_factor = 1

    # (2n - 1)!!
    # n!
    # (2n)!
    # 2^n n!
    for perm1 in perms:
        # if perm1.inv == 0:
        #     continue
        #d = perm1.max_descent
        for perm2 in perms:
            # if perm2.inv == 0:
            #     continue
            # v = uncode([a + 1 for a in perm2.trimcode[:-1]] + [1])
            # prd = DSx(perm1)*DSx(perm2.mul_dominant(),"z")
            # prd2 = DSx(perm1)*DSx(v,"z")
            # count_total = count_nonzero(prd)
            # count_total2 = (n**len((~(perm2.mul_dominant())).trimcode)-1)*count_nonzero(prd2)
            # assert count_total <= count_total2, f"n={n} p1={perm1.trimcode} p2={perm2.trimcode} count_total={count_total} count_total2={count_total2}"

            assert check_with_filter(perm1, perm2, perms)
            print("Stinkbat")

            # #print(f"{perm1,perm2}")
            # count_total2 = count_nonzero(DSx(perm1)*DSx(perm2,"z"))
            # perm_stinkbat = perm2#Permutation.w0(n)
            # ups = [k for k,v in (Sx(perm1) * Sx(perm_stinkbat.mul_dominant())).items() if v != 0]
            # ups2 = [k for k,v in (Sx(perm1) * Sx(perm_stinkbat)).items() if v != 0 and perm1.inversion_set.issubset(k.inversion_set)]
            # count_weak = actual_interval_size(perm1, *ups2)
            # count_bruh = actual_bruhat_interval_size(perm1, *ups)#, descents=perm2.descents() | perm1.descents())
            # #factor = dominant_interval_size(perm2.mul_dominant())
            # #assert count_bruh >= count_total
            # factor = count_bruh/count_weak
            # #assert count_weak <= count_total2
            # assert  factor * count_total2 - count_total > -1e-9, f"n={n} p1={perm1.trimcode} p2={perm2.trimcode} count_total={count_total} count_total2={count_total2} {factor=}"
            # mx_factor = max(mx_factor, factor)
            # print(f"{mx_factor=}")# {count_weak=} {count_bruh=}")
            # mu = perm2.mul_dominant()
            # mu0 = perm2.mul_sortable()
            # c = count_nonzero(DSx(perm1) * DSx(mu, "z"))/count_nonzero(DSx(perm1) * DSx(mu0, "z"))
            # if c >= mx:
            #     mx = c
            #     print(f"n={n}  max count={mx}  p1={perm1.trimcode}  p2={perm2.trimcode} {perm2.mul_dominant().trimcode}  {perm2.mul_sortable().trimcode}")
    #print(mx)
    # for k in range(1, n):
    #     for perm in perms:
    #         initial_max = count_nonzero(DSx(perm) * elem(1, k, "z"))
    #         current_max = initial_max
    #         current_inv = 1
    #         for perm2 in sorted(perms_by_descent[k].values(), key=lambda p: p.inv):
    #             if perm2.inv > current_inv:
    #                 current_inv = perm2.inv
    #                 current_max *= n
    #             prd = DSx(perm) * DSx(perm2, "z")
    #             assert count_nonzero(prd) <= current_max, f"n={n} k={k} p={p} perm={perm} count={count_nonzero(prd)} > {current_max}"
    #             #current_max *= n