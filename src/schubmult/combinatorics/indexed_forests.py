from dataclasses import dataclass


class Node:
    def __init__(self, label):
        self.label = label
        self.insert_label = None
        self.dec_label = None
        self.left = None
        self.right = None

    def __repr__(self):
        return f"Node({self.label}, L:{self.left}, R:{self.right})"

    @property
    def rho(self):
        return self.left.rho if self.left else self.label

    @property
    def inorder_traversal(self):
        if self.left:
            yield from self.left.inorder_traversal
        yield self
        if self.right:
            yield from self.right.inorder_traversal

    def __str__(self):
        return f"Node({self.label})"

    def __eq__(self, other):
        if not isinstance(other, Node):
            return False
        return self.label == other.label and self.left == other.left and self.right == other.right

    def __hash__(self):
        return hash((self.label, self.left, self.right))


class IndexedForest:
    def __init__(self, roots=(), code=None):
        self._roots = tuple(roots)
        self._code = None if code is None else tuple(int(x) for x in code)

    @property
    def roots(self):
        return self._roots

    @property
    def code(self):
        if self._code is not None:
            return self._code
        computed = self._compute_code()
        self._code = computed
        return self._code

    def _compute_code(self):
        if len(self._roots) == 0:
            return ()

        max_rho = 0
        for root in self._roots:
            for node in root.inorder_traversal:
                if node.rho > max_rho:
                    max_rho = node.rho

        if max_rho <= 0:
            return ()

        coeffs = [0] * max_rho

        for root in self._roots:
            for node in root.inorder_traversal:
                coeffs[node.rho - 1] += 1

        return tuple(coeffs)

    def __iter__(self):
        return iter(self._roots)

    def __len__(self):
        return len(self._roots)

    def __getitem__(self, index):
        return self._roots[index]

    def __repr__(self):
        return f"IndexedForest(roots={self._roots}, code={self._code})"

    def __eq__(self, other):
        if not isinstance(other, IndexedForest):
            return False
        return self._roots == other._roots

    def __hash__(self):
        return hash(self._roots)


def _node_support(root):
    return tuple(node.label for node in root.inorder_traversal)


def _forest_intervals(forest_roots):
    intervals = []
    for root in forest_roots:
        support = _node_support(root)
        if len(support) == 0:
            continue
        support_set = set(support)
        intervals.append((min(support_set), max(support_set), support_set, root))
    intervals.sort(key=lambda x: x[0])
    return intervals


@dataclass(frozen=True, order=True)
class ParallelInjLetter:
    primary: int
    secondary: int

    def __repr__(self):
        return f"{self.primary}[{self.secondary}]"


def make_parallel_injective_word(primary_word, secondary_word):
    if len(primary_word) != len(secondary_word):
        raise ValueError("primary_word and secondary_word must have same length")
    return tuple(ParallelInjLetter(int(a), int(b)) for a, b in zip(primary_word, secondary_word))


def _letter_val(letter, val_fn=None):
    if val_fn is not None:
        return int(val_fn(letter))
    if hasattr(letter, "val"):
        return int(letter.val)
    if isinstance(letter, int):
        return abs(letter)
    if isinstance(letter, str):
        token = letter.strip()
        if token.startswith("-") and token[1:].isdigit():
            return int(token[1:])
        if token.isdigit():
            return int(token)
    if isinstance(letter, ParallelInjLetter):
        return letter.primary
    if isinstance(letter, tuple | list) and len(letter) == 2:
        return int(letter[0])
    raise TypeError(f"Cannot extract val(letter) from {letter!r}")


def _letter_order_key(letter):
    if isinstance(letter, ParallelInjLetter):
        return (int(letter.primary), int(letter.secondary))
    if isinstance(letter, tuple | list) and len(letter) == 2:
        return (int(letter[0]), int(letter[1]))
    if isinstance(letter, int):
        return (int(letter), 0)
    if isinstance(letter, str):
        token = letter.strip()
        if token.startswith("-") and token[1:].isdigit():
            return (int(token[1:]), 0)
        if token.isdigit():
            return (int(token), 0)
    if hasattr(letter, "primary") and hasattr(letter, "secondary"):
        return (int(letter.primary), int(letter.secondary))
    return (_letter_val(letter), 0)


def insert_letter_into_indexed_forest(forest_roots, letter, step=None, val_fn=None, q_label=None):
    """
    Insert one letter into an indexed forest using the two-case insertion
    described in Nadeau--Tewari (Section 5.1).

    Parameters
    ----------
    forest_roots : iterable[Node]
        Current forest roots.
    letter : object
        Letter to insert. `val(letter)` is interpreted as `_letter_val(letter)`.
        Ordering comparisons use the LBS pair-letter order key `(a, b)`.
    step : int | None
        Optional insertion index for Dec-labeling.
    val_fn : callable | None
        Optional map to compute val(letter). If omitted, `_letter_val` defaults
        are used.
    q_label : object | None
        Optional explicit recording label for Q. If omitted, `step` is used.
    """
    roots = list(forest_roots)
    intervals = _forest_intervals(roots)
    support = set()
    for _, _, support_set, _ in intervals:
        support.update(support_set)

    i = _letter_val(letter, val_fn=val_fn)

    def interval_index_containing(value):
        for idx, (_, _, support_set, _) in enumerate(intervals):
            if value in support_set:
                return idx
        return None

    new_root = Node(i)
    new_root.insert_label = letter
    new_root.dec_label = step if q_label is None else q_label

    roots_to_remove = []

    if i not in support:
        left_idx = interval_index_containing(i - 1)
        right_idx = interval_index_containing(i + 1)

        if left_idx is not None:
            new_root.left = intervals[left_idx][3]
            roots_to_remove.append(intervals[left_idx][3])
        if right_idx is not None:
            new_root.right = intervals[right_idx][3]
            roots_to_remove.append(intervals[right_idx][3])
    else:
        j = interval_index_containing(i)
        if j is None:
            raise RuntimeError("Internal error: insertion value is in support but no interval found")

        _, _, _, root_j = intervals[j]
        root_j_label = root_j.insert_label if root_j.insert_label is not None else root_j.label

        if _letter_order_key(letter) > _letter_order_key(root_j_label):
            new_root.left = root_j
            roots_to_remove.append(root_j)

            right_neighbor_start = intervals[j][1] + 2
            right_neighbor_idx = interval_index_containing(right_neighbor_start)
            if right_neighbor_idx is not None:
                right_root = intervals[right_neighbor_idx][3]
                new_root.right = right_root
                roots_to_remove.append(right_root)
        else:
            new_root.right = root_j
            roots_to_remove.append(root_j)

            left_neighbor_end = intervals[j][0] - 2
            left_neighbor_idx = interval_index_containing(left_neighbor_end)
            if left_neighbor_idx is not None:
                left_root = intervals[left_neighbor_idx][3]
                new_root.left = left_root
                roots_to_remove.append(left_root)

    remaining = [r for r in roots if r not in roots_to_remove]
    remaining.append(new_root)
    return IndexedForest(remaining)


def word_to_indexed_forest(word, val_fn=None, q_labels=None):
    """
    Build an indexed forest from a word via repeated insertion.

    Returns
    -------
    IndexedForest
        Indexed forest with roots and optional code metadata.
    """
    forest = IndexedForest()
    if q_labels is None:
        q_iter = [None] * len(word)
    else:
        if len(q_labels) != len(word):
            raise ValueError("q_labels must have same length as word")
        q_iter = q_labels

    for step, (letter, q_label) in enumerate(zip(word, q_iter), start=1):
        forest = insert_letter_into_indexed_forest(forest, letter, step=step, val_fn=val_fn, q_label=q_label)
    return forest


def _clone_shape_with_labels(root, label_attr):
    if root is None:
        return None
    new_root = Node(root.label)
    chosen = getattr(root, label_attr)
    new_root.insert_label = chosen
    new_root.dec_label = None
    new_root.left = _clone_shape_with_labels(root.left, label_attr)
    new_root.right = _clone_shape_with_labels(root.right, label_attr)
    return new_root


def split_PQ_labelings(forest_roots):
    """
    Return two forests with the same shape as `forest_roots`:
    - P-forest: insertion labels
    - Q-forest: recording labels
    """
    p_forest = IndexedForest(tuple(_clone_shape_with_labels(root, "insert_label") for root in forest_roots), code=getattr(forest_roots, "code", None))
    q_forest = IndexedForest(tuple(_clone_shape_with_labels(root, "dec_label") for root in forest_roots), code=getattr(forest_roots, "code", None))
    return p_forest, q_forest


def parallel_word_to_indexed_forest(primary_word, secondary_word, val_part="primary"):
    word = make_parallel_injective_word(primary_word, secondary_word)
    if val_part == "primary":
        return word_to_indexed_forest(word, val_fn=lambda x: x.primary)
    if val_part == "secondary":
        return word_to_indexed_forest(word, val_fn=lambda x: x.secondary)
    raise ValueError("val_part must be 'primary' or 'secondary'")


def parallel_word_to_PQ_forests(primary_word, secondary_word, val_part="primary"):
    """
    Build one indexed-forest shape from a parallel injective word and return
    its P/Q labelings as two forests with identical shape.
    """
    word = make_parallel_injective_word(primary_word, secondary_word)
    if val_part == "primary":
        shape = word_to_indexed_forest(word, val_fn=lambda x: x.primary, q_labels=tuple(secondary_word))
    elif val_part == "secondary":
        shape = word_to_indexed_forest(word, val_fn=lambda x: x.secondary, q_labels=tuple(secondary_word))
    else:
        raise ValueError("val_part must be 'primary' or 'secondary'")
    return split_PQ_labelings(shape)


def weak_composition_to_indfor(c):
    """
    Converts a weak composition c+ into an indexed forest.
    Returns a list of root nodes for the trees in the forest.
    """
    c = tuple(int(x) for x in c)
    if any(x < 0 for x in c):
        raise ValueError("weak composition entries must be nonnegative")

    def build_tree_from_code_block(labels, code_block):
        if len(labels) == 0:
            return None
        if len(labels) != len(code_block):
            raise ValueError("labels and code_block must have the same length")

        m = len(labels)
        if sum(code_block) != m:
            raise ValueError(f"Invalid code block {tuple(code_block)} for labels {tuple(labels)}")

        running = 0
        root_pos = None
        for idx, val in enumerate(code_block, start=1):
            running += val
            if running < idx:
                raise ValueError(f"Invalid code block prefix at position {idx}: {tuple(code_block)}")
            if running == idx and root_pos is None:
                root_pos = idx

        if root_pos is None:
            raise ValueError(f"Could not determine root position from code block {tuple(code_block)}")

        root_label = labels[root_pos - 1]
        root = Node(root_label)

        left_labels = labels[: root_pos - 1]
        right_labels = labels[root_pos:]

        if root_pos == 1:
            if code_block[0] != 1:
                raise ValueError(f"Invalid leading code entry for block {tuple(code_block)}")
            left_code = []
        else:
            if code_block[0] < 1:
                raise ValueError(f"Invalid first code entry for block {tuple(code_block)}")
            left_code = [code_block[0] - 1, *code_block[1 : root_pos - 1]]

        right_code = code_block[root_pos:]

        root.left = build_tree_from_code_block(left_labels, left_code)
        root.right = build_tree_from_code_block(right_labels, right_code)
        return root

    forest_roots = []
    n = len(c)
    i = 0

    while i < n:
        while i < n and c[i] == 0:
            i += 1
        if i >= n:
            break

        prefix_sum = 0
        length = 0
        last_valid_j = None
        j = i

        while True:
            val = c[j] if j < n else 0
            prefix_sum += val
            length += 1

            if prefix_sum < length:
                break

            if prefix_sum == length:
                last_valid_j = j

            j += 1

        if last_valid_j is None:
            raise ValueError(f"Cannot parse composition block starting at index {i + 1} for {c}")

        block_len = last_valid_j - i + 1
        labels = list(range(i + 1, i + 1 + block_len))
        code_block = [c[k] if k < n else 0 for k in range(i, i + block_len)]
        forest_roots.append(build_tree_from_code_block(labels, code_block))

        i = last_valid_j + 1

    return IndexedForest(tuple(forest_roots), code=c)


def build_balanced_tree(labels):
    """
    Helper to build a tree where the in-order traversal matches the labels.
    This creates the 'canonical labeling' referenced in the paper.
    """
    if not labels:
        return None

    mid = len(labels) // 2
    root = Node(labels[mid])

    root.left = build_balanced_tree(labels[:mid])
    root.right = build_balanced_tree(labels[mid + 1 :])

    return root


def decreasing_labelings(root, max_val, used_vals=None):
    def _subtree_size(node):
        if node is None:
            return 0
        return 1 + _subtree_size(node.left) + _subtree_size(node.right)

    ret = set()
    used = set() if used_vals is None else set(used_vals)
    required_size = _subtree_size(root)
    available_labels = [val for val in range(1, max_val + 1) if val not in used]
    if len(available_labels) < required_size:
        return ret

    if root.left is None and root.right is None:
        return {(val,) for val in available_labels}

    left_size = _subtree_size(root.left)
    right_size = _subtree_size(root.right)

    for val in available_labels:
        remaining_smaller = sum(1 for a in available_labels if a < val)
        if remaining_smaller < left_size + right_size:
            continue

        vals_used = used | {val}
        if root.left is None:
            right_labelings = decreasing_labelings(root.right, val - 1, used_vals=vals_used)
            for right in right_labelings:
                ret.add((val, *right))
            continue

        if root.right is None:
            left_labelings = decreasing_labelings(root.left, val - 1, used_vals=vals_used)
            for left in left_labelings:
                ret.add((*left, val))
            continue

        right_labelings = decreasing_labelings(root.right, val - 1, used_vals=vals_used)
        for right in right_labelings:
            new_vals_used = vals_used | set(right)
            left_labelings = decreasing_labelings(root.left, val - 1, used_vals=new_vals_used)
            for left in left_labelings:
                ret.add((*left, val, *right))
    return ret
    # right_labelings = set()
    # left_labelings = set()

    # if root.left is not None:
    #     left_labelings.update(decreasing_labelings(root.left, val - 1))
    # if root.right is None:
    #     ret.update(tuple(list(left) + [val]) for left in left_labelings)
    #     return ret
    # if root.left is None:
    #     ret.update(tuple([val] + list(right)) for right in right_labelings)
    #     return ret
    # ret.update(tuple(list(left) + [val] + list(right)) for left in left_labelings for right in right_labelings)
    # return ret


def word_from_labeling(root, labeling):
    from schubmult import Permutation

    trav = list(root.inorder_traversal)
    assert len(trav) == len(labeling)
    perm = Permutation(labeling)
    # print(f"Labeling: {labeling}, Perm: {perm}, Rhos: {[node.rho for node in trav]}")
    return tuple(trav[(~perm)[i] - 1].rho for i in range(len(trav)))

def omega_insertion(word_of_pairs):
    if len(word_of_pairs) == 0:
        F = ()
        P, Q = (), ()
        return F, P, Q
    forest, P, Q = omega_insertion(word_of_pairs[:-1])
    return None




if __name__ == "__main__":
    # Example usage
    comp = (2, 0, 1)
    indfor = weak_composition_to_indfor(comp)
    print("composition:", comp)
    print(indfor)

    comp = (0, 2, 0, 1)
    indfor = weak_composition_to_indfor(comp)
    print("composition:", comp)
    print(indfor)

    comp = (0, 2, 0, 1, 0, 0, 1 ,0, 0, 0, 2)
    indfor = weak_composition_to_indfor(comp)
    print("composition:", comp)
    print(indfor)
