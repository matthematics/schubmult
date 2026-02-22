class Node:
    def __init__(self, label):
        self.label = label
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


def weak_composition_to_indfor(c):
    """
    Converts a weak composition c into an indexed forest.
    Returns a list of root nodes for the trees in the forest.
    """
    forest_roots = []
    n = len(c)
    i = 0

    while i < n:
        while i < n and c[i] == 0:
            i += 1
        if i >= n:
            break

        segment_start = i + 1
        segment_size = 0
        zero_run = 0

        while i < n:
            part = c[i]
            if part == 0:
                zero_run += 1
                if zero_run >= 2:
                    break
            else:
                zero_run = 0
                segment_size += part
            i += 1

        if segment_size > 0:
            tree_support = list(range(segment_start, segment_start + segment_size))
            forest_roots.append(build_balanced_tree(tree_support))

        while i < n and c[i] == 0:
            i += 1

    return tuple(forest_roots)


def build_balanced_tree(labels):
    """
    Helper to build a tree where the in-order traversal matches the labels.
    This creates the 'canonical labeling' referenced in the paper.
    """
    if not labels:
        return None

    # Pick the middle element to keep the tree balanced (standard for binary search trees)
    mid = len(labels) // 2
    root = Node(labels[mid])

    root.left = build_balanced_tree(labels[:mid])
    root.right = build_balanced_tree(labels[mid + 1 :])

    return root
