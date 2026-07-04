from functools import cached_property
from itertools import combinations

from .permutation import Permutation
from .rc_graph import RCGraph
from .wc_graph import WCGraph


def _is_compatible(compat_seq, word):
    """Check if compat_seq is compatible with word, i.e. for each prefix of word, the corresponding prefix of compat_seq has all distinct entries."""
    if len(compat_seq) != len(word):
        return False
    for i in range(1, len(word)):
        if compat_seq[i-1] > compat_seq[i]:
            return False
        if word[i - 1] <= word[i] and compat_seq[i-1] == compat_seq[i]:
            return False
    return True

class InversionsTableau:
    def __init__(self, _dict, *_, **__):
        self._dict = {k: frozenset(v) if isinstance(v, set) else v for k, v in _dict.items()}
        self._reverse_lookup = {}
        for k, v in self._dict.items():
            try:
                v = int(v)
                self._reverse_lookup[v] = self._reverse_lookup.get(v, set())
                self._reverse_lookup[v].add(k)
            except TypeError:
                if isinstance(v, frozenset):
                    for item in v:
                        try:
                            item = int(item)
                        except TypeError:
                            raise ValueError("Values must be integers or sets of integers")
                        self._reverse_lookup[item] = self._reverse_lookup.get(item, set())
                        self._reverse_lookup[item].add(k)
                else:
                    raise ValueError("Values must be integers or sets of integers")
        # print(self._reverse_lookup)

    def __getitem__(self, *key):
        if len(key) == 1:
            key = key[0]
        return self._dict[key]

    def __str__(self):
        return f"InversionsTableau({self._dict})"

    def __repr__(self):
        return f"InversionsTableau({self._dict!r})"

    def iter_keys(self, reverse=False):
        """Iterates in word order"""
        for lookup_key in sorted(self._reverse_lookup.keys(), reverse=reverse):
            yield from iter(sorted(self._reverse_lookup[lookup_key], key=lambda k: (k[0], -k[1]), reverse=reverse))

    def iter_items(self, reverse=False):
        """Iterates in word order"""
        for lookup_key in sorted(self._reverse_lookup.keys(), reverse=reverse):
            for key in sorted(self._reverse_lookup[lookup_key], key=lambda k: (k[0], -k[1]), reverse=reverse):
                yield (key, lookup_key)

    @classmethod
    def from_rc_graph(cls, rc):
        dct = {}
        for i in range(rc.perm.inv):
            dct[rc.left_to_right_inversion(i)] = rc.left_to_right_inversion_coords(i)[0]
        return cls(dct)

    @classmethod
    def from_wc_graph(cls, wc):
        word, seq = wc.perm_word, wc.compatible_sequence
        # wc.to_reduced_compatible_set_sequence()
        #perm = wc.perm
        dct = {}
        perm = Permutation([])
        bangle_word = []

        for i in range(len(word)):
            if perm[word[i] - 1] < perm[word[i]]:
                bangle_word = [*bangle_word, word[i]]
                dct = {Permutation.ref_product(word[i]).act_root(*k): v for k, v in dct.items()}
                root = (word[i], word[i] + 1)
                dct[root] = dct.get(root, set())
                dct[root].add(seq[i])
                perm = perm.swap(word[i] - 1, word[i])
            else:
                root = (bangle_word[-1], bangle_word[-1] + 1)
                dct[root].add(seq[i])
            # root = _abs(perm._right_root_at(i, word=word))
            # dct[root] = set(dct.get(root, set()))
            # dct[root].add(seq[i])
            # dct[root] = frozenset(dct[root])
        if Permutation.ref_product(*bangle_word) != wc.perm:
            raise ValueError("The word and compatible sequence do not correspond to the same permutation")
        return cls(dct)

    @cached_property
    #@property
    def perm_word(self):
        root_dict = {root: set(v) for root, v in self._dict.items()}
        ret_word = []
        while root_dict:
            max_pair = max([(root, v) for root, v in root_dict.items() if root[1] == root[0] + 1], key=lambda x: (max(x[1]), -x[0][1], x[0][0]))
            letter = max_pair[0][0]
            ret_word = [letter, *ret_word]
            root_dict[max_pair[0]].remove(max(max_pair[1]))
            if len(root_dict[max_pair[0]]) == 0:
                root_dict = {Permutation([]).swap(letter - 1, letter).act_root(*root): v for root, v in root_dict.items() if root != max_pair[0]}
        return tuple(ret_word)

    _sv_cache = {}  # noqa: RUF012

    def __eq__(self, other):
        return type(self) is type(other) and self._dict == other._dict

    @property
    def is_valid(self):
        # """Check if the tableau satisfies the necessary conditions to be an inversions tableau.
        # 1. If a root is a simple root, all values assigned to it must be less than or equal to the smaller index.
        # 2. If two roots have the same second index, their value sets must be disjoint.
        # 3. If there exist two roots that have indices (i, j) and (j, k), then necessarily there also exists (i,k) provided
        # i < j < k. In that case, for each such (i,k) and each v in (i,k) set there exist v1 in (i,j) set and v2 in (j,k) set such that either
        # v1 <= v < v2 or v2 < v <= v1.
        # """
        import itertools
        # try:
        #     if not self.to_wc_graph().is_valid:
        #         return False
        # except ValueError:
        #     return False
        # if InversionsTableau.from_wc_graph(self.to_wc_graph()) != self:
        #     return False
        # if False:
        vals = set(self._dict.keys())
        for r1, r2 in itertools.combinations(vals, 2):
            if r1[0] == r1[1] - 1:
                if not InversionsTableau.set_leq(self._dict[r1], {r1[0]}):
                    return False
            if r2[0] == r2[1] - 1:
                if not InversionsTableau.set_leq(self._dict[r2], {r2[0]}):
                    return False
            if r1[1] == r2[1]:
                if len(self._dict[r1].intersection(self._dict[r2])) > 0:
                    return False
                fr1, fr2 = r1, r2
                if fr1[0] < fr2[0]:
                    fr1, fr2 = fr2, fr1
                i, j, k = fr2[0], fr1[1], fr2[1]
                if (i, j) not in self._dict:
                    if not InversionsTableau.set_lt(self._dict[fr2], self._dict[fr1]):
                        return False
                else:
                    if not (InversionsTableau.set_lt(self._dict[fr1], self._dict[fr2]) and InversionsTableau.set_leq(self._dict[fr2], self._dict[(i, j)])) and not (InversionsTableau.set_leq(self._dict[(i,j)], self._dict[fr2]) and InversionsTableau.set_lt(self._dict[fr2], self._dict[fr1])):
                        return False
            if r1[0] == r2[0]:
                fr1, fr2 = r1, r2
                if fr1[1] > fr2[1]:
                    fr1, fr2 = fr2, fr1
                i, j, k = fr1[0], fr1[1], fr2[1]
                if (j, k) not in self._dict:
                    if not InversionsTableau.set_leq(self._dict[fr2], self._dict[fr1]):
                        return False
                else:
                    if not (InversionsTableau.set_leq(self._dict[fr1], self._dict[fr2]) and InversionsTableau.set_lt(self._dict[fr2], self._dict[(j, k)])) and not (InversionsTableau.set_lt(self._dict[(j,k)], self._dict[fr2]) and InversionsTableau.set_leq(self._dict[fr2], self._dict[fr1])):
                        return False
        return True

    @staticmethod
    def set_leq(set1, set2):
        """Check if set1 <= set2, i.e. max(set1) <= min(set2)."""
        if len(set1) == 0:
            return False
        if len(set2) == 0:
            return True
        return max(set1) <= min(set2)

    @staticmethod
    def set_lt(set1, set2):
        """Check if set1 < set2, i.e. max(set1) < min(set2)."""
        if len(set1) == 0:
            return False
        if len(set2) == 0:
            return True
        return max(set1) < min(set2)

    def _snap_min(self):
        return InversionsTableau({k: {min(v)} for k, v in self._dict.items()})

    def _snap_max(self):
        return InversionsTableau({k: {max(v)} for k, v in self._dict.items()})

    @classmethod
    def all_set_valued_inversions_tableaux(cls, perm, max_value=None):
        """Enumerate all set-valued inversions tableaux for ``perm``.

        This enumerates candidate root-label assignments directly from the
        defining axioms in :meth:`is_valid` (no WCGraph construction), then
        filters by validity and target permutation.

        The optional ``max_value`` bounds all labels by ``{1, ..., max_value}``.
        """
        perm = Permutation(perm)
        if max_value is None:
            max_value = perm.max_descent
        if max_value < 0:
            raise ValueError(f"max_value must be nonnegative, got {max_value}")

        key = (perm, int(max_value))
        if key in cls._sv_cache:
            return cls._sv_cache[key]

        ret = set()

        n = len(perm)
        if n <= 1:
            iv = cls({})
            if iv.perm == perm and iv.is_valid:
                ret.add(iv)
            cls._sv_cache[key] = ret
            return ret

        labels = tuple(range(1, max_value + 1))

        # Organize roots by second index. Axiom (2) is column-local and can be
        # enforced while building each column assignment.
        columns = []
        for j in range(2, n + 1):
            roots = [(i, j) for i in range(1, j)]
            columns.append(roots)

        def powerset(vals):
            out = [frozenset()]
            vals = list(vals)
            for r in range(1, len(vals) + 1):
                out.extend(frozenset(c) for c in combinations(vals, r))
            return out

        # Precompute allowed sets per root from axiom (1).
        allowed_for_root = {}
        for j in range(2, n + 1):
            for i in range(1, j):
                if j == i + 1:
                    allowed_labels = tuple(v for v in labels if v <= i)
                else:
                    allowed_labels = labels
                allowed_for_root[(i, j)] = powerset(allowed_labels)

        # Enumerate each column with pairwise-disjoint assigned sets (axiom 2).
        column_assignments = []
        for roots in columns:
            local = []

            def build_col(idx, used, partial):
                if idx == len(roots):
                    local.append(dict(partial))
                    return
                root = roots[idx]
                for st in allowed_for_root[root]:
                    if used.intersection(st):
                        continue
                    partial[root] = st
                    build_col(idx + 1, used.union(st), partial)
                partial.pop(root, None)

            build_col(0, set(), {})
            column_assignments.append(local)

        def backtrack_cols(cidx, current):
            if cidx == len(column_assignments):
                # Drop empty sets from the explicit dictionary representation.
                compact = {r: s for r, s in current.items() if len(s) > 0}
                iv = cls(compact)
                if iv.is_valid:
                    try:
                        if iv.perm == perm:
                            ret.add(iv)
                    except ValueError:
                        # Candidate root assignment may satisfy axioms checked by
                        # is_valid yet fail the greedy simple-root extraction used
                        # by perm_word; such candidates are discarded.
                        pass
                return

            for col_dct in column_assignments[cidx]:
                current.update(col_dct)
                backtrack_cols(cidx + 1, current)
                for r in col_dct:
                    current.pop(r, None)

        backtrack_cols(0, {})

        cls._sv_cache[key] = ret
        return ret

    @cached_property
    def perm(self):
        result = Permutation([])
        for swap in self.perm_word:
            result = result @ Permutation.ref_product(swap)
        return result

    @cached_property
    def compatible_sequence(self):
        result = []
        for key in sorted(self._reverse_lookup.keys()):
            result.extend([key] * len(self._reverse_lookup[key]))
        return tuple(result)

    @cached_property
    def is_reduced(self):
        return self.perm.inv == len(self.perm_word)

    def __hash__(self):
        return hash(tuple(sorted(self._dict.items())))

    @cached_property
    def reduced_word(self):
        word = []
        seq = self.compatible_sequence
        set_seq = []
        working_perm = Permutation([])
        for i, letter in enumerate(self.perm_word):
            if working_perm[letter - 1] > working_perm[letter]:
                for root_index in range(len(word)):
                    root = working_perm.right_root_at(root_index, word=word)
                    if root == (letter, letter + 1):
                        set_seq[root_index].add(seq[i])
                        break
            else:
                working_perm = working_perm.swap(letter - 1, letter)
                word.append(letter)
                set_seq.append({seq[i]})
        return tuple(word)#, tuple(tuple(sorted(s)) for s in set_seq)
        # working_perm = Permutation([])
        # reduced_word = []
        # for a in self.perm_word:
        #     if working_perm[a - 1] < working_perm[a]:
        #         reduced_word.append(a)
        #         working_perm = working_perm.swap(a - 1, a)
        # return tuple(reduced_word)

    @property
    def is_set_valued(self):
        return any(isinstance(v, set | frozenset) for v in self._dict.values())

    def to_rc_graph(self, length=None):
        if self.is_set_valued:
            raise ValueError("Inversions tableau must be reduced to convert to RC graph")
        return RCGraph.from_reduced_compatible(self.perm_word, self.compatible_sequence, length=length)

    def to_wc_graph(self, length=None):
        return WCGraph._from_root_dict({root: set(v) for root, v in self._dict.items()}, length=length)

    def polyvalue(self, x, y=None, *, beta=None, prop_beta=False):
        from schubmult import Gx
        from schubmult.symbolic import S

        if y is not None:
            raise NotImplementedError("This is not implemented yet")
        if beta is None:
            beta = Gx._beta
        result = S.One
        overage = len(self.perm_word) - self.perm.inv
        for v, st in self._reverse_lookup.items():
            result *= x[v] ** len(st)

        if prop_beta:
            result *= beta**overage
        else:
            result *= beta ** len(self.perm_word)
        # print(f"{list(self.iter_keys())=}")
        # print(f"{self.perm_word=}")
        return result

# if __name__ == "__main__":
#     import sys
#     import itertools
#     n = int(sys.argv[1])
#     perms = Permutation.all_permutations(n)
#     for perm in perms:
#         ivs = InversionsTableau.all_set_valued_inversions_tableaux(perm, max_value=perm.max_descent)
#         print(f"PERMUTATION: {perm}  (inv={perm.inv}, descents={perm.descents()}, trimcode={perm.trimcode})")
#         print(f"  Found {len(ivs)} set-valued inversions tableaux:")
#         for iv in sorted(ivs, key=lambda x: (x.perm_word, x.compatible_sequence)):
#             print(f"    {iv._dict}")
