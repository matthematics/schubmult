from functools import cached_property

from .permutation import Permutation
from .rc_graph import RCGraph
from .wc_graph import WCGraph


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
        word, seq = wc.to_reduced_compatible_set_sequence()
        perm = wc.perm
        dct = {}
        for i in range(len(word)):
            dct[perm.right_root_at(i, word=word)] = frozenset(seq[i])
        return cls(dct)

    @cached_property
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
        try:
            self.to_wc_graph()
        except ValueError:
            # import traceback
            # traceback.print_exc()
            # print("Invalid inversion tableau: ", self)
            return False
        return True
        # import itertools
        # if self.perm.inv <= 1:
        #     return True

        # if set(self._dict.keys()) != self.perm.inversion_set:
        #     return False
        # if sum(len(v) for v in self._dict.values()) != len(self.perm_word):
        #     return False
        # #bloated_inversion_set = set()
        # smaller_iv_list = [InversionsTableau({k: min(v) for k, v in self._dict.items()}), InversionsTableau({k: max(v) for k, v in self._dict.items()})]
        # for smaller_iv in smaller_iv_list:
        #     for (root1, val1), (root2, val2) in itertools.combinations(smaller_iv._dict.items(), 2):
        #         if root1[1] == root1[0] + 1 and val1 > root1[0]:
        #             return False
        #         if root2[1] == root2[0] + 1 and val2 > root2[0]:
        #             return False

        #         # if root1 != root2:
        #         #     if root1[1] == root2[1] and val1 == val2:
        #         #         return False
        #         if root1[1] == root2[0]:
        #             root3 = (root1[0], root2[1])
        #             if not (val1 <= smaller_iv[root3] < val2 or val2 < smaller_iv[root3] <= val1):
        #                 return False
        #         #     # good = True
        #         #     for (r3, val3) in bloated_inversion_set:
        #         #         if r3 == root3:
        #         #             if val1 == val3 and val2 > val3:
        #         #                return False
        #                     # if val3 > val1 and val3 > val2:
        #                     #     return False
        #             #     if r3 == root3:
        #             #         if val1 == val3:
        #             #             continue
        #             #         if val2 < val3:
        #             #             if not any(val1_1 >= val3 for r1, val1_1 in bloated_inversion_set if r1 == root1):
        #             #                 good = False
        #             #                 break
        #             #         if val2 > val3:
        #             #             if not any(val1_1 <= val3 for r1, val1_1 in bloated_inversion_set if r1 == root1):
        #             #                 good = False
        #             #                 break
        #             #         if val2 == val3:
        #             #             good = False
        #             #             break
        #             # if not good:
        #             #     return False
        # return True

    @classmethod
    def all_set_valued_inversions_tableaux(cls, perm, max_value=None):
        # import itertools

        # ret = set()
        # if perm.inv == 0:
        #     return {cls({})}
        # if max_value is None:
        #     max_value = len(perm.trimcode)
        # if perm in cls._sv_cache:
        #     return {iv for iv in cls._sv_cache[perm] if max(iv._reverse_lookup.keys(), default=0) <= max_value}
        # for d in perm.descents(zero_indexed=False):
        #     down_perm = perm.swap(d - 1, d)
        #     cls.all_set_valued_inversions_tableaux(down_perm)  # load cache
        #     old_set = cls.all_set_valued_inversions_tableaux(down_perm)
        #     for old_iv in old_set:
        #         max_val = max(old_iv._reverse_lookup.keys(), default=0)
        #         old_perm_word = old_iv.perm_word

        #             new_dct = {**old_iv._dict}
        #             new_dct = {Permutation.ref_product(d).act_root(*key): set(valset) for key, valset in new_dct.items()}
        #             new_dct[(d, d + 1)] = {mxv}
        #             dd = [ddd for ddd in perm.descents(zero_indexed=False) if ddd != d]
        #                     new_iv = cls(new_dct)
        #                     if new_iv.is_valid:
        #                         # print("Warning invalid inversion tableau generated: ", new_iv)
        #                         ret.add(new_iv)

        # for r2 in range(1, max_value + 1 - mxv):
        #     for valset in itertools.combinations(list(range(mxv + 1, max_value + 1)), r2):
        #         new_dct2 = {**new_dct}
        #         new_dct2[(d, d + 1)] = set(valset) | {mxv}
        #         new_iv = cls(new_dct2)
        #         if not new_iv.is_valid:
        #             # print("Warning invalid inversion tableau generated: ", new_iv)
        #             continue
        #         ret.add(new_iv)
        #             new_dct = {}
        #             for key, valset in old_iv._dict.items():
        #                 new_dct[Permutation.ref_product(d).act_root(*key)] = valset

        #             first_new_dct = {**new_dct, (d, d + 1): {mxv}}
        #             new_iv = cls(first_new_dct)
        #             if new_iv.is_valid:
        #                 # print("Warning invalid inversion tableau generated: ", new_iv)
        #                 ret.add(new_iv)
        #             for r in range(1, max_value + 1 - mxv):
        #                 for valset in itertools.combinations(list(range(mxv + 1, max_value + 1)), r):
        #                     new_dct2 = {**new_dct}
        #                     new_dct2[(d, d + 1)] = set(valset) | {mxv}
        #                     new_iv = cls(new_dct2)
        #                     if not new_iv.is_valid:
        #                         # print("Warning invalid inversion tableau generated: ", new_iv)
        #                         continue
        #                     ret.add(new_iv)

        # if max_value == len(perm.trimcode):
        #     cls._sv_cache[perm] = ret
        # return ret
        raise NotImplementedError("This is not implemented yet")

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

    def to_rc_graph(self, length=None):
        if not self.is_reduced:
            raise ValueError("Inversions tableau must be reduced to convert to RC graph")
        return RCGraph.from_reduced_compatible(self.perm_word, self.compatible_sequence, length=length)

    def to_wc_graph(self, length=None):
        set_seq = [self[self.perm.right_root_at(i, word=self.reduced_word)] for i in range(len(self.reduced_word))]
        the_wc = WCGraph.from_reduced_compatible_set_sequence(self.reduced_word, set_seq, length=length)
        if the_wc.perm != self.perm:
            raise ValueError("Inversions tableau is not compatible with its own word")
        return the_wc

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
