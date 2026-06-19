from functools import cached_property

from .permutation import Permutation
from .rc_graph import RCGraph
from .wc_graph import WCGraph


class InversionsTableau:
    def __init__(self, _dict, *_, **__):
        self._dict = {k: frozenset(v) if isinstance(v, set) else v  for k, v in _dict.items()}
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
        dct = {}
        for i in range(len(wc.perm_word)):
            key = wc.left_to_right_inversion(i)
            dct[key] = dct.get(key, set())
            dct[key].add(wc.left_to_right_inversion_coords(i)[0])
        return cls(dct)

    @cached_property
    def perm_word(self):
        keys = list(self.iter_keys())
        keys.reverse()
        build_perm = Permutation([])
        word = []
        for key in keys:
            the_root = (~build_perm).act_root(*key)
            if the_root[1] < the_root[0]:
                word.append(the_root[1])
            else:
                word.append(the_root[0])
            # else:
            if the_root[0] < the_root[1]:
                build_perm = build_perm.swap(the_root[0] - 1, the_root[1] - 1)

        return tuple(reversed(word))
        # if not keys:
        #     return ()
        # print("DEBUG: ", keys)
        # word_roots = []
        # word_roots = [keys[0]]
        # for key in keys[1:]:
        #     word_roots = [*word_roots[:-1], Permutation([]).swap(key[0] - 1, key[1] - 1).act_root(*word_roots[-1]), key]
        # print("DEBUG: ", word_roots)
        # raise NotImplementedError("This is not fully implemented yet, and may be buggy. The idea is to iteratively apply root swaps to the keys to build the permutation word, but it needs more work.")
        # while keys:
        #     # print("Debug: ", keys)
        #     key = keys.pop()
        #     if key[1] != key[0] + 1:
        #         raise ValueError("Invalid key: ", key)
        #     the_swap = Permutation.ref_product(key[0])
        #     # print("Debug: ", key, the_swap)
        #     new_keys = []
        #     swapping = True
        #     for k in reversed(keys):
        #         if k == key:
        #             swapping = not swapping
        #         if swapping:
        #             new_keys = [the_swap.act_root(*k), *new_keys]
        #         else:
        #             new_keys = [k, *new_keys]
        #     keys = new_keys
        #     #keys = [the_swap.act_root(*k) if k != key else k for k in keys]
        #     #keys = [the_swap.act_root(*k) for k in keys]
        #     # print("Debug: ", keys)
        #     result.append(key[0])
        #     # print("Debug: ", result)

        # return tuple(reversed(result))

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

    def to_rc_graph(self):
        if not self.is_reduced:
            raise ValueError("Inversions tableau must be reduced to convert to RC graph")
        return RCGraph.from_reduced_compatible(self.perm_word, self.compatible_sequence)

    def to_wc_graph(self):
        return WCGraph.from_word_compatible(self.perm_word, self.compatible_sequence)

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
            result *= (x[v] ** len(st))

        if prop_beta:
            result *= beta**overage
        else:
            result *= beta ** len(self.perm_word)
        # print(f"{list(self.iter_keys())=}")
        # print(f"{self.perm_word=}")
        return result
