import math
from functools import cache, cached_property

import sympy.combinatorics.permutations as spp

# import schubmult.schub_lib.schub_lib as sl
import schubmult.utils.logging as lg
from schubmult.utils.perm_utils import cyclic_sort, old_code, permtrim_list, sg

# schubmult.poly_lib.variables import GeneratingSet

logger = lg.get_logger(__name__)

zero = 0
n = 100

# TODO: permutations act


class Permutation:
    @classmethod
    def all_permutations(cls, n):
        from itertools import permutations

        return [cls(perm) for perm in list(permutations(list(range(1, n + 1))))]

    def parabolic_reduce(self, *descs):
        descs = set(descs)
        reduced_perm = self
        w_J = Permutation([])
        found = True
        while found:
            found = False
            for d in reduced_perm.descents():
                if d not in descs:
                    w_J = ~((~w_J).swap(d, d + 1))
                    reduced_perm = reduced_perm.swap(d, d + 1)
                    found = True
                    break
        return reduced_perm, w_J

    def __truediv__(self, other):
        """Returns a tuple of (self, other) if other is a permutation. Intended for skew elements."""
        if isinstance(other, Permutation):
            return (self, other)
        raise NotImplementedError("Division by non-permutation is not implemented.")

    def __new__(cls, perm):
        return Permutation.__xnew_cached__(cls, tuple(perm))

    # @classmethod
    # def _af_new(cls, p):
    #     obj = object.__new__(cls)
    #     p = tuple(p)
    #     obj._args = (p,)
    #     # obj._s_perm = tuple([i - 1 for i in p])
    #     obj._perm = p
    #     obj._hash_code = hash(p)
    #     cd = old_code(p)
    #     obj._unique_key = (len(p), sum([cd[i] * math.factorial(len(p) - 1 - i) for i in range(len(cd))]))
    #     return obj

    print_as_code = False

    @classmethod
    def ref_product(cls, *args):
        p = cls([])
        for a in args:
            p = p.swap(a - 1, a)
        return p

    @staticmethod
    @cache
    def __xnew_cached__(_class, perm):
        return Permutation.__xnew__(_class, perm)

    @staticmethod
    def __xnew__(_class, perm):
        p = tuple(permtrim_list([*perm]))
        # s_perm = spp.Permutation([i - 1 for i in p])
        obj = object.__new__(_class)
        obj._args = (p,)
        # obj._s_perm = tuple([i - 1 for i in p])
        obj._perm = p
        obj._hash_code = hash(p)
        cd = old_code(p)
        obj._unique_key = (len(p), sum([cd[i] * math.factorial(len(p) - 1 - i) for i in range(len(cd))]))
        return obj

    @property
    def args(self):
        return self._args

    @classmethod
    def sorting_perm(cls, itera, reverse=False):
        L = [i + 1 for i in range(len(itera))]
        L.sort(key=lambda i: itera[i - 1], reverse=reverse)
        return Permutation(L)

    @classmethod
    def from_code(cls, cd):
        return uncode(cd)

    def _latex(self, printer):
        if Permutation.print_as_code:
            return printer.doprint(trimcode(self))
        return printer.doprint(list(self._perm))

    # pattern is a list, not a permutation
    def has_pattern(self, pattern):
        if self == Permutation(pattern):
            return True
        if len(self._perm) <= len(Permutation(pattern)):
            return False
        expanded = list(self) + list(range(len(self) + 1, len(pattern) + 1))
        for i in range(len(expanded)):
            rmval = expanded[i]
            perm2 = [*expanded[:i], *expanded[i + 1 :]]
            perm2 = tuple([val - 1 if val > rmval else val for val in perm2])
            if Permutation(perm2).has_pattern(pattern):
                return True
        return False

    def _sympystr(self, printer):
        if Permutation.print_as_code:
            return printer.doprint(trimcode(self))
        return printer.doprint(self._perm)

    def __call__(self, i):
        """1-indexed"""
        return self[i - 1]

    def zero_indexed_descents(self):
        desc = set()
        for i in range(len(self._perm) - 1):
            if self[i] > self[i + 1]:
                desc.add(i)
        return desc

    def descents(self, zero_indexed=True):
        if zero_indexed:
            return self.zero_indexed_descents()
        return {i + 1 for i in self.zero_indexed_descents()}

    def get_cycles(self):
        return self.get_cycles_cached()

    @cache
    def get_cycles_cached(self):
        return [tuple(cyclic_sort([i + 1 for i in c])) for c in spp.Permutation([k - 1 for k in self._perm]).cyclic_form]

    @classmethod
    def from_cycles(cls, cycle_iter):
        spoing = spp.Permutation(*cycle_iter)
        return cls([a + 1 for a in spoing.array_form])

    @property
    def code(self):
        return list(self._cached_code())

    @cache
    def _cached_code(self):
        return old_code(self._perm)

    @property
    def trimcode(self):
        if self._perm == ():
            return []
        return self.code[:max(self.descents(False))]

    def mul_dominant(self):
        return (~((~self).minimal_dominant_above()))

    def shiftup(self, k):
        return Permutation.from_code(k * [0] + self.code)

    @cached_property
    def inv(self):
        return sum(self.code)

    def __reduce__(self):
        return (self.__class__, (self._perm,))

    def swap(self, i, j):
        new_perm = [*self._perm]
        # print(f"SWAP {new_perm=}")
        if i > j:
            i, j = j, i
        if j >= len(new_perm):
            new_perm.extend(range(len(new_perm) + 1, j + 2))
        new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
        return Permutation(new_perm)

    def rslice(self, start, stop):
        ttup = [*self._perm, *list(range(len(self._perm) + 1, stop + 2))]
        return ttup[start:stop]

    def __getitem__(self, i):
        try:
            return self._perm[i]
        except Exception:
            if isinstance(i, slice):
                return [self[ii] for ii in range(i.start if i.start is not None else 0, i.stop if i.stop is not None else len(self))]
            if i >= len(self._perm):
                return i + 1

    def __setitem__(self, i, v):
        raise NotImplementedError

    def __hash__(self):
        return self._hash_code

    def __mul__(self, other):
        if len(other._perm) > len(self._perm):
            return Permutation([self[other._perm[i] - 1] for i in range(len(other._perm))])
        return Permutation([*[self._perm[other._perm[i] - 1] for i in range(len(other._perm))], *self._perm[len(other._perm) :]])

    def __iter__(self):
        yield from self._perm.__iter__()

    def __getslice__(self, i, j):
        return self._perm[i:j]

    def __str__(self):
        return str(self._perm)

    def __add__(self, other):
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*self._perm, *other]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def _sympyrepr(self, printer):  # noqa: ARG002
        return f"Permutation({list(self._perm)})"

    def __radd__(self, other):
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*other, *self._perm]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def __eq__(self, other):
        if isinstance(other, Permutation):
            # print(f"{other._perm= } {self._perm=} {type(self._perm)=}")
            # return other._perm == self._perm
            return other._unique_key == self._unique_key
        if isinstance(other, list):
            # print(f"{[*self._perm]= } {other=}")
            return [*self._perm] == other
        if isinstance(other, tuple):
            # print(f"{self._perm=} {other=}")
            return self._perm == other
        return False

    def __len__(self):
        # print("REMOVE THIS")
        return max(len(self._perm), 2)

    def __invert__(self):
        new_perm = [0] * len(self._perm)
        for i in range(len(self._perm)):
            new_perm[self[i] - 1] = i + 1
        return Permutation(new_perm)

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return tuple(self) < tuple(other)

    def pattern_at(self, *indices):
        indices = sorted(indices)
        seq = [self[i] for i in indices]
        return ~Permutation.sorting_perm(seq)

    def minimal_dominant_above(self):
        return uncode(theta(self))


def inv(perm):
    return perm.inv


def code(perm):
    return perm.code


def mulperm(perm1, perm2):
    return perm1 * perm2


def uncode(cd):
    cd2 = [*cd]
    if cd2 == []:
        return Permutation([])
    max_required = max([cd2[i] + i for i in range(len(cd2))])
    cd2 += [0 for i in range(len(cd2), max_required)]
    fullperm = [i + 1 for i in range(len(cd2) + 1)]
    perm = []
    for i in range(len(cd2)):
        perm += [fullperm.pop(cd2[i])]
    perm += [fullperm[0]]
    return Permutation(perm)


def inverse(perm):
    return ~perm


def permtrim(perm):
    return Permutation(perm)


def strict_theta(u):
    ret = [*trimcode(u)]
    did_one = True
    while did_one:
        did_one = False
        for i in range(len(ret) - 2, -1, -1):
            if ret[i + 1] != 0 and ret[i] <= ret[i + 1]:
                ret[i], ret[i + 1] = ret[i + 1] + 1, ret[i]
                did_one = True
                break
    while len(ret) > 0 and ret[-1] == 0:
        ret.pop()
    return ret


def longest_element(indices):
    perm = Permutation([1, 2])
    did_one = True
    while did_one:
        did_one = False
        for i in range(len(indices)):
            j = indices[i] - 1
            if sg(j, perm) == 0:
                perm = perm.swap(j, j + 1)
                did_one = True
    return permtrim(perm)


def theta(perm):
    cd = code(perm)
    for i in range(len(cd) - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            if cd[j] < cd[i]:
                cd[i] += 1
    cd.sort(reverse=True)
    return cd


def trimcode(perm):
    cd = perm.code
    while len(cd) > 0 and cd[-1] == 0:
        cd.pop()
    return cd


def cycle(p, q):
    return Permutation(list(range(1, p)) + [i + 1 for i in range(p, p + q)] + [p])


def phi1(u):
    c_star = (~u).code
    c_star.pop(0)
    # print(f"{uncode(c_star)=}")
    return ~(uncode(c_star))


def one_dominates(u, w):
    c_star_u = (~u).code
    c_star_w = (~w).code

    a = c_star_u[0]
    b = c_star_w[0]

    for i in range(a, b):
        if i >= len(u) - 1:
            return True
        if u[i] > u[i + 1]:
            return False
    return True


def dominates(u, w):
    u2 = u
    w2 = w
    while inv(u2) > 0 and one_dominates(u2, w2):
        u2 = phi1(u2)
        w2 = phi1(w2)
    if inv(u2) == 0:
        return True
    return False


def medium_theta(perm):
    cd = code(perm)
    found_one = True
    while found_one:
        found_one = False
        for i in range(len(cd) - 1):
            if cd[i] < cd[i + 1]:
                found_one = True
                cd[i], cd[i + 1] = cd[i + 1] + 1, cd[i]
                break
            if cd[i] == cd[i + 1] and cd[i] != 0 and i > 0 and cd[i - 1] <= cd[i] + 1:
                cd[i] += 1
                found_one = True
                break
    return cd


def split_perms(perms):
    perms2 = [perms[0]]
    for perm in perms[1:]:
        cd = code(perm)
        index = -1
        not_zero = False
        did = False
        for i in range(len(cd)):
            if cd[i] != 0:
                not_zero = True
            elif not_zero and cd[i] == 0:
                not_zero = False
                index = i
                num_zeros_to_miss = 0
                for j in range(index):
                    if cd[j] != 0:
                        num_zeros_to_miss = max(num_zeros_to_miss, cd[j] - (index - 1 - j))
                num_zeros = 0
                for j in range(index, len(cd)):
                    if cd[j] != 0:
                        break
                    num_zeros += 1
                if num_zeros >= num_zeros_to_miss:
                    cd1 = cd[:index]
                    cd2 = [0 for i in range(index)] + cd[index:]
                    perms2 += [
                        uncode(cd1),
                        uncode(cd2),
                    ]
                    did = True
                    break
        if not did:
            perms2 += [perm]
    return perms2


class NilPlactic:
    def __init__(self, word):
        self._word = tuple(word)

    @staticmethod
    def from_word(word):
        if len(word) <= 1:
            return NilPlactic(word)
        return NilPlactic.from_word(word[:-1]).ed_insert(word[-1])

    def __hash__(self):
        return hash(self._word)

    @staticmethod
    def _ed_insert(word, letter):
        if len(word) == 0:
            return ([letter],)
        # print(word)
        row = word[0]
        index = 0
        while index < len(row) and row[index] < letter:
            index += 1
        if index == len(row):
            return [(*row, letter), *word[1:]]
        if row[index] == letter:
            if index < len(row) - 1:
                if row[index + 1] == letter + 1:
                    # return tuple([*NilPlactic._ed_insert(word[:first_row_start], letter + 1), *word[first_row_start:]])
                    return [row, *NilPlactic._ed_insert(word[1:], letter + 1)]
            while index < len(row) and row[index] == letter:
                index += 1
            if index == len(row):
                return [(*row, letter), *word[1:]]
        new_word = [*word]
        bump = row[index]
        new_word[0] = [*new_word[0]]
        new_word[0][index] = letter
        new_word[0] = (*new_word[0],)
        # return tuple([*NilPlactic._ed_insert(new_word[:first_row_start], bump), *new_word[first_row_start:]])
        return [new_word[0], *NilPlactic._ed_insert(new_word[1:], bump)]

    @staticmethod
    def ed_insert_rsk(word, word2, letter, letter2):
        if len(word) == 0:
            return (letter,), (letter2,)
        first_row_start = len(word) - 1
        if len(word) > 1:
            while first_row_start > 0 and word[first_row_start - 1] < word[first_row_start]:
                first_row_start -= 1
        index = first_row_start
        while index < len(word) and word[index] < letter:
            index += 1
        if index == len(word):
            return (*word, letter), (*word2, letter2)
        if word[index] == letter:
            word0, word2_0 = NilPlactic.ed_insert_rsk(word[:first_row_start], word2[:first_row_start], letter + 1, letter2)
            return (
                (*word0, *word[first_row_start:]),
                (*word2_0, *word2[first_row_start:]),
            )
        new_word = [*word]
        new_word2 = [*word2]
        bump = new_word[index]
        # bump2 = new_word2[index]
        new_word[index] = letter
        # new_word2[index] = letter2
        word0, word2_0 = NilPlactic.ed_insert_rsk(new_word[:first_row_start], new_word2[:first_row_start], bump, letter2)
        return (*word0, *new_word[first_row_start:]), (*word2_0, *new_word2[first_row_start:])

    @staticmethod
    def standardize(word):
        if len(word) == len(set(word)) and set(word) == set(range(1, len(word) + 1)):
            return word
        index = 0
        while index > 0 and word[index - 1] <= word[index]:
            index -= 1
        subword = [*word[:index], *word[index + 1 :]]
        subword = [i + 1 for i in NilPlactic.standardize(subword)]
        subword.insert(index, 1)
        return subword

    @staticmethod
    def inverse_ed_insert_rsk(word, word2):
        # print(f"DBG: {word=}, {word2=}")
        if len(word) != len(word2):
            raise ValueError("Words must be of the same length for inverse ed insert.")
        if len(word) == 0:
            return (), ()
        sputnik = NilPlactic.standardize(word2)
        # print(f"DBG: {sputnik=}")
        index = sputnik.index(max(sputnik))
        new_word = [*word]
        new_word2 = [*word2]
        a, b = new_word.pop(index), new_word2.pop(index)
        index2 = index
        while index2 < len(word2) - 1 and word[index2] < word[index2 + 1]:
            index2 += 1
        if index2 == len(word2) - 1:
            new_word, new_word2 = NilPlactic.inverse_ed_insert_rsk(new_word, new_word2)
            return (*new_word, a), (*new_word2, b)
        row_start = index2 + 1
        index3 = row_start
        while index3 < len(word2) - 1 and word[index3] < word[index3 + 1] and word[index3] < a:
            index3 += 1
        if index3 < len(word2) - 1:
            if word[index3] < a:
                # print("pangolin")
                # print(f"{word=}, {word2=}, {index3=}, {a=}, {new_word=}, {new_word2=} {b=}")
                if word[index3] == a - 1:
                    new_word[index3 - 1] = a
                    new_word, new_word2 = NilPlactic.inverse_ed_insert_rsk(new_word, new_word2)
                    return (*new_word, a - 1), (*new_word2, b)
                raise ValueError(f"Cannot perform inverse ed insert: word is not nilplactic. DBG {index3=}, {word=}, {word2=}, {a=}")
            if word[index3] == a:
                # print("pangolin2")
                # print(f"{word=}, {word2=}, {index3=}, {a=}, {new_word=}, {new_word2=} {b=}")
                new_word, new_word2 = NilPlactic.inverse_ed_insert_rsk(new_word, new_word2)
                return (*new_word, a - 1), (*new_word2, b)

        a2 = word[index3]
        new_word[index3 - 1] = a
        new_word, new_word2 = NilPlactic.inverse_ed_insert_rsk(new_word, new_word2)
        return (*new_word, a2), (*new_word2, b)

    @staticmethod
    def reverse_insert_rsk(word, word2, letter, letter2):
        if len(word2) == 0:
            return (letter,), (letter2,)
        first_row_start = len(word2) - 1
        if len(word2) > 1:
            while first_row_start > 0 and word2[first_row_start - 1] < word2[first_row_start]:
                first_row_start -= 1
        index = first_row_start
        while index < len(word2) and word2[index] <= letter2:
            index += 1
        if index == len(word2):
            return (*word, letter), (*word2, letter2)
        new_word = [*word]
        new_word2 = [*word2]
        bump = new_word2[index]
        # bump2 = new_word2[index]
        new_word2[index] = letter2
        # new_word2[index] = letter2
        word0, word2_0 = NilPlactic.reverse_insert_rsk(new_word[:first_row_start], new_word2[:first_row_start], letter, bump)
        return (*word0, *new_word[first_row_start:]), (*word2_0, *new_word2[first_row_start:])

    def ed_insert(self, letter):
        """Insert a letter into the nilplactic word."""
        return NilPlactic(NilPlactic._ed_insert(self._word, letter))

    def inverse_insert(self, position):
        if position >= len(self._word):
            raise IndexError("Position out of bounds for nilplactic word.")
        row_start = position
        row_end = position
        if len(self._word) > 1:
            while row_start > 0 and self._word[row_start - 1] < self._word[row_start]:
                row_start -= 1
            while row_end < len(self._word) - 1 and self._word[row_end] < self._word[row_end + 1]:
                row_end += 1
        row = self._word[row_start : row_end + 1]
        pos = position - row_start
        return row.pop(pos)

    def __eq__(self, other):
        if isinstance(other, NilPlactic):
            return self._word == other._word
        return False

    def __repr__(self):
        return f"{self._word}"

    def __str__(self):
        return f"{self._word}"


class Plactic:
    def __init__(self, word):
        self._word = tuple(word)

    @staticmethod
    def from_word(word):
        if len(word) <= 1:
            return Plactic(word)
        return Plactic.from_word(word[:-1]).rs_insert(word[-1])

    def __hash__(self):
        return hash(self._word)

    @staticmethod
    def _rs_insert(word, letter):
        if len(word) == 0:
            return (letter,)
        first_row_start = len(word) - 1
        if len(word) > 1:
            while first_row_start > 0 and word[first_row_start - 1] <= word[first_row_start]:
                first_row_start -= 1
        index = first_row_start
        while index < len(word) and word[index] <= letter:
            index += 1
        if index == len(word):
            return (*word, letter)
        new_word = [*word]
        bump = new_word[index]
        new_word[index] = letter
        return (*Plactic._rs_insert(new_word[:first_row_start], bump), *new_word[first_row_start:])

    def rs_insert(self, letter):
        """Insert a letter into the nilplactic word."""
        return Plactic(Plactic._rs_insert(self._word, letter))

    def __eq__(self, other):
        if isinstance(other, Plactic):
            return self._word == other._word
        return False

    def __repr__(self):
        return f"{self._word}"

    def __str__(self):
        return f"{self._word}"


bad_classical_patterns = [Permutation([1, 4, 2, 3]), Permutation([1, 4, 3, 2]), Permutation([4, 1, 3, 2]), Permutation([3, 1, 4, 2])]


# def perm_to_key(perm):
#     if len(perm) == 0:
#         return {NilPlactic(()): S.One}

#     ret = {}
#     stack = [(perm, NilPlactic(()), 1, S.One)]

#     while len(stack) > 0:
#         current_perm, word, index, poly = stack.pop()
#         if current_perm.inv == 0:
#             np_elem = word
#             ret[np_elem] = ret.get(np_elem, S.Zero) + poly
#             continue
#         # L = sl.pull_out_var(1, current_perm)
#         for index_list, new_perm in L:
#             index_list.sort(reverse=True)
#             new_word = word
#             for index2 in index_list:
#                 new_word = new_word.ed_insert(index + index2 - 1)
#             stack.append((new_perm, new_word, index + 1, poly * prod([x[index] - y[a] for a in index_list])))
#     return ret

ID_PERM = Permutation([])


@cache
def s(i):
    return Permutation([*list(range(1, i)), i + 1, i])
