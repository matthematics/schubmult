import math
from functools import cache, cached_property

import sympy.combinatorics.permutations as spp
from sympy.printing.defaults import Printable

import schubmult.utils.logging as lg
import schubmult.utils.perm_utils as sl

logger = lg.get_logger(__name__)

zero = 0
n = 100


class Permutation(Printable):
    """Permutation class representing permutations of positive integers."""

    def act_root(self, a, b):
        return self[a - 1], self[b - 1]

    @property
    def antiperm(self):
        w0 = Permutation.w0(len(self))
        return w0 * (self) * w0

    def coset_decomp(self, *descs):
        descs = set(descs)
        reduced_perm = self
        w_J = Permutation([])
        found = True
        while found:
            found = False
            for d in reduced_perm.descents():
                if d + 1 in descs:
                    w_J = ~((~w_J).swap(d, d + 1))
                    reduced_perm = reduced_perm.swap(d, d + 1)
                    found = True
                    break
        return reduced_perm, w_J

    def min_coset_rep(self, *descs):
        return self.coset_decomp(*descs)[0]

    def max_coset_rep(self, *descs):
        red, w_J = self.coset_decomp(*descs)
        return red * Permutation.longest_element(*descs)

    @classmethod
    def longest_element(cls, *descs):
        perm = Permutation([])
        did_one = True
        while did_one:
            did_one = False
            for i in range(len(descs)):
                j = descs[i] - 1
                if perm[j] < perm[j + 1]:
                    perm = perm.swap(j, j + 1)
                    did_one = True
        return perm

    @classmethod
    def w0(cls, n):
        return cls.from_code([n - 1 - i for i in range(n - 1)])

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

    @property
    def code_word(self):
        cd = self.trimcode
        word = []
        for i in range(len(cd)):
            word += list(range(i + cd[i], i, -1))
        return tuple(word)

    def root_swap(self, root):
        return self.swap(root[0] - 1, root[1] - 1)

    @classmethod
    def reflection(cls, root):
        return cls([]).swap(root[0] - 1, root[1] - 1)

    def right_root_at(self, index, word=None):
        if word is None:
            word = [*self.code_word]
        word_piece = word[index + 1 :]
        # print(f"{word=}")
        # print(f"{word_piece=}")
        apply = ~Permutation.ref_product(*word_piece)
        root = apply.act_root(word[index], word[index] + 1)
        # print(f"{root=}")
        return root

    def code_index_of_index(self, index):
        running_sum = 0
        running_code_index = 0
        for code_index, code_elem in enumerate(self.trimcode):
            if code_elem == 0:
                continue
            running_sum += code_elem
            if running_sum > index:
                return running_code_index
            running_code_index += 1
        return len(self.trimcode)

    @staticmethod
    def cycle(p, q):
        """
        Construct the cycle permutation used elsewhere in the code.
        Kept as a staticmethod on Permutation for call sites like Permutation.cycle(p,q).
        """
        return Permutation(list(range(1, p)) + [i + 1 for i in range(p, p + q)] + [p])

    @staticmethod
    @cache
    def __xnew_cached__(_class, perm):
        return Permutation.__xnew__(_class, perm)

    @staticmethod
    def __xnew__(_class, perm):
        p = tuple(sl.permtrim_list([*perm]))
        # s_perm = spp.Permutation([i - 1 for i in p])
        obj = object.__new__(_class)
        obj._args = (p,)
        # obj._s_perm = tuple([i - 1 for i in p])
        obj._perm = p
        obj._hash_code = hash(p)
        cd = sl.old_code(p)
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

    def bruhat_leq(perm, perm2):
        if perm.inv == perm2.inv:
            return perm == perm2
        if perm.inv > perm2.inv:
            return False
        ml = max(len(perm), len(perm2))
        full_perm = [perm[i] for i in range(ml)]
        full_perm2 = [perm2[i] for i in range(ml)]
        for i in range(1, ml):
            arr1 = list(full_perm[:i])
            arr2 = list(full_perm2[:i])
            arr1.sort()
            arr2.sort()
            if not (arr1 <= arr2):
                return False
        return True

    @classmethod
    def from_code(cls, cd):
        return uncode(cd)

    # def _latex(self, printer):
    #     if Permutation.print_as_code:
    #         return printer._print(trimcode(self))
    #     return printer._print(list(self._perm))

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

    def _pretty(self, printer=None):
        return printer._print_Tuple(tuple(self))

    def _sympystr(self, printer=None):
        from sympy.printing.str import StrPrinter

        if printer is None:
            printer = StrPrinter()
        if Permutation.print_as_code:
            return printer.doprint(trimcode(self))
        return printer.doprint(tuple(self._perm))

    def __call__(self, *tup):
        if len(tup) == 1:
            if isinstance(tup[0], (list, tuple)):
                tup = tup[0]
            else:
                return self._perm[tup[0] - 1]
        return tuple(self[i - 1] for i in tup)

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
        return [tuple(sl.cyclic_sort([i + 1 for i in c])) for c in spp.Permutation([k - 1 for k in self._perm]).cyclic_form]

    @classmethod
    def from_cycles(cls, cycle_iter):
        spoing = spp.Permutation(*cycle_iter)
        return cls([a + 1 for a in spoing.array_form])

    @property
    def code(self):
        return list(self._cached_code())

    @cache
    def _cached_code(self):
        return sl.old_code(self._perm)

    @property
    def graph(self):
        return {(i + 1, self[i]) for i in range(len(self._perm))}

    @property
    def shape(self):
        return tuple(sorted(self.code, reverse=True))

    @property
    def inversion_set(self):
        inv_set = set()
        for i in range(len(self._perm)):
            for j in range(i + 1, len(self._perm)):
                if self[i] > self[j]:
                    inv_set.add((i + 1, j + 1))
        return inv_set

    @property
    def diagram(self):
        diag = set()
        for i in range(len(self._perm)):
            for j in range(len(self._perm)):
                if self[i] > j + 1 and (~self)[j] > i + 1:
                    diag.add((i + 1, j + 1))
        return diag

    @property
    def maximal_corner(self):
        maxd = len(self.trimcode)
        end_spot = max(self[i] for i in range(maxd, len(self)) if self[i] < self[maxd - 1])
        return (maxd, end_spot)

    @classmethod
    def from_partial(cls, partial_perm):
        max_required = max([a for a in partial_perm if a is not None], default=len(partial_perm))
        search_space = {i for i in partial_perm if i is not None}

        # Need enough values to fill all None positions
        full_perm = [i + 1 for i in range(max(max_required, len(partial_perm))) if i + 1 not in search_space]
        perm = [*partial_perm]
        j = 0
        for i in range(len(perm)):
            if perm[i] is None:
                perm[i] = full_perm[j]
                j += 1
        return cls(perm)

    def pivots(self, a, b):
        piv = set()
        for i in range(1, len(self) + 1):
            j = self[i - 1]
            if i >= a or j >= b:
                continue
            good = True
            for i_prime in range(i, a + 1):
                if not good:
                    break
                for j_prime in range(j, b + 1):
                    if (i == i_prime and j == j_prime) or (i == a and j == b):
                        continue
                    if self[i_prime - 1] == j_prime:
                        good = False
                        break
            if good:
                piv.add((i, j))
        return piv


    @property
    def trimcode(self):
        if self._perm == ():
            return []
        return self.code[: max(self.descents(False))]

    def mul_dominant(self):
        return ~((~self).minimal_dominant_above())

    def shiftup(self, k):
        return Permutation.from_code(k * [0] + self.code)

    @cached_property
    def inv(self):
        return sum(self.code)

    @property
    def is_dominant(self):
        return self.minimal_dominant_above() == self

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

    # def __str__(self):
    #     return str(self._perm)

    def __add__(self, other):
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*self._perm, *other]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    # def _sympyrepr(self, printer):
    #     return f"Permutation({list(self._perm)})"

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

    def __le__(self, other):
        return self.bruhat_leq(other)

    # def __lt__(self, other):
    #     return self != other and self.bruhat_leq(other)

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

    @property
    def foundational_root(self):
        if self.inv == 0:
            return None
        mx = -1
        k = max(self.descents()) + 1
        for i in range(k + 1, len(self) + 1):
            if self[k - 1] > self[i - 1]:
                mx = i
        return (k, mx)


def inv(perm):
    return perm.inv


def code(perm):
    return perm.code


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
            if sl.sg(j, perm) == 0:
                perm = perm.swap(j, j + 1)
                did_one = True
    return perm


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
    # keep a thin module-level wrapper for backwards compatibility
    return Permutation.cycle(p, q)


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
