import math
from functools import cache, cached_property

import sympy.combinatorics.permutations as spp
from symengine import sympify
from sympy import Basic, Tuple

import schubmult.utils.logging as lg
from schubmult.utils.perm_utils import cyclic_sort, permtrim_list, sg

# schubmult.poly_lib.variables import GeneratingSet

logger = lg.get_logger(__name__)

zero = sympify(0)
n = 100

# TODO: permutations act


class Permutation(Basic):
    def __new__(cls, perm):
        return Permutation.__xnew_cached__(cls, tuple(perm))

    print_as_code = False

    @staticmethod
    @cache
    def __xnew_cached__(_class, perm):
        return Permutation.__xnew__(_class, perm)

    @staticmethod
    def __xnew__(_class, perm):
        p = tuple(permtrim_list([*perm]))
        s_perm = spp.Permutation._af_new([i - 1 for i in p])
        obj = Basic.__new__(_class, Tuple(*perm))
        obj._s_perm = s_perm
        obj._perm = p
        obj._hash_code = hash(p)
        cd = s_perm.inversion_vector()
        obj._unique_key = (len(p), sum([cd[i] * math.factorial(len(p) - 1 - i) for i in range(len(cd))]))
        return obj

    @classmethod
    def sorting_perm(cls, itera):
        L = [i + 1 for i in range(len(itera))]
        L.sort(key=lambda i: itera[i - 1])
        return Permutation(L)

    def _sympystr(self, printer):
        if Permutation.print_as_code:
            return printer.doprint(trimcode(self))
        return printer.doprint(self._perm)

    def __call__(self, i):
        """1-indexed"""
        return self[i - 1]

    def descents(self, zero_indexed=True):
        if zero_indexed:
            return self._s_perm.descents()
        return {i + 1 for i in self._s_perm.descents()}

    def get_cycles(self):
        return self.get_cycles_cached()

    @cache
    def get_cycles_cached(self):
        return [tuple(cyclic_sort([i + 1 for i in c])) for c in self._s_perm.cyclic_form]

    @property
    def code(self):
        return list(self.cached_code())

    @cache
    def cached_code(self):
        return self._s_perm.inversion_vector()

    @cached_property
    def inv(self):
        return self._s_perm.inversions()

    def swap(self, i, j):
        new_perm = [*self._perm]
        # print(f"SWAP {new_perm=}")
        if i > j:
            i, j = j, i
        if j >= len(new_perm):
            # print(f"SWAP {j}>={new_perm=}")
            new_perm += list(range(len(new_perm) + 1, j + 2))
            # print(f"SWAP extended {new_perm=}")
        new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
        # print(f"SWAP iddle {new_perm=}")
        return Permutation(new_perm)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return [self[ii] for ii in range(i.start if i.start is not None else 0, i.stop if i.stop is not None else len(self))]
        if i >= len(self._perm):
            return i + 1
        return self._perm[i]

    def __setitem__(self, i, v):
        raise NotImplementedError

    def __hash__(self):
        return self._hash_code

    def __mul__(self, other):
        new_sperm = other._s_perm * self._s_perm
        new_perm = permtrim_list([new_sperm.array_form[i] + 1 for i in range(new_sperm.size)])
        return Permutation(new_perm)

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
        new_sperm = ~(self._s_perm)
        new_perm = [new_sperm.array_form[i] + 1 for i in range(new_sperm.size)]
        return Permutation(new_perm)

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return tuple(self) < tuple(other)


def ensure_perms(func):
    def wrapper(*args):
        return func(*[Permutation(arg) if (isinstance(arg, list) or isinstance(arg, tuple)) else arg for arg in args])

    return wrapper


@ensure_perms
def inv(perm):
    return perm.inv


@ensure_perms
def code(perm):
    return perm.code


@ensure_perms
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


@ensure_perms
def inverse(perm):
    return ~perm


def permtrim(perm):
    return Permutation(perm)


@ensure_perms
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


@ensure_perms
def theta(perm):
    cd = code(perm)
    for i in range(len(cd) - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            if cd[j] < cd[i]:
                cd[i] += 1
    cd.sort(reverse=True)
    return cd


@ensure_perms
def trimcode(perm):
    cd = perm.code
    while len(cd) > 0 and cd[-1] == 0:
        cd.pop()
    return cd


def cycle(p, q):
    return Permutation(list(range(1, p)) + [i + 1 for i in range(p, p + q)] + [p])


@ensure_perms
def phi1(u):
    c_star = (~u).code
    c_star.pop(0)
    # print(f"{uncode(c_star)=}")
    return ~(uncode(c_star))


@ensure_perms
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
