from schubmult.perm_lib import Permutation, uncode
from schubmult.symbolic import sympify
from schubmult.utils.perm_utils import add_perm_dict

from .free_algebra import FreeAlgebra, FreeAlgebraElement
from .free_algebra_basis import SchubertBasis, WordBasis
from .tensor_ring import TensorRingElement

FAS = FreeAlgebra(basis=SchubertBasis)

class RCGraph(tuple):

    def __new__(cls, *args):
        obj = tuple.__new__(cls, *args)
        return obj

    @property
    def perm(self):
        perm = Permutation([])
        for row in self:
            for p in row:
                perm = perm.swap(p - 1, p)
        return perm

    def act(self, p):
        pm = self.perm
        elem = FAS(pm, len(self))
        bumpup = FAS(uncode([p]), 1)* elem
        ret = set()
        for k, v in bumpup.items():
            perm2 = k[0]
            new_row = [pm[i] for i in range(max(len(pm), len(perm2))) if pm[i] == perm2[i+1]]
            new_row.sort(reverse=True)
            nrc = RCGraph([tuple(new_row), *[tuple([row[i] + 1 for i in range(len(row))]) for row in self]])
            assert nrc.perm == perm2
            # except AssertionError:
            #     print(self)
            #     print(perm2)
            #     print(nrc.perm)
            #     print(nrc)
            #     print(f"{new_row=}")
            #     raise
            ret.add(nrc)
        return ret

    def as_str_lines(self):
        lines = []
        for i in range(len(self)):
            row = self[i]
            splug_row = [*row, i]
            row = [str(splug_row[j]) + ("  ")*(splug_row[j] - splug_row[j+1] - 1) for j in range(len(splug_row)-1)]
            line = ""
            line += " ".join([str(r) for r in row])
            lines += [line]
        if len(lines) == 0:
            lines = ["" ]
        lines2 = []
        ml = max([len(line) for line in lines])
        if len(lines[0]) < ml:
            lines[0] = " " * (ml - len(lines[0])) + lines[0]
        for line in lines:
            lines2 += [" "*(ml - len(line)) + line]
        return lines2

    def __str__(self):
        return "\n".join(self.as_str_lines())

    def __hash__(self):
        return hash((tuple(self), "RCGRAPH"))




# rc graph module is a module of tuples
# backwards inserted, left action by free algebra


class RCGraphModule(dict):

    def __add__(self, other):
        if isinstance(other, RCGraphModule):
            return RCGraphModule(add_perm_dict(self, other))
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, FreeAlgebraElement):
            wd_dict = other.change_basis(WordBasis)
            ret = {}
            for k0, v0 in self.items():
                for k, v in wd_dict.items():
                    addup = RCGraphModule({k0: v0})
                    for a in reversed(k):
                        new_addup = RCGraphModule()
                        for kr, v in addup.items():
                            new_addup += RCGraphModule({r: v for r in kr.act(a)})
                        addup = new_addup

                    ret = add_perm_dict(ret, new_addup)
            return RCGraphModule(ret)
        try:
            other = sympify(other)
            return RCGraphModule({k: v * other for k, v in self.items()})
        except Exception:
            return NotImplemented

    def as_str_lines(self):
        lines = [""]
        first = True
        for k, v in self.items():
            lines2 = k.as_str_lines()
            if len(lines) < len(lines2):
                upstr = ""
                if len(lines[0]) > 0:
                    upstr = " " * len(lines[0])
                lines += [upstr] * (len(lines2) - len(lines))
            padlen = 0
            for i in range(len(lines2)):
                coeffstr = ""
                if not first:
                    if i == 0:
                        coeffstr += " + "
                    else:
                        coeffstr += "   "
                if i == 0:
                    if v == -1:
                        coeffstr += "-"
                    elif v != 1:
                        coeffstr += str(v) + " * "
                    padlen = len(coeffstr)
                else:
                    coeffstr = " " * padlen

                lines[i] += coeffstr + lines2[i]
            first = False
        return lines

    def __str__(self):
        return "\n".join(self.as_str_lines())

class RCGraphTensor(tuple):

    def __new__(cls, graph1, graph2):
        obj = tuple.__new__(cls, (graph1, graph2))
        return obj

    def as_str_lines(self):
        lines1 = self[0].as_str_lines()
        lines2 = self[1].as_str_lines()

        ml = max(len(lines1), len(lines2))
        if len(lines1) < ml:
            lines1 = [" " * len(lines1[0])] * (ml - len(lines1)) + lines1
        if len(lines2) < ml:
            lines2 = [" " * len(lines2[0])] * (ml - len(lines2)) + lines2

        lines = [lines1[0] + "  #  " + lines2[0]]
        for i in range(1, ml):
            lines += [lines1[i] + "     " + lines2[i]]

        return lines

    def __str__(self):
        return "\n".join(self.as_str_lines())

    def __hash__(self):
        return hash((tuple(self), "RCGRAPHTENSOR"))

class TensorModule(RCGraphModule):
    def __add__(self, other):
        if isinstance(other, TensorModule):
            return TensorModule(add_perm_dict(self, other))
        return NotImplemented

    @classmethod
    def ext_multiply(cls, elem1, elem2):
        ret = cls()
        for key, val in elem1.items():
            for key2, val2 in elem2.items():
                ret += cls({RCGraphTensor(key, key2): val * val2})
        return ret

    def __rmul__(self, other):
        if isinstance(other, TensorRingElement):
            ret = {}
            for k0, v0 in self.items():
                for k, v in other.items():
                    addup = TensorModule({k0: v0})
                    new_addup = TensorModule()
                    for kr, v in addup.items():
                        elem1 = other.ring.rings[0](*k[0]) * RCGraphModule({kr[0]: 1})
                        elem2 = other.ring.rings[1](*k[1]) * RCGraphModule({kr[1]: 1})
                        new_addup += v * TensorModule.ext_multiply(elem1, elem2)
                    addup = new_addup

                    ret = add_perm_dict(ret, addup)
            return TensorModule(ret)
        try:
            other = sympify(other)
            return TensorModule({k: v * other for k, v in self.items()})
        except Exception:
            return NotImplemented



if __name__ == "__main__":
    FA = FreeAlgebra(basis=SchubertBasis)

    spug = TensorModule({RCGraphTensor(RCGraph(()), RCGraph(())): 1})
    
    # spug1 = FA(uncode([0,1]), 2).coproduct() * spug
    # print(spug1)

    # spug2 = FA(uncode([1,0]), 2).coproduct() * spug
    # print(spug2)
    res = TensorModule({})
    spurg = (FA @ FA).zero
    for i in range(2):
        for j in range(2):
            oil =  FA(uncode([i]),1).coproduct() * FA(uncode([j]),1).coproduct()
            spug2 =oil * spug
            print(oil)
            print(spug2)
            res += spug2
            spurg += oil
    # print(res)
    print(spurg)