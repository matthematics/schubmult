from schubmult.perm_lib import Permutation, uncode
from schubmult.symbolic import S, sympify
from schubmult.utils.perm_utils import add_perm_dict

from .free_algebra import FreeAlgebra, FreeAlgebraElement
from .free_algebra_basis import SchubertBasis, WordBasis
from .tensor_ring import TensorRing, TensorRingElement

FAS = FreeAlgebra(basis=SchubertBasis)

class RCGraph(tuple):

    def __new__(cls, *args):
        obj = tuple.__new__(cls, *args)
        return obj

    def polyvalue(self, x):
        ret = S.One
        for i, row in enumerate(self):
            ret *= x[i+1]**len(row)
        return ret

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

    def __new__(cls, *args):
        obj = dict.__new__(cls)
        dct = dict(*args)
        obj.update({k: v for k, v in dct.items() if v != 0})
        return obj

    def polyvalue(self, x):
        ret = S.Zero
        for k, v in self.items():
            ret += v * k.polyvalue(x)
        return ret

    def __init__(self, *args):
        pass

    def __add__(self, other):
        if isinstance(other, RCGraphModule):
            return RCGraphModule(add_perm_dict(self, other))
        return NotImplemented

    def __neg__(self):
        return RCGraphModule({k: -v for k, v in self.items()})

    def __sub__(self, other):
        if isinstance(other, RCGraphModule):
            return RCGraphModule(add_perm_dict(self, -other))
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, FreeAlgebraElement):
            wd_dict = other.change_basis(WordBasis)
            ret = {}
            for k0, v0 in self.items():
                for k, v in wd_dict.items():
                    addup = RCGraphModule({k0: v0*v})
                    for a in reversed(k):
                        new_addup = RCGraphModule()
                        for kr, vv in addup.items():
                            new_addup += RCGraphModule(dict.fromkeys(kr.act(a), vv))
                        addup = new_addup

                    ret = add_perm_dict(ret, new_addup)
            return RCGraphModule(ret)
        try:
            other = sympify(other)
            return self.__class__({k: v * other for k, v in self.items()})
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

    def polyvalue(self, x):
        return self[0].polyvalue(x) * self[1].polyvalue(x)


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

    def __new__(cls, *args):
        return RCGraphModule.__new__(cls, *args)

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
                    addup = TensorModule({k0: v*v0})
                    new_addup = TensorModule()
                    for kr, vv in addup.items():
                        elem1 = other.ring.rings[0](*k[0]) * RCGraphModule({kr[0]: 1})
                        if isinstance(other.ring.rings[1], TensorRing):
                            elem2 = other.ring.rings[1](k[1]) * TensorModule({kr[1]: 1})
                        else:
                            elem2 = other.ring.rings[1](*k[1]) * RCGraphModule({kr[1]: 1})

                        new_addup += vv * TensorModule.ext_multiply(elem1, elem2)
                    addup = new_addup

                    ret = add_perm_dict(ret, addup)
            return TensorModule(ret)
        try:
            other = sympify(other)
            return TensorModule({k: v * other for k, v in self.items()})
        except Exception:
            return NotImplemented



if __name__ == "__main__":
    from schubmult import Sx
    from schubmult.abc import x
    from schubmult.rings.free_algebra_basis import *
    FA = FreeAlgebra(basis=SchubertBasis)
    FAW = FreeAlgebra(basis=WordBasis)

    #spug = TensorModule({RCGraphTensor(RCGraph(()), RCGraph(())): 1})
    spug = RCGraphModule({RCGraph(()): 1})
    
    # spug1 = FA(uncode([0,1]), 2).coproduct() * spug
    # print(spug1)

    # spug2 = FA(uncode([1,0]), 2).coproduct() * spug
    # print(spug2)
    #resspug = TensorModule({RCGraphTensor(RCGraph(()), RCGraphTensor(RCGraph(), RCGraph(()))): 1})

    res = TensorModule({RCGraphTensor(RCGraph(()), RCGraph(())): 1})
    respo = TensorModule()
    #spurg = (FA @ FA).zero
    bungalo = {}
    fratboy = {}
    ring = FAW @ (FAW @ FAW)
    mx = 2
    for i in range(mx + 1):
        for j in range(mx + 1):
            #oil =  FA(uncode([i]),1).coproduct() * FA(uncode([j]),1).coproduct()
            #oil = ring.ext_multiply(FAW(i,j),FAW(i,j).coproduct())
            oil = FAW(i,j).coproduct()
            #spud_dog = oil.change_basis(SchubertBasis)
            respo += oil * res
            fratboy[(i,j)] = fratboy.get((i,j), RCGraphModule()) + oil * res
            # for k, v in spud_dog.items():
            #     fratboy[k] = fratboy.get(k, RCGraphModule()) + v*spug2
            #print(f"{FreeAlgebraBasis.change_tensor_basis(oil,WordBasis,WordBasis)=}")
            # print(spug2)
            # print(f"{dict(spug2)=}")
            #res += spug2
            #spurg += oil
    foingle = {}
    for k, v in fratboy.items():
        foingle2 = FAW(*k).change_basis(SchubertBasis)
        for kk, vv in foingle2.items():
            foingle[kk] = foingle.get(kk, RCGraphModule()) + vv * v
    # print(res)
    # for bacon, v in fratboy.items():
    #     print(f"{bacon=}")
    #     print(v)
    Permutation.print_as_code = True
    for k, v in foingle.items():
        if any(a > mx for a in k[0].trimcode):
            continue
        print(k[0].trimcode)
        print(v)
        print(v.polyvalue(x))
        sm0 = S.Zero
        sm1 = S.Zero
        side0 = {}
        side1 = {}
        plathbucket = {}
        for kk, vv in v.items():
            if not kk[0].perm.bruhat_leq(k[0]) or not kk[1].perm.bruhat_leq(k[0]):
                continue
            #     print("EXCLUDE")
            #     print(kk)
            #     print("END EXCLUDE")
            #     continue
            print("GOOD GRAPH")
            print(kk)
            
            key = (kk[0].perm, kk[1].perm)
            side0[key[0]] = side0.get(key[0], 0) + vv*kk[0].polyvalue(x)
            side1[key[1]] = side1.get(key[1], 0) + vv*kk[1].polyvalue(x)
            plathbucket[key] = plathbucket.get(key, 0) + vv*kk[0].polyvalue(x)*kk[1].polyvalue(x)
            # if kk[0].perm.inv ==0:
            #     sm1 += vv * kk[1].polyvalue(x)
            # elif kk[1].perm.inv ==0:
            #     sm0 += vv * kk[0].polyvalue(x)
        for perm, v in side0.items():
            print(f"  L {perm} : {v} {Sx(v)}")
        for perm, v in side1.items():
            print(f"  R {perm} : {v} {Sx(v)}")
        for (p0, p1), v in plathbucket.items():
            print(f"  P {p0}, {p1} : {v} {Sx(v)}")
        
        print("-----")
