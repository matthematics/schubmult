import symengine
import sympy

from schubmult import *
from schubmult.abc import *
from schubmult.rings import DoubleSchubertRing, SingleSchubertRing
from schubmult.rings.symmetric_polynomials.functions import coeffvars, degree, is_of_func_type
from schubmult.schub_lib.positivity import compute_positive_rep
from schubmult.symbolic import Add, Mul, Pow, S, expand

r = GeneratingSet("r")
znz = {z[i]: -z[i] for i in range(20)}
subs_dict = {y[1]: S.Zero}

for i in range(1,len(y[1:])):
    subs_dict[y[i+1]] = subs_dict[y[i]] + r[i]

subs_dict2 = {r[i]: y[i+1]-y[i] for i in range(1,20)}
# symengine.var("y_(1:100)")
# symengine.var("z_(1:100)")
# bargain = E(1, 1, y_4, z_1)*E(1, 1, y_4, z_2)*E(1, 3, y_1, y_4, y_5, z_1, z_2, z_3) + E(1, 1, y_4, z_1)*E(1, 1, y_5, z_4)*E(1, 3, y_1, y_4, y_5, z_1, z_2, z_3) - E(1, 1, y_4, z_1)*E(2, 3, y_1, y_4, y_5, z_1, z_2) + E(1, 1, y_5, z_2)*E(1, 1, y_5, z_4)*E(1, 3, y_1, y_4, y_5, z_1, z_2, z_3) - E(1, 1, y_5, z_4)*E(2, 3, y_1, y_4, y_5, z_1, z_2) + E(3, 3, y_1, y_4, y_5, z_1)
sympy.init_printing(pretty_print=False)

z_ring = SingleSchubertRing(z)
zy_ring = DoubleSchubertRing(z,y)

bagel1 = [4, 1, 5, 2, 3]
porn = [5, 1, 4, 2, 3]
for k, bargain in (DSx(bagel1, elem_sym=True) * DSx(porn, "z")).items():
    if expand(bargain,func=True) == S.Zero:
        continue
    bacon = canonicalize_elem_syms(bargain)
    if isinstance(bacon, Add):
        dct = {}
        for arg in bacon.args:
            if is_of_func_type(arg, FactorialElemSym):
                monom = coeffvars(arg)[0] ** degree(arg)
            elif isinstance(arg, Mul):
                monom = S.One
                for arg2 in arg.args:
                    if is_of_func_type(arg2, FactorialElemSym):
                        monom *= coeffvars(arg2)[0] ** degree(arg2)
                    if isinstance(arg2, Pow) and is_of_func_type(arg2.args[0], FactorialElemSym):
                        monom *= coeffvars(arg2.args[0])[0] ** (degree(arg2.args[0]) * int(arg2.args[1]))
            elif isinstance(arg, Pow):
                monom = coeffvars(arg.args[0])[0] ** (degree(arg.args[0]) * int(arg.args[1]))
            dct[monom] = dct.get(monom, S.Zero) + arg
        #schuber = z_ring.from_expr(Add(*[expand(efficient_subs(vl,{y[i]: S.Zero for i in range(10)}),func=True,mul=False) for vl in dct.values()]))
        #schuber = z_ring.from_expr(expand(efficient_subs(bargain,znz),func=True,mul=False))
        # if schuber.expand() == S.Zero:
        #     print("Dodged a bullet")
        #     continue
        print(f"{k}")
        #print(schuber.ring.from_dict({k: abs(v) for k, v in schuber.items()}))
        dctyep = {}
        dctnope = {}
        dctbool = {}
        pos_part = S.Zero
        neg_part = S.Zero
        anyn = False
        for monom, bargle in dct.items():
            if expand(bargle, func=True, deep=True, mul=True) != 0:
                voib_bo = efficient_subs(expand(bargle, func=False, mul=True),znz)
                voib=expand(voib_bo,func=True,mul=False)
                plop = z_ring(voib)
                plop2 = z_ring.from_dict({k: expand(efficient_subs(v,subs_dict),mul=False) for k,v in plop.items()})
                try:
                    compute_positive_rep(expand(bargle, func=True, mul=False), y, z, False, False)
                    print(f"Yep: {monom}")#: {bargle=}")
                    #dctbool[monom] = True
                    pos_part += plop2
                    dctyep[monom] = voib_bo
                except Exception:
                    neg_part += plop2
                    print(f"Nope {monom}")#: {bargle=}")
                    anyn = True
                    dctnope[monom] = voib_bo
        if anyn:
            pos_part = pos_part.expand(deep=False)
            neg_part = neg_part.expand(deep=False)
            if pos_part == S.Zero:
                continue
            dctyep2 = {}
            dctnope2 = {}

            def forple(voib2):
                plop = zy_ring(efficient_subs(voib2,znz))
                plop2 = zy_ring.from_dict({k: expand(efficient_subs(v,subs_dict),mul=False) for k,v in plop.items()})
                return plop2
            for monom, blarp in dctyep.items():
                if isinstance(blarp, Add):
                    dctyep2[monom] = [{k: expand(efficient_subs(v,subs_dict),func=True) for k,v in z_ring.from_expr(arg).items() if k in pos_part} for arg in blarp.args]
                else:
                    dctyep2[monom] = {k: expand(efficient_subs(v, subs_dict),func=True) for k,v in z_ring.from_expr(blarp).items() if k in pos_part}
            for monom, blarp in dctnope.items():
                if isinstance(blarp, Add):
                    dctnope2[monom] = [{k: expand(efficient_subs(v,subs_dict),func=True) for k,v in z_ring.from_expr(arg).items() } for arg in blarp.args]
                else:
                    dctnope2[monom] = {k: expand(efficient_subs(v, subs_dict),func=True) for k,v in z_ring.from_expr(blarp).items()}
            print(f"{dctyep2=}")
            print(f"{dctnope2=}")

        
      # Anything that isn't Schub will disappear
