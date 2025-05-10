import symengine
import sympy

from schubmult import *
from schubmult.abc import *
from schubmult.rings import SingleSchubertRing
from schubmult.rings.symmetric_polynomials.functions import coeffvars, degree, is_of_func_type
from schubmult.schub_lib.positivity import compute_positive_rep
from schubmult.symbolic import Add, Mul, Pow, S, expand

# symengine.var("y_(1:100)")
# symengine.var("z_(1:100)")
# bargain = E(1, 1, y_4, z_1)*E(1, 1, y_4, z_2)*E(1, 3, y_1, y_4, y_5, z_1, z_2, z_3) + E(1, 1, y_4, z_1)*E(1, 1, y_5, z_4)*E(1, 3, y_1, y_4, y_5, z_1, z_2, z_3) - E(1, 1, y_4, z_1)*E(2, 3, y_1, y_4, y_5, z_1, z_2) + E(1, 1, y_5, z_2)*E(1, 1, y_5, z_4)*E(1, 3, y_1, y_4, y_5, z_1, z_2, z_3) - E(1, 1, y_5, z_4)*E(2, 3, y_1, y_4, y_5, z_1, z_2) + E(3, 3, y_1, y_4, y_5, z_1)
sympy.init_printing(pretty_print=False)

z_ring = SingleSchubertRing(z)

bagel1 = [4, 1, 5, 2, 3]
porn = [5, 1, 4, 2, 3]
for k, bargain in (DSx(bagel1, elem_sym=True) * DSx(porn, "z")).items():
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
        
        print(f"{k}: {dct}")
        schuber = z_ring.from_expr(Add(*[expand(efficient_subs(vl,{y[i]: S.Zero for i in range(10)}),func=True,mul=False) for vl in dct.values()]))
        print(schuber)
        for monom, bargle in dct.items():
            if expand(bargle, func=True, deep=True, mul=True) != 0:
                try:
                    compute_positive_rep(expand(bargle, func=True, mul=False), y, z, False, False)
                    print(f"Yep: {monom}: {bargle=}")
                except Exception:
                    print(f"Nope {monom}: {bargle=}")

        
      # Anything that isn't Schub will disappear
