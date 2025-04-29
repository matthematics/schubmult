from schubmult.schub_lib.quantum_double import divide_out_diff
from schubmult import *
from symengine import *

# pol = sympify("(y_3 + z_1)*((y_1 + z_1)*(y_8 + z_1) + (y_3 + z_2)*(y_1 + y_8 + z_1 + z_2)) - (y_1 + z_1)*(y_3 + z_1)*(y_8 + z_1)")
# pol = sympify("(y_3 + z_1)*(y_1 + y_3 + z_1 + z_2) - (y_1 + z_1)*(y_3 + z_1)")
pol = sympify("(y_3 + z_1)*((y_1 + z_1)*(y_8 + z_1) + (y_3 + z_2)*(y_1 + y_8 + z_1 + z_2)) + (y_8 + z_3)*((y_1 + z_1)*(y_8 + z_1) + (y_3 + z_2)*(y_1 + y_8 + z_1 + z_2)) - (y_1 + z_1)*(y_3 + z_1)*(y_8 + z_1)")

def proc_pos(poly):
    print(f"bingdang {poly=}")
    if str(poly).find("*y") != -1:
        raise Exception("bong")
    if str(poly).find("*z") != -1:
        raise Exception("bing")
    if str(poly).find("- y") != -1:
        raise Exception("bong")
    if str(poly).find("- z") != -1:
        raise Exception("bing")
    if str(poly).find("-") == -1:
        print("bongaofsin")
        return poly
    if isinstance(poly, Add):
        for arg in poly.args:
            for arg2 in arg.args:
                print(f"{arg2=}")
                if isinstance(arg2, Add) and len(arg2.args) == 2 and all(a in arg2.free_symbols for a in arg2.args):
                    print(f"divibble {arg2=}")
                    try:
                        beffo = proc_pos(divide_out_diff(poly,-arg2.args[1],arg2.args[0]))
                        return beffo
                    except Exception:
                        print("Welp that didn't work")
                        pass                                    
        return Add(*[proc_pos(arg) for arg in poly.args])
    #         print(f"{arg=}")
    #         if arg.args[0] == -1:
    #             for arg2 in arg.args[1:]:
    #                 print(f"{arg2=}")
    #                 if isinstance(arg2, Add) and len(arg2.args) == 2 and all(a in arg2.free_symbols for a in arg2.args):
    #                     print(f"divibble {arg2=}")
    #                     try:
    #                         beffo = proc_pos(divide_out_diff(poly,-arg2.args[1],arg2.args[0]))
    #                     except Exception:
    #                         print("Welp that didn't work")
    #                         pass                                    
    #         try:
    #             new_args += [proc_pos(arg)]
    #         except Exception:
    #             new_args += [arg]
    #     print(f"boing f{poly=}")
    #     return proc_pos(Add(*new_args))
    if isinstance(poly, Mul):
        print("bizzle")
        if poly.args[0] == -1:
            raise Exception("boingdoopy")
        return Mul(*[proc_pos(a) for a in poly.args])
    raise Exception("bongdiddly")
dang = proc_pos(pol)
print(f"{dang=}")
y = GeneratingSet("y")
z = GeneratingSet("z")
pol2 = pol.args[0]*divide_out_diff(pol.args[1],-z[1],y[1])
print(pol)
print(pol2)
print(expand(pol-dang))
# pol3 = divide_out_diff(pol2.args[0],-z[2],y[3])
# print(pol3)
