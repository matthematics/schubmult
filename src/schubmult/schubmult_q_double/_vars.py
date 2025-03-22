from symengine import symarray
from schubmult.perm_lib import q_var

var2 = symarray("y", 100)
var3 = symarray("z", 100)

var_y = var2.tolist()
var_z = var3.tolist()
var_x = symarray("x", 100).tolist()

x = var_x
y = var_y
z = var_z

q_var2 = q_var.tolist()

var2_t = tuple(var2.tolist())
var3_t = tuple(var3.tolist())

n = 100

