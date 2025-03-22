from symengine import symarray

n = 100
fvar = 0

var = tuple(symarray("x", n).tolist())
var2 = tuple(symarray("y", n).tolist())
var3 = tuple(symarray("z", n).tolist())

var_x = symarray("x", 100).tolist()
var_y = var2
var_z = var3

x = var_x
y = var_y
z = var_z

var_r = symarray("r", 100)