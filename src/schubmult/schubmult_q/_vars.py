from symengine import symarray

n = 100

var = symarray("x", n)
var2 = symarray("y", n)
var3 = var2
var_r = symarray("r", n)

subs_dict = {}

var_x = symarray("x", 100).tolist()

for i in range(1, n):
    sm = var_r[0]
    for j in range(1, i):
        sm += var_r[j]
    subs_dict[var2[i]] = sm