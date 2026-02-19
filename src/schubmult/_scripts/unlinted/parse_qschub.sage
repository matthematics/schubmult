import sys

QSym = QuasiSymmetricFunctions(ZZ)
QS = QSym.QS()
dI = QSym.dI()
M = QSym.M()

lines = []

with open(sys.argv[1], "r") as f:
    lines = f.readlines()

for line in lines:
    name, val = line.split(" = ")
    val = eval(val)
    print(f"{name} = {dI(val)}")