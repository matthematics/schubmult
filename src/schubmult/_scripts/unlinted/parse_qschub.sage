import sys

QSym = QuasiSymmetricFunctions(ZZ)
QS = QSym.QS()
dI = QSym.dI()
F = QSym.F()
M = QSym.M()

lines = []

with open(sys.argv[1], "r") as f:
    lines = f.readlines()

for line in lines:
    name, val = line.split(" = ")
    val = eval(val)
    dingbat = QS(val)
    print(dingbat.is_symmetric())
    print(f"{name} = {dingbat}")