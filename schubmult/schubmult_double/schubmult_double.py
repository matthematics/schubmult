from symengine import *
from functools import cache
from itertools import chain
from schubmult.perm_lib import *
from schubmult.schubmult_yz import schubmult
import sys

n = 100

var = symarray('x', n)
var2 = symarray('y',n)
var3 = var2
var_r = symarray('r',n)

def trimcode(perm):
	cd = code(perm)
	while len(cd)>0 and cd[-1] == 0:
		cd.pop()
	return cd

subs_dict = {}

for i in range(1,n):
	sm = var_r[0]
	for j in range(1,i):
		sm += var_r[j]
	subs_dict[var2[i]] = sm


def main():
	perms=[]
	curperm = []
	
	pr = True
	ascode = False
	
	try:
		for s in sys.argv[1:]:
			if s == "-np":
				pr = False
				continue
			if s == "-code":
				ascode = True
				continue
			if s == "-":
				perms += [curperm]
				curperm = []
				continue
			curperm += [int(s)]
	except Exception:
		print("Usage: schubmult_double <-np> <perm1> - <perm2>")
		exit(1)
	
	perms += [curperm]
	
	if ascode:
		for i in range(len(perms)):
			perms[i] = uncode(perms[i])
	
	coeff_dict = {tuple(permtrim([*perms[0]])): 1}
	
	for perm in perms[1:]:
		coeff_dict = schubmult(coeff_dict,tuple(permtrim([*perm])),var2,var2)
		
	if pr:
		if ascode:
			width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys()])
		else:
			width = max([len(str(perm)) for perm in coeff_dict.keys()])
		
		coeff_perms = list(coeff_dict.keys())
		coeff_perms.sort(key=lambda x: (inv(x),*x))
		
		for perm in coeff_perms:
			val = sympify(coeff_dict[perm]).subs(subs_dict).expand()
			if val != 0:
				if ascode:					
					print(f"{str(trimcode(perm)):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
				else:
					print(f"{str(perm):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
			
if __name__ == "__main__":
	main()