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
	
	try:
		for s in sys.argv[1:]:
			if s == "-np":
				pr = False
				continue
			if s == "-":
				perms += [tuple(permtrim(curperm))]
				curperm = []
				continue
			curperm += [int(s)]
	except Exception:
		print("Usage: python3 schubmult_double.py <-np> <perm1> - <perm2>")
		exit(1)
	
	perms += [tuple(permtrim(curperm))]
	coeff_dict = {perms[0]: 1}
	
	for perm in perms[1:]:
		coeff_dict = schubmult(coeff_dict,perm,var2,var2)
		
	if pr:
		width = max([len(str(perm)) for perm in coeff_dict.keys()])
		
		coeff_perms = list(coeff_dict.keys())
		coeff_perms.sort(key=lambda x: (inv(x),*x))
		
		for perm in coeff_perms:
			val = sympify(coeff_dict[perm]).subs(subs_dict).expand()
			if val != 0:
				print(f"{str(perm):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
			
if __name__ == "__main__":
	main()