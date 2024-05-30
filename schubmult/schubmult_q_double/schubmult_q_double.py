from schubmult.perm_lib import *
from schubmult.schubmult_q_yz import schubmult, mult_poly
from symengine import *
import sys


q_var = symarray("q",100)

var2 = symarray("y",100)
var3 = symarray("z",100)

var3 = var2
var_r = symarray('r',100)

subs_dict = {}

for i in range(1,100):
	sm = var2[1]
	for j in range(1,i):
		sm += var_r[j]
	subs_dict[var2[i]] = sm
	
	

def main():
	global var2
	try:
		sys.setrecursionlimit(1000000)
	
		perms=[]
		curperm = []
		
		pr = True
		display_positive = False
		ascode = False
		coprod = False
		check = True
		msg = False
		mult = False
		mulstring = ""
		try:
			for s in sys.argv[1:]:
				if s == "-np" or s == "--no-print":
					pr = False
					continue
				if mult:
					mulstring += s
					continue
				if s == "-mult":
					mult = True
					continue
				if s == "-nocheck":
					check = False
					continue
				if s == "--display-positive":
					display_positive = True
					continue
				if s == "--optimizer-message":
					msg = True
					continue
				if s == "--version":
					print(f"Python version {sys.version}")
					exit(0)
				if s == "-code":
					ascode = True
					continue
				if s == "--usage":
					print("**** schubmult_q_double ****")
					print("Purpose: Compute equivariant Gromov-Witten invariants, structure constants of quantum double Schubert polynomials")
					print("Usage: schubmult_q_double <-np|--no-print> <-code> <--display-positive> <--optimizer-message> perm1 - perm2 < - perm3 .. >")
					#print("Alternative usage: schubmult_yz <-code> <--display-positive> -coprod perm - indexlist")
					exit(0)
				if s == "-":
					perms += [curperm]
					curperm = []
					continue
				curperm += [int(s)]
		except Exception:
			print("**** schubmult_q_double ****")
			print("Purpose: Compute equivariant Gromov-Witten invariants, structure constants of quantum double Schubert polynomials")
			print("Usage: schubmult_q_double <-np|--no-print> <-code> <--display-positive> <--optimizer-message> perm1 - perm2 < - perm3 .. >")
			exit(1)
				
		perms += [curperm]
		
		
		if ascode:
			for i in range(len(perms)):
				perms[i] = tuple(permtrim(uncode(perms[i])))
		else:
			for i in range(len(perms)):
				perms[i] = tuple(permtrim([*perms[i]]))
		
		size = 0
		L = len(perms)
		
		coeff_dict = {perms[0]: 1}
		for perm in perms[1:]:
			coeff_dict = schubmult(coeff_dict,perm,var2,var2)
		if mult:
			mul_exp = sympify(mulstring)
			coeff_dict = mult_poly(coeff_dict,mul_exp)
		
		if pr:
			if ascode:
				width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys() if expand(sympify(coeff_dict[perm]).subs(subs_dict))!=0])
			else:
				width = max([len(str(perm)) for perm in coeff_dict.keys() if expand(sympify(coeff_dict[perm]).subs(subs_dict))!=0])
			
			coeff_perms = list(coeff_dict.keys())
			coeff_perms.sort(key=lambda x: (inv(x),*x))
			
			for perm in coeff_perms:
				val = expand(sympify(coeff_dict[perm]).subs(subs_dict)).simplify()
				if val != 0:
					notint = False
					try:
						int(val)
					except Exception:
						notint = True
					if val!=0:
						if ascode:
							print(f"{str(trimcode(perm)):>{width}}  {str(val).replace('**','^').replace('*',' ')}")	
						else:
							print(f"{str(perm):>{width}}  {str(val).replace('**','^').replace('*',' ')}")	
	except BrokenPipeError:
		pass
		
if __name__ == "__main__":
	main()