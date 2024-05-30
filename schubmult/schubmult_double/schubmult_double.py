from symengine import *
from functools import cache
from itertools import chain
from schubmult.perm_lib import *
from schubmult.schubmult_yz import schubmult, mult_poly
import sys

n = 100

var = symarray('x', n)
var2 = symarray('y',n)
var3 = var2
var_r = symarray('r',n)

var_q = Symbol("q")

subs_dict = {}

for i in range(1,n):
	sm = var2[1]
	for j in range(1,i):
		sm += var_r[j]
	subs_dict[var2[i]] = sm



def main():
	try:			
		perms=[]
		curperm = []
		
		pr = True
		ascode = False
		coprod = False
		mult = False
		mulstring = ""
		try:
			for s in sys.argv[1:]:
				if s == "-np" or s == "--no-print":
					pr = False
					continue
				if mult:
					mulstring+=s
					continue
				if s == "-code":
					ascode = True
					continue
				if s == "-coprod":
					coprod = True
					continue
				if s == "-":
					perms += [curperm]
					curperm = []
					continue
				if s == "-mult":
					mult = True
					continue
				curperm += [int(s)]
		except Exception:
			print("**** schubmult_double ****")
			print("Purpose: Compute products (and coproducts) of double Schubert polynomials in the same set of variables")
			print("Usage: schubmult_double <-np|--no-print> <-code> perm1 - perm2 < - perm 3 ... > <-mult poly>")
			print("Alternative usage: schubmult_double <-code> -coprod perm - indexlist")
			exit(1)
		
		perms += [curperm]
		
		
		if coprod:
			subs_dict_coprod = {}
			if ascode:
				mperm = tuple(permtrim(uncode(perms[0])))
			else:
				mperm = tuple(permtrim(perms[0]))
	
			cd = code(mperm)
			perms[0] = mperm
			pos = perms[1]
	
			while cd[-1] == 0:
				cd.pop()
			k = len(pos)
			n = len(perms[0])
			kcd = [pos[i]-i-1 for i in range(len(pos))] + [n+1-k for i in range(k,n)]
			N = len(kcd)
					
			for i in range(1,100):
				if i<=N:
					subs_dict_coprod[var[i]] = var2[i]
				else:
					subs_dict_coprod[var[i]] = var2[i-N]
			
			kperm = inverse(uncode(kcd))
			coeff_dict = {tuple(kperm): 1}
			coeff_dict = schubmult(coeff_dict,perms[0],var,var2)
			
			inv_perm0 = inv(perms[0])
			inv_kperm = inv(kperm)
			inverse_kperm = inverse(kperm)
			if pr:
				coeff_perms = list(coeff_dict.keys())
				coeff_perms.sort(key=lambda x: (inv(x),*x))
			
				perm_pairs = []
				
				for perm in coeff_perms:
					downperm = mulperm(list(perm),inverse_kperm)
					if inv(downperm) == inv(perm) - inv_kperm:
						flag = True			
						for i in range(N):
							if downperm[i] > N:
								flag = False
								break
						if not flag:
							continue
						firstperm = downperm[0:N]
						secondperm = [downperm[i]-N for i in range(N,len(downperm))]
						perm_pairs += [[permtrim(firstperm),permtrim(secondperm)]]
			
				if ascode:
					width = max([len(str(trimcode(perm[0]))+" "+str(trimcode(perm[1]))) for perm in perm_pairs])
				else:
					width = max([len(str(perm[0])+" "+str(perm[1])) for perm in perm_pairs])
			
				for perm in coeff_perms:
					val = coeff_dict[perm]
					downperm = mulperm(list(perm),inverse_kperm)
					if inv(downperm) == inv(perm) - inv_kperm:
						flag = True			
						for i in range(N):
							if downperm[i] > N:
								flag = False
								break
						if not flag:
							continue
						firstperm = downperm[0:N]
						secondperm = [downperm[i]-N for i in range(N,len(downperm))]
						val = sympify(val).subs(subs_dict_coprod)
						val = sympify(val).subs(subs_dict).expand()				
						if val != 0:
							if not ascode:
								width2 = width - len(str(permtrim(firstperm))) - len(str(permtrim(secondperm)))
								print(f"{tuple(permtrim(firstperm))}{' ':>{width2}}{tuple(permtrim(secondperm))}  {str(val).replace('**','^').replace('*',' ')}")
							else:
								width2 = width - len(str(trimcode(firstperm))) - len(str(trimcode(secondperm)))
								print(f"{trimcode(firstperm)}{' ':>{width2}}{trimcode(secondperm)}  {str(val).replace('**','^').replace('*',' ')}")
		else:
			if ascode:
				for i in range(len(perms)):
					perms[i] = uncode(perms[i])
			
			coeff_dict = {tuple(permtrim([*perms[0]])): 1}
			
			for perm in perms[1:]:
				coeff_dict = schubmult(coeff_dict,tuple(permtrim([*perm])),var2,var2)
			if mult:
				mul_exp = sympify(mulstring)
				coeff_dict = mult_poly(coeff_dict,mul_exp)
			
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
	except BrokenPipeError:
		pass
			
if __name__ == "__main__":
	main()