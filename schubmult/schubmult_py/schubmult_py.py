import sys
from functools import cache
from itertools import chain
from schubmult.perm_lib import *
from symengine import *

var_x = symarray("x",100).tolist()

def single_variable(coeff_dict,varnum):
	ret = {}
	for u in coeff_dict:
		new_perms_k = elem_sym_perms(u,1,varnum)
		new_perms_km1 = []
		if varnum > 1:
			new_perms_km1 = elem_sym_perms(u,1,varnum-1)
		for perm, udiff in new_perms_k:
			if udiff == 1:
				ret[perm] = ret.get(perm,0) + coeff_dict[u]
		for perm, udiff in new_perms_km1:
			if udiff == 1:
				ret[perm] = ret.get(perm,0) - coeff_dict[u]
	return ret

def mult_poly(coeff_dict,poly):
	if poly in var_x:
		return single_variable(coeff_dict,var_x.index(poly))
	elif isinstance(poly,Mul):
		ret = coeff_dict
		for a in poly.args:
			ret = mult_poly(ret,a)
		return ret
	elif isinstance(poly,Pow):
		base = poly.args[0]
		exponent = int(poly.args[1])
		ret = coeff_dict
		for i in range(int(exponent)):
			ret = mult_poly(ret,base)
		return ret
	elif isinstance(poly,Add):
		ret = {}
		for a in poly.args:
			ret = add_perm_dict(ret,mult_poly(coeff_dict,a))
		return ret
	else:
		ret = {}
		for perm in coeff_dict:
			ret[perm] = poly*coeff_dict[perm]
		return ret

def schubmult(perm_dict,v):
	vn1 = inverse(v)
	th = theta(vn1)
	if th[0]==0:
		return perm_dict		
	mu = permtrim(uncode(th))
	vmu = permtrim(mulperm(list(v),mu))
	inv_vmu = inv(vmu)
	inv_mu = inv(mu)
	ret_dict = {}
	vpaths = [([(vmu,0)],1)]
	while th[-1] == 0:
		th.pop()
	vpathdicts = compute_vpathdicts(th,vmu)
	
	mx_th = [0 for i in range(len(th))]
	for index in range(len(th)):
		for vp in vpathdicts[index]:
			for v2,vdiff,s in vpathdicts[index][vp]:
				if th[index]-vdiff > mx_th[index]:
					mx_th[index] = th[index] - vdiff
	
	for u,val in perm_dict.items():
		inv_u = inv(u)		
		vpathsums = {u: {(1,2): val}}		
		
		for index in range(len(th)):
			newpathsums = {}
			for up in vpathsums:
				inv_up = inv(up)
				newperms = elem_sym_perms(up,min(mx_th[index],inv_mu-inv_vmu-(inv_up-inv_u)),th[index])
				for vp in vpathsums[up]:
					sumval = vpathsums[up][vp]
					if sumval == 0:
						continue
					for v2,vdiff,s in vpathdicts[index][vp]:
						addsumval = s*sumval
						for up2, udiff in newperms:					
							if vdiff + udiff == th[index]:								
								if up2 not in newpathsums:
									newpathsums[up2]={}
								newpathsums[up2][v2] = newpathsums[up2].get(v2,0)+addsumval							
			vpathsums = newpathsums
		toget = tuple(vmu)
		ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget,0) for ep in vpathsums},ret_dict)
	return ret_dict

def main():
	try:
		perms=[]
		curperm = []
		
		pr = True
		coprod = False
		ascode = False
		mult = False
		mulstring = ""
		
		try:
			for s in sys.argv[1:]:
				if mult:
					mulstring += s
					continue
				if s == "-mult":
					mult = True
					continue
				if s == "-np" or s == "--no-print":
					pr = False
					continue
				if s == "-coprod":
					coprod = True
					continue
				if s == "-code":
					ascode = True
					continue
				if s == "-":
					perms += [tuple(curperm)]
					curperm = []
					continue				
				else:
					curperm += [int(s)]
		except Exception:
			print("**** schubmult_py ****")
			print("Purpose: Compute products (and coproducts) of ordinary Schubert polynomials")
			print("Usage: schubmult_py <-np|--no-print> <-code> perm1 - perm2 <- perm3...>")
			print("Alternative usage: schubmult_py -coprod <-np> <perm> - <index list>")
			exit(1)
		
		perms += [tuple(curperm)]
		
		if coprod:
			if ascode:
				perms[0] = tuple(permtrim(uncode(perms[0])))
			pos = [*perms[1]]
			pos.sort()
			mperm = perms[0]
	
			cd = code(mperm)
			perms[0] = mperm
	
			while cd[-1] == 0:
				cd.pop()
			k = len(pos)
			n = len(perms[0])
			kcd = [pos[i]-i-1 for i in range(len(pos))] + [n+1-k for i in range(k,n)]
			N = len(kcd)
			kperm = inverse(uncode(kcd))
			coeff_dict = {tuple(permtrim(kperm)): 1}
			coeff_dict = schubmult(coeff_dict,tuple(permtrim([*perms[0]])))
		
			inv_perm0 = inv(perms[0])
			inv_kperm = inv(kperm)
			inverse_kperm = inverse(kperm)
			if pr:
				for perm, val in coeff_dict.items():
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
						if val != 0:
							#firstcode = code(firstperm)
							#while len(firstcode)>0 and firstcode[-1] == 0:
							#	firstcode.pop()
							#secondcode = code(secondperm)
							#while len(secondcode)>0 and secondcode[-1] == 0:
							#	secondcode.pop()
							if ascode:
								print(f"{val} {trimcode(firstperm)} {trimcode(secondperm)}")
							else:
								print(f"{val} {tuple(permtrim(firstperm))} {tuple(permtrim(secondperm))}")
		else:
			if ascode:
				for i in range(len(perms)):
					perms[i] = tuple(permtrim(uncode(perms[i])))
		
		
			perms.sort(reverse=True,key=lambda x: sum(theta(inverse(x)))-inv(x))
					
			coeff_dict = {tuple(permtrim([*perms[0]])): 1}
					
			for perm in perms[1:]:
				coeff_dict = schubmult(coeff_dict,tuple(permtrim([*perm])))
			if mult:
				mul_exp = sympify(mulstring)
				coeff_dict = mult_poly(coeff_dict,mul_exp)
				
			if pr:
				for perm, val in coeff_dict.items():
					if val!= 0:
						if ascode:						
							print(f"{val}  {trimcode(perm)}")
						else:
							print(f"{val}  {perm}")					
	except BrokenPipeError:
		pass

if __name__ == "__main__":
	main()