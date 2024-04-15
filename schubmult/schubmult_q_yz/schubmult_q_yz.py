from schubmult.perm_lib import *
from schubmult.schubmult_yz import compute_positive_rep
from symengine import *
import sys


#q_var = symarray("q",100)

var2 = symarray("y",100)
var3 = symarray("z",100)

	

	

	
def schubmult(perm_dict,v,var2=var2,var3=var3):
	vn1 = inverse(v)
	th = [len(v)-i for i in range(1,len(v))]
	mu = permtrim(uncode(th))
	vmu = permtrim(mulperm([*v],mu))
	inv_vmu = inv(vmu)
	inv_mu = inv(mu)
	ret_dict = {}
	vpaths = [([(vmu,0)],1)]
	while th[-1] == 0:
		th.pop()
	thL = len(th)
	vpathdicts = compute_vpathdicts(th,vmu,True)
	for u,val in perm_dict.items():
		inv_u = inv(u)
		vpathsums = {u: {(1,2): val}}
		for index in range(thL):			
			mx_th = 0
			for vp in vpathdicts[index]:
				for v2,vdiff,s in vpathdicts[index][vp]:
					if th[index]-vdiff > mx_th:
						mx_th = th[index] - vdiff					
			newpathsums = {}
			for up in vpathsums:
				inv_up = inv(up)
				newperms = elem_sym_perms_q(up,mx_th,th[index])
				for up2, udiff,mul_val in newperms:
					if up2 not in newpathsums:
						newpathsums[up2]={}
					for v in vpathdicts[index]:
						sumval = vpathsums[up].get(v,zero)*mul_val
						if sumval == 0:
							continue
						for v2,vdiff,s in vpathdicts[index][v]:
							#print(f"{code(up2)=} {elem_sym_func_q(th[index],index+1,up,up2,v,v2,udiff,vdiff,var2,var3)=} {mul_val=} {sumval=}")							
							newpathsums[up2][v2] = newpathsums[up2].get(v2,zero)+s*sumval*elem_sym_func_q(th[index],index+1,up,up2,v,v2,udiff,vdiff,var2,var3)
			vpathsums = newpathsums
		toget = tuple(vmu)
		ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget,0) for ep in vpathsums},ret_dict)
	return ret_dict

def factor_out_q(poly):
	coeff_dict = expand(poly).as_coefficients_dict()
	ret = {}
	q_var2 = q_var.tolist()
	for key in coeff_dict:
		coeff = coeff_dict[key]
		if coeff == 0:
			continue
		q_part = 1
		yz_part = coeff
		if isinstance(key,Mul):
			for var_maybe_pow in key.args:
				if isinstance(var_maybe_pow,Pow):
					real_var = var_maybe_pow.args[0]
					if real_var in q_var2:
						q_part*=var_maybe_pow
					else:
						yz_part*=var_maybe_pow
				else:
					real_var = var_maybe_pow
					if real_var in q_var2:
						q_part*=var_maybe_pow
					else:
						yz_part*=var_maybe_pow
		elif isinstance(key,Pow):
			real_var = key.args[0]
			if real_var in q_var2:
				q_part*=key
			else:
				yz_part*=key
		else:
			if key in q_var2:
				q_part *= key
			else:
				yz_part*=key
			
		ret[q_part] = ret.get(q_part,0) + yz_part
	return ret
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
		try:
			for s in sys.argv[1:]:
				if s == "-np" or s == "--no-print":
					pr = False
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
					print("**** schubmult_q_yz ****")
					print("Purpose: Compute Molev-Sagan coefficients of quantum double Schubert polynomials")
					print("Usage: schubmult_q_yz <-np|--no-print> <-code> <--display-positive> <--optimizer-message> perm1 - perm2 < - perm3 .. >")
					exit(0)
				if s == "-":
					perms += [curperm]
					curperm = []
					continue
				curperm += [int(s)]
		except Exception:
			print("**** schubmult_q_yz ****")
			print("Purpose: Compute Molev-Sagan coefficients of quantum double Schubert polynomials")
			print("Usage: schubmult_q_yz <-np|--no-print> <-code> <--display-positive> <--optimizer-message> perm1 - perm2 < - perm3 .. >")
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
			coeff_dict = schubmult(coeff_dict,perm)
		
		if pr:
			if ascode:
				width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys() if expand(coeff_dict[perm])!=0])
			else:
				width = max([len(str(perm)) for perm in coeff_dict.keys() if expand(coeff_dict[perm])!=0])
			
			coeff_perms = list(coeff_dict.keys())
			coeff_perms.sort(key=lambda x: (inv(x),*x))
			
			for perm in coeff_perms:
				val = sympify(coeff_dict[perm]).simplify()
				if val != 0:
					notint = False
					try:
						int(val)
					except Exception:
						notint = True
						val2 = 0
						if display_positive:
							q_dict = factor_out_q(val)
							for q_part in q_dict:
								#print(f"{q_part=} {q_dict[q_part]=}")
								try:
									val2 += q_part*int(q_dict[q_part])
								except Exception:
									try:
										val2 += q_part*compute_positive_rep(q_dict[q_part],var2,var3,msg,False)
									except TypeError:
										print(f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {check_coeff_dict.get(perm,0)=}")
										exit(1)
							if check and expand(val - val2)!=0:
								print(f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {check_coeff_dict.get(perm,0)=}")
								exit(1)
							val = val2
					if val!=0:
						if ascode:
							print(f"{str(trimcode(perm)):>{width}}  {str(val).replace('**','^').replace('*',' ')}")	
						else:
							print(f"{str(perm):>{width}}  {str(val).replace('**','^').replace('*',' ')}")	
	except BrokenPipeError:
		pass
		
if __name__ == "__main__":
	main()