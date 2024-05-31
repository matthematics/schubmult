from schubmult.perm_lib import *
from schubmult.schubmult_yz import compute_positive_rep, posify
from symengine import *
import sys

var2 = symarray("y",100)
var3 = symarray("z",100)

var_y = var2.tolist()
var_z = var3.tolist()
var_x = symarray("x",100).tolist()

x = var_x
y = var_y
z = var_z

def E(p,k,varl=var_y[1:]):
	return elem_sym_poly_q(p,k,var_x[1:],varl)

def single_variable(coeff_dict,varnum):
	ret = {}
	for u in coeff_dict:
		if varnum -1 < len(u):
			ret[u] = ret.get(u,0) + var2[u[varnum-1]]*coeff_dict[u]
		else:
			ret[u] = ret.get(u,0) + var2[varnum]*coeff_dict[u]
		new_perms_k = elem_sym_perms_q(u,1,varnum)
		new_perms_km1 = []
		if varnum > 1:
			new_perms_km1 = elem_sym_perms_q(u,1,varnum-1)
		for perm, udiff, mul_val in new_perms_k:
			if udiff == 1:
				ret[perm] = ret.get(perm,0) + coeff_dict[u]*mul_val
		for perm, udiff, mul_val in new_perms_km1:
			if udiff == 1:
				ret[perm] = ret.get(perm,0) - coeff_dict[u]*mul_val
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
		

def nil_hecke(perm_dict,v,n,var2=var2,var3=var3):
	if v == (1,2):
		return perm_dict
	th = strict_theta(inverse(v))
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
				newperms = elem_sym_perms_q_op(up,mx_th,th[index],n)
				for up2, udiff,mul_val in newperms:
					if up2 not in newpathsums:
						newpathsums[up2]={}
					for v in vpathdicts[index]:
						sumval = vpathsums[up].get(v,zero)*mul_val
						if sumval == 0:
							continue
						for v2,vdiff,s in vpathdicts[index][v]:
							#print(f"{code(up2)=} {elem_sym_func_q(th[index],index+1,up,up2,v,v2,udiff,vdiff,var2,var3)=} {mul_val=} {sumval=}")							
							newpathsums[up2][v2] = newpathsums[up2].get(v2,zero)+s*sumval*elem_sym_func_q(th[index],index+1,up2,up,v,v2,udiff,vdiff,var2,var3)
			vpathsums = newpathsums
		toget = tuple(vmu)
		ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget,0) for ep in vpathsums},ret_dict)
	return ret_dict


	
def schubmult(perm_dict,v,var2=var2,var3=var3):
	if v == (1,2):
		return perm_dict
	th = strict_theta(inverse(v))
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
				newperms = elem_sym_perms_q(up,min(mx_th,(inv_mu-(inv_up-inv_u))-inv_vmu),th[index])
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

def schubmult_db(perm_dict,v,var2=var2,var3=var3):
	if v == (1,2):
		return perm_dict
	th = medium_theta(inverse(v))
	#print(f"{th=}")
	while th[-1] == 0:
		th.pop()
	#if len(set(th))!=len(th):
	#	print(f"medium theta {th=}")
	mu = permtrim(uncode(th))
	vmu = permtrim(mulperm([*v],mu))
	inv_vmu = inv(vmu)
	inv_mu = inv(mu)
	ret_dict = {}
	vpaths = [([(vmu,0)],1)]
	
	thL = len(th)
	#if thL!=2 and len(set(thL))!=1:
	#	raise ValueError("Not what I can do")
	vpathdicts = compute_vpathdicts(th,vmu,True)
	#print(f"{vpathdicts=}")
	for u,val in perm_dict.items():
		inv_u = inv(u)
		vpathsums = {u: {(1,2): val}}
		for index in range(thL):
			if index>0 and th[index-1] == th[index]:
				continue				
			mx_th = 0
			for vp in vpathdicts[index]:
				for v2,vdiff,s in vpathdicts[index][vp]:
					if th[index]-vdiff > mx_th:
						mx_th = th[index] - vdiff
			if index<len(th)-1 and th[index] == th[index+1]:
				mx_th1 = 0
				for vp in vpathdicts[index+1]:
					for v2,vdiff,s in vpathdicts[index+1][vp]:
						if th[index+1]-vdiff > mx_th1:
							mx_th1 = th[index+1] - vdiff				
				newpathsums = {}
				for up in vpathsums:
					newpathsums0 = {}
					inv_up = inv(up)
					newperms = double_elem_sym_q(up,mx_th,mx_th1,th[index])
					#for up1, up2, udiff1,udiff2,mul_val1,mul_val2 in newperms:
					for v in vpathdicts[index]:
						sumval = vpathsums[up].get(v,zero)
						if sumval == 0:
							continue
						for v2,vdiff2,s2 in vpathdicts[index][v]:
							for up1, udiff1, mul_val1 in newperms:
								esim1 = elem_sym_func_q(th[index],index+1,up,up1,v,v2,udiff1,vdiff2,var2,var3)*mul_val1*s2
								mulfac = sumval*esim1
								if (up1,udiff1,mul_val1) not in newpathsums0:
									newpathsums0[(up1,udiff1,mul_val1)] = {}
								#newpathsums0[(up1, udiff1, mul_val1
								newpathsums0[(up1,udiff1,mul_val1)][v2] = newpathsums0[(up1,udiff1,mul_val1)].get(v2,0) + mulfac
					
					for up1, udiff1, mul_val1 in newpathsums0:
						for v in vpathdicts[index+1]:
							sumval = newpathsums0[(up1,udiff1,mul_val1)].get(v,zero)
							if sumval == 0:
								continue
							for v2,vdiff2,s2 in vpathdicts[index+1][v]:
								for up2, udiff2, mul_val2 in newperms[(up1,udiff1,mul_val1)]:
									esim1 = elem_sym_func_q(th[index+1],index+2,up1,up2,v,v2,udiff2,vdiff2,var2,var3)*mul_val2*s2
									mulfac = sumval*esim1
									if up2 not in newpathsums:
										newpathsums[up2] = {}
									newpathsums[up2][v2] = newpathsums[up2].get(v2,0) + mulfac
											#for up2, udiff2, mul_val2 in newperms[(up1,udiff1,mul_val1)]:
											#	if up2 not in newpathsums:
											#		newpathsums[up2]={}
											#	for v3,vdiff3,s3 in vpathdicts[index+1][v2]:
											#			newpathsums[up2][v3] = newpathsums[up2].get(v3,zero)+s3*mul_val2*mulfac*elem_sym_func_q(th[index+1],index+2,up1,up2,v2,v3,udiff2,vdiff3,var2,var3)
			else:
				newpathsums = {}
				for up in vpathsums:
					inv_up = inv(up)
					newperms = elem_sym_perms_q(up,min(mx_th,(inv_mu-(inv_up-inv_u))-inv_vmu),th[index])
					for up2, udiff,mul_val in newperms:
						if up2 not in newpathsums:
							newpathsums[up2]={}
						for v in vpathdicts[index]:
							sumval = vpathsums[up].get(v,zero)*mul_val
							if sumval == 0:
								continue
							for v2,vdiff,s in vpathdicts[index][v]:
								newpathsums[up2][v2] = newpathsums[up2].get(v2,zero)+s*sumval*elem_sym_func_q(th[index],index+1,up,up2,v,v2,udiff,vdiff,var2,var3)
			vpathsums = newpathsums
		toget = tuple(vmu)
		ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget,0) for ep in vpathsums},ret_dict)
	return ret_dict

q_var2 = q_var.tolist()

def sum_q_dict(q_dict1,q_dict2):
	ret = {**q_dict1}
	for key in q_dict2:
		ret[key] = ret.get(key,0) + q_dict2[key]
	return ret

def mul_q_dict(q_dict1,q_dict2):
	ret = {}
	for key1 in q_dict1:
		for key2 in q_dict2:
			key3 = key1*key2
			ret[key3] = ret.get(key3,0) + q_dict1[key1]*q_dict2[key2]
	return ret

def factor_out_q_keep_factored(poly):
	ret = {}
	if str(poly).find("q") == -1:
		ret[1] = poly
		return ret
	elif poly in q_var2:
		ret[poly] = 1
		return ret
	elif isinstance(poly,Add):
		ag = poly.args
		ret = factor_out_q_keep_factored(ag[0])
		for i in range(1,len(ag)):
			ret = sum_q_dict(ret,factor_out_q_keep_factored(ag[i]))
		return ret
	elif isinstance(poly,Mul):
		ag = poly.args
		ret = factor_out_q_keep_factored(ag[0])
		for i in range(1,len(ag)):
			ret = mul_q_dict(ret,factor_out_q_keep_factored(ag[i]))
		return ret
	elif isinstance(poly,Pow):
		base = poly.args[0]
		exponent = int(poly.args[1])
		#print(f"exponent {exponent}")
		work_val = factor_out_q_keep_factored(base)
		ret = {1: 1}
		while exponent > 0:
			if exponent % 2 == 1:
				if ret == {1: 1}:
					ret = {**work_val}
				else:
					ret = mul_q_dict(ret,work_val)
				exponent -= 1
			else:
				work_val = mul_q_dict(work_val,work_val)
				exponent //= 2
		return ret
	return ret

def factor_out_q(poly):
	coeff_dict = expand(poly).as_coefficients_dict()
	ret = {}	
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

var2_t = tuple(var2.tolist())
var3_t = tuple(var3.tolist())	

def print_usage():
	print("**** schubmult_q_yz ****")
	print("Purpose: Compute Molev-Sagan coefficients of quantum double Schubert polynomials")
	print("Usage: schubmult_q_yz <-np|--no-print> <-code> <--display-positive> <--optimizer-message> perm1 - perm2 < - perm3 .. > <-mult poly_expression>")
	print("       *** Computes products")
	print("Alternative usage: schubmult_q_yz -nil-hecke n <-code> <--display-positive> perm")
	print("       *** Computes nil-Hecke representation of quantum double Schubert polynomial, limiting to the group S_n")

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
		nilhecke = False
		check = True
		msg = False
		just_nil = False
		mult = False
		slow = False
		
		nil_N = 0
		
		mulstring = ""
		norep = False
		expa = False
		
		try:
			for s in sys.argv[1:]:
				if just_nil:
					just_nil = False
					nil_N = int(s)
					continue
				if s == "--slow":
					slow = True
					continue
				if s == "--norep":
					norep = True
					continue
				if s == "--expand":
					expa = True
					continue
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
				if s == "-nil-hecke":
					nilhecke = True
					just_nil = True
					continue
				if s == "--version":
					print(f"Python version {sys.version}")
					exit(0)
				if s == "-code":
					ascode = True
					continue
				if s == "--usage" or s == "--help":
					print_usage()
					exit(0)				
				if s == "-":
					perms += [curperm]
					curperm = []
					continue				
				curperm += [int(s)]
		except Exception:
			print_usage()
			exit(1)
				
		perms += [curperm]
		
		
		if ascode:
			for i in range(len(perms)):
				perms[i] = tuple(permtrim(uncode(perms[i])))
		else:
			for i in range(len(perms)):
				if len(perms[i])<2 and (len(perms[i])==0 or perms[i][0]==1):
					perms[i] = (1,2)
				perms[i] = tuple(permtrim([*perms[i]]))
		
		size = 0
		L = len(perms)
		
		if nilhecke:
			coeff_dict = nil_hecke({(1,2): 1},perms[0],nil_N)			
			rep = ("y","x")		
		else:
			coeff_dict = {perms[0]: 1}
			for perm in perms[1:]:
				if not slow:
					coeff_dict = schubmult_db(coeff_dict,perm)
				else:
					coeff_dict = schubmult(coeff_dict,perm)
			if mult:
				for v in var2:
					globals()[str(v)] = v
				for v in var3:
					globals()[str(v)] = v	
				for v in var_x:
					globals()[str(v)] = v	
				for v in q_var:
					globals()[str(v)] = v
				q = q_var
				mul_exp = eval(mulstring)
				coeff_dict = mult_poly(coeff_dict,mul_exp)
			rep = ("","")			
		
		if pr:
			#if ascode:
			#	width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys() if expand(coeff_dict[perm])!=0])
			#else:
			#	width = max([len(str(perm)) for perm in coeff_dict.keys() if expand(coeff_dict[perm])!=0])
			
			coeff_perms = list(coeff_dict.keys())
			coeff_perms.sort(key=lambda x: (inv(x),*x))
			
			for perm in coeff_perms:
				val = coeff_dict[perm]
				if expand(val) != 0:
					notint = False
					try:
						int(val)
					except Exception:
						notint = True
						val2 = 0
						if display_positive:
							q_dict = factor_out_q_keep_factored(val)
							for q_part in q_dict:
								#print(f"{q_part=} {q_dict[q_part]=}")
								try:
									val2 += q_part*int(q_dict[q_part])
								except Exception:
									try:
										if len(perms) == 2 and q_part == 1 and not mult:
											u = permtrim([*perms[0]])
											v = permtrim([*perms[1]])
											val2 += posify(q_dict[q_part],tuple(u),tuple(v),perm,var2_t,var3_t,msg,False)
										elif len(perms) == 2 and q_part in q_var2 and not mult:
											i = q_var2.index(q_part)
											u = permtrim([*perms[0]])
											v = permtrim([*perms[1]])											
											#print(f"{u=} {v=} {q_part=} {q_dict[q_part]=}")
											if i<len(u) and i<len(v) and u[i-1]>u[i] and v[i-1]>v[i]:
												u[i], u[i-1] = u[i-1], u[i]
												v[i], v[i-1] = v[i-1], v[i]
												#print(f"new {u=} {v=}")
												val2 += q_part*posify(q_dict[q_part],tuple(permtrim(u)),tuple(permtrim(v)),perm,var2_t,var3_t,msg,False)
										else:
											val2 += q_part*compute_positive_rep(q_dict[q_part],var2_t,var3_t,msg,False)
									except Exception as e:
										if mult:
											print("warning; --display-positive is on but result is not positive",file=sys.stderr)
											val2 = val
											break
										else:
											print(f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {coeff_dict.get(perm,0)=}")
											print(f"Exception: {e}")
											exit(1)
							if check and expand(val - val2)!=0:
								if mult:
									val2 = val
								else:
									print(f"error: value not equal; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {coeff_dict.get(perm,0)=}")
									exit(1)
							val = val2
					if expa:
						val = expand(val)
					if val!=0:
						if ascode:
							if norep:
								print(f"{str(trimcode(perm))}  {str(val).replace(*rep)}")	
							else:
								print(f"{str(trimcode(perm))}  {str(val).replace('**','^').replace('*',' ').replace(*rep)}")	
						else:
							if norep:
								print(f"{str(perm)}  {str(val).replace(*rep)}")	
							else:
								print(f"{str(perm)}  {str(val).replace('**','^').replace('*',' ').replace(*rep)}")	
	except BrokenPipeError:
		pass
		
if __name__ == "__main__":
	main()