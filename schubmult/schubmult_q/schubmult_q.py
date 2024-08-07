from symengine import *
from functools import cache
from itertools import chain
from schubmult.perm_lib import *
from schubmult.schubmult_q_yz import factor_out_q_keep_factored
import schubmult.schubmult_yz as yz
import sys

var = symarray('x', n)
var2 = symarray('y',n)
var3 = var2
var_r = symarray('r',n)

var_q = Symbol("q")

subs_dict = {}

var_x = symarray("x",100).tolist()


def single_variable(coeff_dict,varnum):
	ret = {}
	for u in coeff_dict:
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

for i in range(1,n):
	sm = var_r[0]
	for j in range(1,i):
		sm += var_r[j]
	subs_dict[var2[i]] = sm

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
					for v in vpathdicts[index]:
						sumval = vpathsums[up].get(v,zero)
						if sumval == 0:
							continue
						for v2,vdiff2,s2 in vpathdicts[index][v]:
							for up1, udiff1, mul_val1 in newperms:
								if (up1,udiff1,mul_val1) not in newpathsums0:
									newpathsums0[(up1,udiff1,mul_val1)] = {}
								if udiff1 + vdiff2 == th[index]:
									newpathsums0[(up1,udiff1,mul_val1)][v2] = newpathsums0[(up1,udiff1,mul_val1)].get(v2,zero)+s2*sumval*mul_val1
					
					for up1, udiff1, mul_val1 in newpathsums0:
						for v in vpathdicts[index+1]:
							sumval = newpathsums0[(up1,udiff1,mul_val1)].get(v,zero)
							if sumval == 0:
								continue
							for v2,vdiff2,s2 in vpathdicts[index+1][v]:
								for up2, udiff2, mul_val2 in newperms[(up1,udiff1,mul_val1)]:
									if up2 not in newpathsums:
										newpathsums[up2]={}
									if udiff2 + vdiff2 == th[index+1]:
										newpathsums[up2][v2] = newpathsums[up2].get(v2,zero)+s2*sumval*mul_val2
			else:
				newpathsums = {}
				for up in vpathsums:
					inv_up = inv(up)
					newperms = elem_sym_perms_q(up,min(mx_th,(inv_mu-(inv_up-inv_u))-inv_vmu),th[index])
					for up2, udiff, mul_val in newperms:
						if up2 not in newpathsums:
							newpathsums[up2]={}
						for v in vpathdicts[index]:
							sumval = vpathsums[up].get(v,zero)
							if sumval == 0:
								continue
							for v2,vdiff,s in vpathdicts[index][v]:
								if udiff+vdiff==th[index]:
									newpathsums[up2][v2] = newpathsums[up2].get(v2,zero)+s*sumval*mul_val
			vpathsums = newpathsums
		toget = tuple(vmu)
		ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget,0) for ep in vpathsums},ret_dict)
	return ret_dict


def schubmult(perm_dict,v):
	th = strict_theta(inverse(v))
	mu = permtrim(uncode(th))
	vmu = permtrim(mulperm([*v],mu))
	#print(f"{th=} {mu=} {vmu=}")
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
				for up2, udiff, mul_val in newperms:
					if up2 not in newpathsums:
						newpathsums[up2]={}
					for v in vpathdicts[index]:
						sumval = vpathsums[up].get(v,zero)
						if sumval == 0:
							continue
						for v2,vdiff,s in vpathdicts[index][v]:
							if udiff+vdiff==th[index]:
								newpathsums[up2][v2] = newpathsums[up2].get(v2,zero)+s*sumval*mul_val
			vpathsums = newpathsums
		toget = tuple(vmu)
		ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget,0) for ep in vpathsums},ret_dict)
	return ret_dict

var_q = Symbol("q")

def grass_q_replace(perm,k,d,n):
	if k-d<0:
		return None
	#print("Here")
	ret = []
	cd = code(perm)
	for i in range(k-d,k):
		if i>=len(cd) or cd[i]<d:
			return None
	grass_rep = [0 for i in range(n)]
	perm2 = [*perm] + [i+1 for i in range(len(perm),n)]
	for i in range(k,n):
		grass_rep[perm2[i]-1] = 2
	num_0 = 0
	#print(f"{grass_rep=} {d=}")
	for i in range(len(grass_rep)-1,-1,-1):
		if num_0 == d:
			break
		if grass_rep[i] == 0:
			grass_rep[i] = 1
			num_0 += 1
	num_2 = 0
	for i in range(len(grass_rep)):
		if num_2 == d:
			break
		if grass_rep[i] == 2:
			grass_rep[i] = 1
			num_2 += 1
	#print(f"New {grass_rep=}")
	k1 = k - d
	k2 = k + d
	pos_1 = 0
	pos_2 = 0
	pos_3 = 0
	new_perm = [0 for i in range(n)]
	for i in range(len(grass_rep)):
		if grass_rep[i] == 0:
			new_perm[pos_1] = i+1
			pos_1 += 1
		if grass_rep[i] == 1:
			new_perm[k1+pos_2] = i+1
			pos_2 += 1
		if grass_rep[i] == 2:
			new_perm[k2+pos_3] = i+1
			pos_3 += 1
	return tuple(permtrim(new_perm))

def to_two_step(perm,k1,k2,n):
	rep = [0 for i in range(n)]
	perm2 = [*perm] + [i+1 for i in range(len(perm),n)]	
	for i in range(n):
		if i<k1:
			rep[perm2[i]-1] = 0
		elif i<k2:
			rep[perm2[i]-1] = 1
		else:
			rep[perm2[i]-1] = 2
	return rep

def main():
	try:			
		perms=[]
		curperm = []
		
		pr = True
		ascode = False
		grass = False
		grass_q_n = 0
		equiv = False
		mult = False
		mulstring = ""
		slow = False
		parabolic = False
		just_parabolic = False
		parabolic_index = []
		
		try:
			for s in sys.argv[1:]:
				if s == "-np" or s == "--no-print":
					pr = False
					continue
				if s == "--slow":
					slow = True
					continue
				if s == "-parabolic":
					just_parabolic = True
					parabolic = True
					continue
				if just_parabolic:
					just_parabolic = False
					parabolic_index = [int(i) for i in s.split(",")]
					parabolic_index.sort()
					continue
				if mult:
					mulstring += s
					continue
				if s == "-code":
					ascode = True
					continue
				if s == "-grass":
					grass = None
					continue
				if s == "-equiv":
					equiv = True
					continue
				if grass is None:
					grass = True
					grass_q_n = int(s)
					continue
				if s == "-mult":
					mult = True
					continue
				if s == "-":
					perms += [curperm]
					curperm = []
					continue
				curperm += [int(s)]
		except Exception:
			print("**** schubmult_q ****")
			print("Purpose: Compute the structure constants of quantum Schubert polynomials")
			print("Usage: schubmult_q <-np|--no-print> <-parabolic descent1,descent2,...> <-code> <-grass n> <-equiv> perm1 - perm2 < - perm 3 ... >")
			#print("For the -grass option, must use Grassmannian permutations. -equiv only works together with -grass.")			
			exit(1)
		
		perms += [curperm]
		
		if parabolic and len(perms) != 2:
			print("Only two permutations supported for parabolic.")
			exit(1)
							
		
		if grass:
			perms_t = []
			if ascode:
				perms_t = [tuple(permtrim(uncode(perms[i]))) for i in range(len(perms))]
			else:
				perms_t = [tuple(permtrim(perms[i])) for i in range(len(perms))]
			
			k = -1
			
			for perm in perms_t:
				desc = -1
				for i in range(len(perm)-1):
					if desc != -1 and perm[i]>perm[i+1]:
						print("Error: permutations must have one descent")
						exit(1)
					elif desc == -1:
						if perm[i]>perm[i+1]:
							if k != -1 and k != i+1:
								print("Error: permutations must all have the same descent")
								exit(1)
							k = i+1
							desc = i+1
			#perms1 = [perms_t[0]]
			#perms2 = [perms_t[1]]
			perms1 = []
			perms2 = []
			for d in range(100,-1,-1):
				#print(f"{k=} {d=}")
				grass_rep_1 = grass_q_replace(perms_t[0],k,d,grass_q_n)
				grass_rep_2 = grass_q_replace(perms_t[1],k,d,grass_q_n)
				#print(f"{grass_rep_1} {grass_rep_2}")
				if grass_rep_1 is None or grass_rep_2 is None:
					continue
				perms1 += [grass_rep_1]
				perms2 += [grass_rep_2]
			for i in range(len(perms1)):
				d = len(perms1) - i - 1
				#print(f"{d=} {perms1[i]=} {perms2[i]=}")
				if equiv:
					coeff_dict = yz.schubmult({perms1[i]: 1},perms2[i],var2,var2)
				else:
					coeff_dict = yz.schubmult({perms1[i]: 1},perms2[i],[0 for i in range(100)],[0 for i in range(100)])
				
				#if ascode:
				#	width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys()])
				#else:
				#	width = max([len(str(perm)) for perm in coeff_dict.keys()])
				k1 = k - d
				k2 = k + d
				coeff_perms = list([key for key in coeff_dict.keys() if len(key)<=grass_q_n])
				
				coeff_perms.sort(key=lambda x: (inv(x),*x))
				
				coeff_perms2 = []
				
				for perm in coeff_perms:
					two_step = to_two_step(perm,k1,k2,grass_q_n)					
					num_1 = len([i for i in range(len(two_step)) if two_step[i] == 1])
					if num_1 != 2*d:
						coeff_perms2 += [None]
						continue
					one_step = [*two_step]
					one_step.reverse()
					two_step_r = [*two_step]
					two_step_r.reverse()
					for i in range(len(one_step)-1,-1,-1):
						if one_step[i] == 0:
							break
						if one_step[i] == 1:
							one_step[i] = 0
					no_good = False
					for i in range(len(one_step)):
						if one_step[i] == 1:
							one_step[i] = 2
						elif one_step[i] == 2:
							if len([j for j in range(len(one_step)) if one_step[j]==1]) != 0:
								no_good = True
								break
					num_0 = len([i for i in range(len(one_step)) if one_step[i] == 0])
					if num_0 != k or no_good:
						#print(f"No good {one_step=} {two_step_r=} {no_good=}")
						coeff_perms2 += [None]
						continue
					pos_0 = 0
					pos_1 = 0
					one_step.reverse()
					#print(f"{two_step=}")
					#print(f"{one_step=}")
					grass_perm = [0 for i in range(len(one_step))]
					for i in range(len(one_step)):
						if one_step[i] == 0:
							grass_perm[pos_0] = i+1
							pos_0 += 1
						else:
							grass_perm[k+pos_1] = i+1
							pos_1 += 1
					coeff_perms2 += [tuple(permtrim(grass_perm))]
				
				try:
					if ascode:
						width = max([len(str(trimcode(perm))) for perm in coeff_perms2 if perm is not None])
					else:
						width = max([len(str(perm)) for perm in coeff_perms2 if perm is not None])
				except ValueError:
					continue
				
				for i in range(len(coeff_perms)):
					if coeff_perms2[i] is None:
						continue
					perm = coeff_perms[i]
					val = (var_q**d) * sympify(coeff_dict[perm]).subs(subs_dict).expand()
					
					if val != 0:
						if ascode:					
							print(f"{str(trimcode(coeff_perms2[i])):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
						else:
							print(f"{str(coeff_perms2[i]):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
						#print(f"{str(two_step):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
						#else:
						#	print(f"{str(perm):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
		else:		
			if ascode:
				for i in range(len(perms)):
					perms[i] = uncode(perms[i])
		
			if parabolic:
				for i in range(len(parabolic_index)):
					index = parabolic_index[i] - 1
					if sg(index,perms[0]) == 1 or sg(index,perms[1]) == 1:
						print("Parabolic given but elements are not minimal length coset representatives.")
						exit(1)
		
			coeff_dict = {tuple(permtrim([*perms[0]])): 1}
			
			if not slow:
				for perm in perms[1:]:
					coeff_dict = schubmult_db(coeff_dict,tuple(permtrim([*perm])))
			else:
				for perm in perms[1:]:
					coeff_dict = schubmult(coeff_dict,tuple(permtrim([*perm])))
			
			if mult:
				mul_exp = sympify(mulstring)
				coeff_dict = mult_poly(coeff_dict,mul_exp)
				
			if pr:						
				if ascode:
					width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys() if expand(coeff_dict[perm])!=0])
				else:
					width = max([len(str(perm)) for perm in coeff_dict.keys() if expand(coeff_dict[perm])!=0])
				
				if parabolic:
					w_P = longest_element(parabolic_index)
					w_P_prime = [1,2]
					coeff_dict_update = {}
					for w_1 in coeff_dict:
						val = coeff_dict[w_1]
						q_dict = factor_out_q_keep_factored(val)
						for q_part in q_dict:
							qv = q_vector(q_part)
							w = [*w_1]
							good = True
							parabolic_index2 = []
							for i in range(len(parabolic_index)):							
								if omega(parabolic_index[i],qv) == 0:
									parabolic_index2 += [parabolic_index[i]]
								elif omega(parabolic_index[i],qv) != -1:
									good = False
									break
							if not good:
								continue
							w_P_prime = longest_element(parabolic_index2)
							if not check_blocks(qv,parabolic_index):
								continue
							w = permtrim(mulperm(mulperm(w,w_P_prime),w_P))
							if not is_parabolic(w,parabolic_index):							
								continue
					
							w = tuple(permtrim(w))
							
							new_q_part = np.prod([q_var[index+1-count_less_than(parabolic_index,index+1)]**qv[index] for index in range(len(qv)) if index+1 not in parabolic_index])
							
							try:
								new_q_part = int(new_q_part)
							except Exception:
								pass
							q_val_part = q_dict[q_part]
							coeff_dict_update[w] = coeff_dict_update.get(w,0) + new_q_part*q_val_part
					coeff_dict = coeff_dict_update
				
				coeff_perms = list(coeff_dict.keys())
				coeff_perms.sort(key=lambda x: (inv(x),*x))
				
				for perm in coeff_perms:
					val = sympify(coeff_dict[perm]).expand()
					if val != 0:
						if ascode:					
							print(f"{str(trimcode(perm))}  {str(val).replace('**','^').replace('*',' ')}")
						else:
							print(f"{str(perm)}  {str(val).replace('**','^').replace('*',' ')}")
	except BrokenPipeError:
		pass
			
if __name__ == "__main__":
	main()