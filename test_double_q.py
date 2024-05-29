from symengine import *
import sys
from schubmult.perm_lib import *
from schubmult.schubmult_q_yz import schubmult, var2, var3, factor_out_q_keep_factored
from schubmult.schubmult_yz import compute_positive_rep
from itertools import permutations
from joblib import Parallel, delayed, Memory
from functools import cache


n = int(sys.argv[1])

def identify_k_path(perm0,perm1,perm2,k):
	perms1 = [u for u,udiff,mul_val in elem_sym_perms_q(perm0,k,k)]
	if perm1 in perms1:
		perms2 = [u for u,udiff,mul_val in elem_sym_perms_q(perm1,k,k)]
		if perm2 in perms2:
			return True
	return False

def get_cycles(perm):
	cycle_set = []
	done_vals = set()
	for i in range(len(perm)):
		p = i + 1
		if perm[i] == p:
			continue
		if p in done_vals:
			continue
		cycle = []
		m = -1
		max_index = -1
		while p not in done_vals:
			cycle += [p]
			done_vals.add(p)
			if p>m:
				m = p
				max_index = len(cycle) - 1
			p = perm[p-1]
		cycle = tuple(cycle[max_index+1:] + cycle[:max_index+1])
		cycle_set += [cycle]
	return cycle_set
	
#up, down
def get_reflections(perm1,perm1i,perm2):
	cycle_agg = mulperm(perm1i,[*perm2])
	cycles = get_cycles(cycle_agg)
	up_ref = []
	for cycle in cycles:
		for i in range(len(cycle)-2,-1,-1):
			if getpermval(perm1,cycle[i]-1)<getpermval(perm1,cycle[i+1]-1):
				up_ref+=[(cycle[i],cycle[-1])]			
	return up_ref

def double_elem_sym_q(u,p1,p2,k):
	ret_list = {}
	perms1 = elem_sym_perms_q(u,p1,k)
	iu = inverse(u)
	for perm1, udiff1, mul_val1 in perms1:
		perms2 = elem_sym_perms_q(perm1,p2,k)
		#up_ref1 = get_reflections(u,iu,perm1)
		cycles1 = get_cycles(mulperm(iu,[*perm1]))
		#cycles1.sort(key=lambda c: c[-1])
		#cycles1_set = [set(c) for c in cycles1]
		cycles1_dict = {}
		for c in cycles1:
			if c[-1] not in cycles1_dict:
				cycles1_dict[c[-1]] = []
			cycles1_dict[c[-1]]+= [set(c)]
		#print(f"{cycles1_dict=}")
		ip1 = inverse(perm1)
		#tp1 = tuple(permtrim([*perm1]))
		#perm1i = inverse(perm1)
		for perm2, udiff2, mul_val2 in perms2:
			#up_ref2 = get_reflections(perm1,ip1,perm2)
			cycles2 = get_cycles(mulperm(ip1,[*perm2]))
			#cycles2_set = [set(c) for c in cycles2]
			good = True
			for i in range(len(cycles2)):
				c2 = cycles2[i]
				#c1_s = cycles1_set[i]
				if c2[-1] not in cycles1_dict:
					continue
				for c1_s in cycles1_dict[c2[-1]]:
					#print(f"{c1_s=}")
					for a in range(len(c2)-2,-1,-1):
						if c2[a] in c1_s:
							good = False
							break
					if not good:
						break														
				if not good:
					break
			
			#for a,b in up_ref1:
			#	for cycle in cycles2:
			#		if a in cycle and b in cycle:
			#			good = False
			#			break
			#	if not good:
			#		break
			#if not good:
			#	continue
			#for a,b in up_ref2:
			#	for cycle in cycles1:
			#		if a in cycle and b in cycle:
			#			good = False
			#			break
			#	if not good:
			#		break
			if good:				
				#ret[perm2] = ret.get(perm2,0) + val*mul_val1*mul_val2*elem_sym_func_q(k,1,u,perm1,(1,2),(1,2),udiff1,0,var2,var3)*elem_sym_func_q(k,2,perm1,perm2,(1,2),(1,2),udiff2,0,var2,var3)
				if (perm1,udiff1,mul_val1) not in ret_list:
					ret_list[(perm1,udiff1,mul_val1)] = []
				ret_list[(perm1,udiff1,mul_val1)] += [(perm2,udiff2,mul_val2)]
	return ret_list

def double_mul(coeff_dict,k):
	ret = {}
	s_k = [i+1 for i in range(k-1)] + [k+1] + [k]
	for u in coeff_dict:
		val = coeff_dict[u]
		perms1 = elem_sym_perms_q(u,k,k)
		for perm1, udiff1, mul_val1 in perms1:
			perms2 = elem_sym_perms_q(perm1,k,k)
			up_ref1, down_ref1 = get_reflections(u,perm1)
			cycles1 = get_cycles(mulperm(inverse(u),[*perm1]))
			down_cycles1 = set()
			for cycle in cycles1:
				for a,b in down_ref1:
					if a in cycle:
						down_cycles1.add(cycle)
			for perm2, udiff2, mul_val2 in perms2:
				up_ref2, down_ref2 = get_reflections(perm1,perm2)
				cycles2 = get_cycles(mulperm(inverse(perm1),[*perm2]))
				down_cycles2 = set()
				for cycle in cycles2:
					for a,b in down_ref2:
						if a in cycle:
							down_cycles2.add(cycle)
				good = True
				for a,b in up_ref1:
					for cycle in cycles2:
						if a in cycle and b in cycle:
						#if b in cycle:
							good = False
							break
					if not good:
						break
				
				for a,b in up_ref2:
					for cycle in cycles1:
						if a in cycle and b in cycle:
						#if b in cycle:
							good = False
							break
					if not good:
						break
				if good:				
					ret[perm2] = ret.get(perm2,0) + val*mul_val1*mul_val2*elem_sym_func_q(k,1,u,perm1,(1,2),(1,2),udiff1,0,var2,var3)*elem_sym_func_q(k,2,perm1,perm2,(1,2),(1,2),udiff2,0,var2,var3)
	return ret


def medium_theta(perm):
	cd = code(perm)
	found_one = True
	while found_one:
		found_one = False
		for i in range(len(cd)-1):
			if cd[i]<cd[i+1]: 
				found_one = True
				cd[i], cd[i+1] = cd[i+1]+1, cd[i]
				break
			if cd[i]==cd[i+1] and i>0 and cd[i-1]<=cd[i]+1:
				cd[i]+=1
				found_one = True
				break
	return cd

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
					inv_up = inv(up)
					newperms = double_elem_sym_q(up,mx_th,mx_th1,th[index])
					#for up1, up2, udiff1,udiff2,mul_val1,mul_val2 in newperms:
					for v in vpathdicts[index]:
						sumval = vpathsums[up].get(v,zero)
						if sumval == 0:
							continue
						for v2,vdiff2,s2 in vpathdicts[index][v]:
							for up1, udiff1, mul_val1 in newperms:
								esim1 = elem_sym_func_q(th[index],index+1,up,up1,v,v2,udiff1,vdiff2,var2,var3)*mul_val1*sumval*s2
								for up2, udiff2, mul_val2 in newperms[(up1,udiff1,mul_val1)]:
									if up2 not in newpathsums:
										newpathsums[up2]={}																
									for v3,vdiff3,s3 in vpathdicts[index+1][v2]:									
											newpathsums[up2][v3] = newpathsums[up2].get(v3,zero)+s3*mul_val2*esim1*elem_sym_func_q(th[index+1],index+2,up1,up2,v2,v3,udiff2,vdiff3,var2,var3)
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
	
	
def single_variable(coeff_dict,varnum):
	ret = {}
	for u in coeff_dict:
		if varnum -1 < len(u):
			ret[u] = ret.get(u,0) + (var2[u[varnum-1]]-var3[1])*coeff_dict[u]
		else:
			ret[u] = ret.get(u,0) + (var2[varnum] - var3[1])*coeff_dict[u]
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

#memory = Memory("/mnt/c/Math/pieri_cache", verbose=0)

glob_pieri = {}

def test_pieri(perm,k):
	if perm in glob_pieri:
		arr = glob_pieri[perm]
		if k in arr:
			return arr[k]
	coeff_dict = {perm: 1}
	if k == 0:
		return coeff_dict
	if k < 0:
		return {}
	ret_dict = {}
	
	#ret_dict = add_perm_dict(ret_dict,test_pieri(single_variable(coeff_dict,k),k-1))
	new_dict0 = single_variable(coeff_dict,k)
	for u in new_dict0:	
		#new_dict = test_pieri(coeff_dict,k-2)
		new_dict1 = test_pieri(u,k-1)
		#new_dict2 = test_pieri(u,k-2)
		#ret_dict[u] = ret_dict.get(u,0) + new_dict0[u]*test_pieri(u,k-1)
		for v in new_dict1:
			ret_dict[v] = ret_dict.get(v,0) + new_dict0[u]*new_dict1[v]
		#for v in new_dict2:
		#	ret_dict[v] = ret
	new_dict2 = test_pieri(perm,k-2)
	for u in new_dict2:
		ret_dict[u] = ret_dict.get(u,0)+q_var[k-1]*new_dict2[u]
	#glob_pieri[(perm,k)]=ret_dict
	arr = glob_pieri.get(perm,{})
	arr[k] = ret_dict
	glob_pieri[perm] = arr
	return ret_dict

def test_pieri_dict(coeff_dict,k):
	ret_dict = {}
	for u in coeff_dict:
		new_dict = test_pieri(u,k)
		for v in new_dict:
			ret_dict[v] = ret_dict.get(v,0) + coeff_dict[u]*new_dict[v]
	return ret_dict

def pieri_cycle(k):
	return tuple(permtrim(uncode([1 for i in range(k)])))
	
perms = list(permutations([i+1 for i in range(n)]))

new_perms = []

for perm in perms:
	new_perms += [tuple(permtrim([*perm]))]

var2_t = tuple(var2.tolist())
var3_t = tuple(var3.tolist())

def posify_dict(coeff_dict):
	ret = {}
	for perm in coeff_dict:
		val = coeff_dict[perm]
		val2 = 0
		q_dict = factor_out_q_keep_factored(val)
		for q_part in q_dict:
			try:
				val2 += q_part*int(q_dict[q_part])
			except Exception:
				val2 += q_part*compute_positive_rep(q_dict[q_part],var2_t,var3_t,False,False)				
		if val2 != 0:
			ret[perm] = val2
	return ret

def domul(perm):
	global permos
	for v in permos:
		#to_mul = uncode([2 for i in range(k)])
		coeff_dict = {perm: 1}
		#coeff_dict_test = schubmult(coeff_dict,v)
		coeff_dict_try = schubmult_db(coeff_dict,v)
		#fail = False
		#for u in coeff_dict_try:
		#	#print(f"{coeff_dict_try=}")
		#	if expand(coeff_dict_test.get(u,0) - coeff_dict_try.get(u,0)) != 0:
		#		#coeff_dict_test = posify_dict(coeff_dict_test)
		#		#coeff_dict_try = posify_dict(coeff_dict_try)
		#		print(f"Fail {perm} {v} {theta(inverse(v))=}")
		#		#print(f"Fail {perm} {k}",file=sys.stderr)				
		#		fail = True
		#		print(f"test={perm}: {u}: {coeff_dict_test.get(u,0)}")
		#		print(f"try={perm}: {u}: {coeff_dict_try.get(u,0)}")
		#		print(f"Fail {perm} {v} {theta(inverse(v))=}",file=sys.stderr)
		#		exit(1)
		#for u in coeff_dict_test:
		#	if expand(coeff_dict_test.get(u,0) - coeff_dict_try.get(u,0)) != 0:
		#		#coeff_dict_test = posify_dict(coeff_dict_test)
		#		#coeff_dict_try = posify_dict(coeff_dict_try)
		#		print(f"Fail {perm} {v} {theta(inverse(v))=}")
		#		#print(f"Fail {perm} {k}",file=sys.stderr)
		#		print(f"test={perm}: {u}: {coeff_dict_test.get(u,0)}")
		#		print(f"try={perm}: {u}: {coeff_dict_try.get(u,0)}")
		#		print(f"Fail {perm} {v} {theta(inverse(v))=}",file=sys.stderr)
		#		fail = True
		#		exit(1)
		#if not fail:
		#	print(f"Success {perm} {v}")
		#	print(f"Success {perm} {v}",file=sys.stderr)		

permos=new_perms
#for perm in new_perms:
#	#th = theta(inverse(perm))
#	#while len(th)>0 and th[-1] == 0:
#	#	th.pop()
#	#if len(th)==2 and th[0] == th[1]:
#	#	permos += [perm]
#	#good = True
#	#for i in range(len(th)):
#	#	if th.count(th[i])>2:
#	#		good = False
#	#		break
#	#for i in range(len(th)):
#	#	if i<len(th)-1 and th[i]==th[i+1] and i>0 and th[i-1]<=th[i]+1:
#	#		good = False
#	#		break
#	#if good and len(set(th))!=len(th):		
#	permos+=[perm]

#Parallel(n_jobs=11)(delayed(domul)(perm) for perm in new_perms)		
[domul(perm) for perm in new_perms]