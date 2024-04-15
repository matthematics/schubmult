from symengine import *
from functools import cache
from itertools import chain
from schubmult.perm_lib import *
import schubmult.schubmult_yz as yz
import sys

var = symarray('x', n)
var2 = symarray('y',n)
var3 = var2
var_r = symarray('r',n)

var_q = Symbol("q")

subs_dict = {}

for i in range(1,n):
	sm = var_r[0]
	for j in range(1,i):
		sm += var_r[j]
	subs_dict[var2[i]] = sm

def schubmult(perm_dict,v):
	vn1 = inverse(v)
	th = [len(v)-i for i in range(1,len(v)+1)]	
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
				newperms = elem_sym_perms_q(up,mx_th,th[index])
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
		try:
			for s in sys.argv[1:]:
				if s == "-np" or s == "--no-print":
					pr = False
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
				if s == "-":
					perms += [curperm]
					curperm = []
					continue
				curperm += [int(s)]
		except Exception:
			print("**** schubmult_q ****")
			print("Purpose: Compute the structure constants of quantum Schubert polynomials")
			print("Usage: schubmult_q <-np|--no-print> <-code> <-grass n> <-equiv> perm1 - perm2 < - perm 3 ... >")
			print("For the -grass option, must use Grassmannian permutations. -equiv only works together with -grass.")
			exit(1)
		
		perms += [curperm]
		
		
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
		
			coeff_dict = {tuple(permtrim([*perms[0]])): 1}
			
			for perm in perms[1:]:
				coeff_dict = schubmult(coeff_dict,tuple(permtrim([*perm])))
				
			if pr:
				if ascode:
					width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys() if expand(coeff_dict[perm])!=0])
				else:
					width = max([len(str(perm)) for perm in coeff_dict.keys() if expand(coeff_dict[perm])!=0])
				
				coeff_perms = list(coeff_dict.keys())
				coeff_perms.sort(key=lambda x: (inv(x),*x))
				
				for perm in coeff_perms:
					val = sympify(coeff_dict[perm]).expand()
					if val != 0:
						if ascode:					
							print(f"{str(trimcode(perm)):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
						else:
							print(f"{str(perm):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
	except BrokenPipeError:
		pass
			
if __name__ == "__main__":
	main()