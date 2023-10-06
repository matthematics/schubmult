import sys
from functools import cache
from itertools import chain
from schubmult.perm_lib import *

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

def will_formula_work(u,v):
	muv = uncode(theta(v))
	vn1muv = mulperm(inverse(v),muv)
	while True:
		found_one = False
		for i in range(len(vn1muv)-1):
			if vn1muv[i]>vn1muv[i+1]:				
				found_one = True
				if i<len(u)-1 and u[i]>u[i+1]:
					return False
				else:
					vn1muv[i],vn1muv[i+1] = vn1muv[i+1],vn1muv[i]
					break
		if not found_one:
			return True


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
		print("Usage: schubmult_py <-np> <perm1> - <perm2>")
		exit(1)
	
	perms += [tuple(permtrim(curperm))]
	
	#try:		
	#	u, v = perms
	#
	#	if will_formula_work(u,v) or will_formula_work(v,u):
	#		if will_formula_work(v,u):
	#			u, v = v, u
	#		print("Doing samuel formula",file=sys.stderr)
	#		inv_u = inv(u)
	#		inv_v = inv(v)	
	#		
	#		th = theta(v)
	#		muv = uncode(th)
	#		lm = p_trans(th)
	#		coeff_dict = schubmult({tuple(u): 1},muv)
	#		muvn1v = mulperm(inverse(muv),v)
	#		
	#		if pr:
	#			for perm, val in coeff_dict.items():
	#				w = mulperm(list(perm),muvn1v)
	#				if inv_u+inv_v == inv(w):
	#					print(f"{val}  {tuple(permtrim(w))}")
	#		sys.exit(0)
	#except Exception:
	#	pass
	
	perms.sort(reverse=True,key=lambda x: sum(theta(inverse(x)))-inv(x))
	
	coeff_dict = {perms[0]: 1}
	
	for perm in perms[1:]:
		coeff_dict = schubmult(coeff_dict,perm)
		
	if pr:
		for perm, val in coeff_dict.items():
			if val!= 0:
				print(f"{val}  {perm}")

if __name__ == "__main__":
	main()