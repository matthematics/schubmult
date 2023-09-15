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
	for u,val in perm_dict.items():
		inv_u = inv(u)
		vpathsums = {u: {(1,2): val}}
		for index in range(len(th)):			
			newpathsums = {}
			for up in vpathsums:
				inv_up = inv(up)
				newperms = elem_sym_perms(up,min(th[index],(inv_mu-(inv_up-inv_u))-inv_vmu),th[index])
				for up2, udiff in newperms:					
					for v in vpathsums[up]:
						for v2,vdiff,s in vpathdicts[index][v]:
							if vdiff + udiff == th[index]:								
								if up2 not in newpathsums:
									newpathsums[up2]={}
								newpathsums[up2][v2] = newpathsums[up2].get(v2,0)+s*vpathsums[up][v]
			vpathsums = newpathsums
		ret_dict = add_perm_dict({ep: vpathsums[ep].get(tuple(vmu),0) for ep in vpathsums},ret_dict)
	return ret_dict

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
		print("Usage: python3 schubmult_py.py <-np> <perm1> - <perm2>")
		exit(1)
	
	perms += [tuple(permtrim(curperm))]
	
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