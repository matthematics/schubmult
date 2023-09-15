from symengine import *
import sys
from functools import cache
from itertools import chain
from schubmult.perm_lib import *

n = 100

var = symarray('x', n)
var2 = symarray('y',n)
var3 = symarray('z',n)


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
	thL = len(th)
	vpathdicts = compute_vpathdicts(th,vmu)
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
				newperms = elem_sym_perms(up,min(mx_th,(inv_mu-(inv_up-inv_u))-inv_vmu),th[index])
				for up2, udiff in newperms:
					if up2 not in newpathsums:
						newpathsums[up2]={}
					for v in vpathdicts[index]:
						for v2,vdiff,s in vpathdicts[index][v]:
							newpathsums[up2][v2] = newpathsums[up2].get(v2,0)+s*vpathsums[up][v]*elem_sym_func(th[index],index+1,up,up2,v,v2,udiff,vdiff,var2,var3)
			vpathsums = newpathsums
		ret_dict = add_perm_dict({ep: vpathsums[ep][tuple(vmu)] for ep in vpathsums},ret_dict)
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
		print("Usage: python3 schubmult_yz.py <-np> <perm1> - <perm2>")
		exit(1)
	
	perms += [tuple(permtrim(curperm))]
	coeff_dict = {perms[0]: 1}
	
	for perm in perms[1:]:
		coeff_dict = schubmult(coeff_dict,perm)
		
	if pr:
		width = max([len(str(perm)) for perm in coeff_dict.keys()])
		
		coeff_perms = list(coeff_dict.keys())
		coeff_perms.sort(key=lambda x: (inv(x),*x))
		
		for perm in coeff_perms:
			val = coeff_dict[perm]
			if val != 0:
				print(f"{str(perm):>{width}}  {str(val).replace('**','^').replace('*',' ')}")
				
if __name__ == "__main__":
	main()
