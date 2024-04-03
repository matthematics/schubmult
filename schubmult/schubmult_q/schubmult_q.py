from symengine import *
from functools import cache
from itertools import chain
from schubmult.perm_lib import *
import sys

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


def main():
	try:			
		perms=[]
		curperm = []
		
		pr = True
		ascode = False
		try:
			for s in sys.argv[1:]:
				if s == "-np" or s == "--no-print":
					pr = False
					continue
				if s == "-code":
					ascode = True
					continue
				if s == "-":
					perms += [curperm]
					curperm = []
					continue
				curperm += [int(s)]
		except Exception:
			print("Usage: schubmult_q <-np|--no-print> <-code> perm1 - perm2 < - perm 3 ... >")
			exit(1)
		
		perms += [curperm]
		
		
		if ascode:
			for i in range(len(perms)):
				perms[i] = uncode(perms[i])
		
		coeff_dict = {tuple(permtrim([*perms[0]])): 1}
		
		for perm in perms[1:]:
			coeff_dict = schubmult(coeff_dict,tuple(permtrim([*perm])))
			
		if pr:
			if ascode:
				width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys()])
			else:
				width = max([len(str(perm)) for perm in coeff_dict.keys()])
			
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