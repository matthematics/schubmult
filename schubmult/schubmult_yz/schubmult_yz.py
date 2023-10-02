from symengine import *
import sys
from functools import cache
from itertools import chain
from schubmult.perm_lib import *
import numpy as np
import resource

n = 100

var = symarray('x', n)
var2 = symarray('y',n)
var3 = symarray('z',n)

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

def schubmult(perm_dict,v,var2=var2,var3=var3):
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
				newperms = elem_sym_perms(up,min(mx_th,(inv_mu-(inv_up-inv_u))-inv_vmu),th[index])
				for up2, udiff in newperms:
					if up2 not in newpathsums:
						newpathsums[up2]={}
					for v in vpathdicts[index]:
						sumval = vpathsums[up].get(v,0)
						if sumval == 0:
							continue
						for v2,vdiff,s in vpathdicts[index][v]:
							newpathsums[up2][v2] = newpathsums[up2].get(v2,zero)+s*sumval*elem_sym_func(th[index],index+1,up,up2,v,v2,udiff,vdiff,var2,var3)
			vpathsums = newpathsums
		toget = tuple(vmu)
		ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget,0) for ep in vpathsums},ret_dict)
	return ret_dict



dimen = 0
monom_to_vec = {}

def main():
	sys.setrecursionlimit(1000000)
	resource.setrlimit(resource.RLIMIT_STACK, (-1, -1))

	perms=[]
	curperm = []
	
	pr = True
	display_positive = False
	
	try:
		for s in sys.argv[1:]:
			if s == "-np" or s == "--no-print":
				pr = False
				continue
			if s == "--display-positive":
				display_positive = True
				continue
			if s == "--version":
				print(f"Python version {sys.version}")
				exit(0)
			if s == "--usage":
				print("Usage: schubmult_yz <-np|--no-print> <--display-positive> <perm1> - <perm2>")
				exit(0)
			if s == "-":
				perms += [tuple(permtrim(curperm))]
				curperm = []
				continue
			curperm += [int(s)]
	except Exception:
		print("Usage: schubmult_yz <-np|--no-print> <--display-positive> <perm1> - <perm2>")
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
				if display_positive:
					notint = False
					try:
						int(val)
					except Exception:
						notint = True
					if notint:
						try:
							import pulp as pu
						except ImportError:
							print("Package \"pulp\" required for --display-positive",file=sys.stderr)
							exit(1)
							
						import itertools as it
						try:
							import sympy
						except ImportError:
							print("Package \"sympy\" required for --display-positive",file=sys.stderr)
							exit(1)					
																			
						sval = str(val)
						val = expand(val)
						import re
						
						regex = re.compile("y_([0-9]+)")
											
						n1 = max([int(index) for index in regex.findall(sval)])+1
						n2 = max([len(L) for L in perms[1:]])
						var22 = sympy.symarray('y',n1)
						var33 = sympy.symarray('z',n2)
						
						def poly_to_vec(poly,vec0=None):	
							global dimen, monom_to_vec
							dc = poly.as_coefficients_dict()
							mn = list(dc.keys())
							
							if vec0 is None:
								init_basevec(mn)
							vec = np.zeros(dimen)
							try:
								for i in range(len(mn)-1,-1,-1):
									cf = dc[mn[i]]
									index = monom_to_vec[mn[i]]
									if vec0 is not None and abs(vec0[index])<abs(cf):
										return None
									vec[index] = cf
							except Exception:
								return None
							return vec
							
						def init_basevec(monoms):
							global dimen, monom_to_vec
							monom_to_vec = {}
							index = 0
							for mn in monoms:
								monom_to_vec[mn] = index
								index+=1
							dimen = index
						val_poly = sympy.poly(sval,*var22[1:],*var33[1:])
					
						test_degree = inv(perms[0])+inv(perms[1])-inv(perm)
						did = False						
						if sympy.total_degree(val_poly)<test_degree or val_poly == 0:
							continue
						
						list0 = []
						list1 = []
						mn = val_poly.monoms()
						
						L0 = [0 for i in range(1,n1)]
						L1 = [0 for i in range(1,n2)]
						mn1L = []
						mnn1L = []
						for mm0 in mn:
							if list(mm0[:n1-1]) == L0:
								mnn1L += [mm0]
							if list(mm0[n1-1:]) == L1:
								mn1L += [mm0]
						
						base = []
						base_vectors = []
						base_monoms = []
						
						vec = poly_to_vec(val)
						
						set1 = set()
						set2 = set()
						full_varlist = []
						for mnn1 in mnn1L:
							for i in range(len(mnn1)):
								if mnn1[i]!=0:
									set2.add(i+1-n1+1)
						for index1 in range(1,n1):
							full_varlist += [[]]
							for index2 in set2:
								full_varlist[index1-1] += [var2[index1]-var3[index2]]
						comblists = []
						for mn1 in mn1L:				
							comblistmn1 = [1]
							for i in range(n1-1):					
								if mn1[i]!=0:
									arr = np.array(comblistmn1)
									comblistmn12 = []
									newset = it.combinations(full_varlist[i],mn1[i])
									bobo = list(newset)
									for bob in bobo:
										boblot = np.prod(bob)
										comblistmn12 += (arr*boblot).tolist()
									comblistmn1 = comblistmn12
							for b1 in comblistmn1:
								vec0 = poly_to_vec(expand(b1),vec)
								if vec0 is not None:
									base_vectors += [vec0]
									base_monoms += [b1]
						vrs = [pu.LpVariable(name=f"a{i}",lowBound=0,cat='Integer') for i in range(len(base_vectors))]
						lp_prob = pu.LpProblem('Problem',pu.LpMinimize)
						lp_prob += 0
						for i in range(len(vec)):
							eq = 0
							for j in range(len(base_vectors)):
								if base_vectors[j][i]!=0:
									eq += int(base_vectors[j][i])*vrs[j]
							lp_prob += eq == int(vec[i])
						try:
							status = lp_prob.solve(pu.XPRESS_PY(msg=0))
						except pu.PulpSolverError:
							status = lp_prob.solve(pu.PULP_CBC_CMD(msg=0))
						if status != 1:
							print(f"Counterexample {perms[0]} {perms[1]} {perm} {val}")
							exit(1)
						val = 0
						for i in range(len(base_vectors)):
							x = vrs[i].value()
							b1 = base_monoms[i]
							if x!=0:
								val += int(x)*b1
				if val!=0:
					print(f"{str(perm):>{width}}  {str(val).replace('**','^').replace('*',' ')}")	
if __name__ == "__main__":
	main()
