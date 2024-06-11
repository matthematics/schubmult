from schubmult.perm_lib import *
from schubmult.schubmult_q_yz import nil_hecke, var2, var3
from schubmult.schubmult_yz import schubpoly
from itertools import permutations
import sys

n = int(sys.argv[1])

perms = list(permutations([i+1 for i in range(n)]))

perms.sort()

dom = 0
num = 0

def pieri_cycle(p,k):
	return tuple(permtrim(uncode([0 for i in range(k-p)] + [1 for i in range(k-p,k)])))

def elem_sym_nil(p,k,v,n,var2=var2,var3=var3):
	#coeff_dict = nil_hecke({(1,2): 1},pieri_cycle(p,k),var2,var2)
	eperms = elem_sym_perms_q_op((1,2),p,k,n)
	
	#print(f"{p=} {k=}")
	ret = 0	
	#print(f"{eperms=}")
	
	doneperms = set()
	
	for u, udiff, mul_val in eperms:		
		print(f"{u=} {udiff=}")
		pc = pieri_cycle(p-udiff,k-udiff)
		#print(f"{code(pc)=}")
		tomul = inverse(pc)
		u2 = tuple(permtrim(mulperm(tomul,u)))		
		print(f"{u2=}")
		if u2 in doneperms:
			print("Already done")
			continue		
		inv_u2 = inv(u2)
		if inv_u2 == inv(u) - p + udiff:
			# do div diff
			doneperms.add(u2)
			v2 = mulperm(v,inverse(u2))
			if inv(v2) == inv(v) - inv_u2:
				spoly = sympify(schubpoly(v2,var2,var3))
				spoly2 = spoly.subs({var2[i+1]: var2[getpermval(pc,i)] for i in range(20)})
				print(f"Adding {spoly2=} from {spoly=} {pc=} {mul_val=}")
				ret += mul_val*spoly2
	return ret

for perm in perms:
	v = permtrim([*perm])
	if v == [1,2]:
		continue
	nh = nil_hecke({(1,2): 1},tuple(v),n,var2,var3)
	for k in range(1,n):
		for p in range(k,k+1):
			pc = pieri_cycle(p,k)
			test_val = nh.get(pc,0)
			new_val = elem_sym_nil(p,k,v,n)
			if expand(test_val - new_val) != 0:
				print(f"Fail {p=} {k=} {v=} {pc=} {test_val=} {new_val=}")
				#exit(1)
			else:
				print(f"Success {v=} {p=} {k=}")
	
		
