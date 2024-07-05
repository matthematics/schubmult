from schubmult.perm_lib import *
from schubmult.schubmult_q_yz import schubmult, schubmult_db, mult_poly, factor_out_q_keep_factored
from symengine import *
import sys


q_var = symarray("q",100)

var2 = symarray("y",100)
var3 = symarray("z",100)

var3 = var2
var_r = symarray('r',100)

subs_dict = {}

for i in range(1,100):
	sm = var2[1]
	for j in range(1,i):
		sm += var_r[j]
	subs_dict[var2[i]] = sm
	
	

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
		check = True
		msg = False
		mult = False
		slow = False
		mulstring = ""
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
				if s == "--version":
					print(f"Python version {sys.version}")
					exit(0)
				if s == "-code":
					ascode = True
					continue
				if s == "--usage":
					print("**** schubmult_q_double ****")
					print("Purpose: Compute equivariant Gromov-Witten invariants, structure constants of quantum double Schubert polynomials")
					print("Usage: schubmult_q_double <-np|--no-print> <-parabolic descent1,descent2,...> <-code> <--display-positive> <--optimizer-message> perm1 - perm2 < - perm3 .. >")
					#print("Alternative usage: schubmult_yz <-code> <--display-positive> -coprod perm - indexlist")
					exit(0)
				if s == "-":
					perms += [curperm]
					curperm = []
					continue
				curperm += [int(s)]
		except Exception:
			print("**** schubmult_q_double ****")
			print("Purpose: Compute equivariant Gromov-Witten invariants, structure constants of quantum double Schubert polynomials")
			print("Usage: schubmult_q_double <-np|--no-print> <-code> <-parabolic descent1,descent2,...> <--display-positive> <--optimizer-message> perm1 - perm2 < - perm3 .. >")
			exit(1)
				
		perms += [curperm]
		
		
		if ascode:
			for i in range(len(perms)):
				perms[i] = tuple(permtrim(uncode(perms[i])))
		else:
			for i in range(len(perms)):
				perms[i] = tuple(permtrim([*perms[i]]))
		
		size = 0
		L = len(perms)
		
		coeff_dict = {perms[0]: 1}
		for perm in perms[1:]:
			if slow:
				coeff_dict = schubmult(coeff_dict,perm,var2,var2)
			else:
				coeff_dict = schubmult_db(coeff_dict,perm,var2,var2)
		if mult:
			mul_exp = sympify(mulstring)
			coeff_dict = mult_poly(coeff_dict,mul_exp)
		
		if pr:
		
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
		
			if ascode:
				width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys() if expand(sympify(coeff_dict[perm]).subs(subs_dict))!=0])
			else:
				width = max([len(str(perm)) for perm in coeff_dict.keys() if expand(sympify(coeff_dict[perm]).subs(subs_dict))!=0])
			
			coeff_perms = list(coeff_dict.keys())
			coeff_perms.sort(key=lambda x: (inv(x),*x))
			
			for perm in coeff_perms:
				val = expand(sympify(coeff_dict[perm]).subs(subs_dict)).simplify()
				if val != 0:
					notint = False
					try:
						int(val)
					except Exception:
						notint = True
					if val!=0:
						if ascode:
							print(f"{str(trimcode(perm))}  {str(val).replace('**','^').replace('*',' ')}")	
						else:
							print(f"{str(perm)}  {str(val).replace('**','^').replace('*',' ')}")	
	except BrokenPipeError:
		pass
		
if __name__ == "__main__":
	main()