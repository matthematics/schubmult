from symengine import *
import sys
from functools import cache
from itertools import chain
from schubmult.perm_lib import *
import numpy as np
import pulp as pu
import itertools as it
import sympy
import re
		
n = 100

var = symarray('x',n)
var2 = symarray('y',n)
var3 = symarray('z',n)

def forwardcoeff(u,v,perm,var2=var2,var3=var3):
	th = theta(v)
	muv = uncode(th)
	vmun1 = mulperm(inverse(list(v)),muv)
	
	w = mulperm(list(perm),vmun1)	
	if inv(w) == inv(vmun1) + inv(perm):
		coeff_dict = {tuple(permtrim([*u])): 1}
		coeff_dict = schubmult(coeff_dict,muv,var2,var3)
		return coeff_dict.get(tuple(permtrim(w)),0)
	return 0

def dualcoeff(u,v,perm,var2=var2,var3=var3):
	th = theta(u)
	muu = uncode(th)
	umun1 = mulperm(inverse(list(u)),muu)
	if u == (1,2):
		vp = mulperm(list(v),inverse(perm))
		if inv(vp) == inv(v) - inv(perm):
			val = schubpoly(vp,var2,var3)
		else:
			val = 0
	elif dominates(u,perm):
		ret = 0
		dpret = dualpieri(list(u),list(v),list(perm))
		for vlist,vp in dpret:
			toadd = 1
			for i in range(len(vlist)):										
				for j in range(len(vlist[i])):
					toadd *= var2[i+1]-var3[vlist[i][j]]
			toadd *= schubpoly(vp,var2,var3,len(vlist)+1)
			ret += toadd		
		val = ret
	else:
		ret = 0
		w = mulperm(list(perm),umun1)
		if inv(w) == inv(umun1) + inv(perm):
			ret = 0
			dpret = dualpieri(muu,list(v),w)
			for vlist,vp in dpret:
				toadd = 1
				for i in range(len(vlist)):										
					for j in range(len(vlist[i])):
						toadd *= var2[i+1]-var3[vlist[i][j]]
				toadd *= schubpoly(vp,var2,var3,len(vlist)+1)
				ret += toadd		
		val = ret
	return val


def dualpieri(mu,v,w):
	lm = code(inverse(mu))
	cn1w = code(inverse(w))
	while len(lm)>0 and lm[-1] == 0:
		lm.pop()
	while len(cn1w)>0 and cn1w[-1]==0:
		cn1w.pop()
	if len(cn1w)<len(lm):
		return []
	for i in range(len(lm)):
		if lm[i]>cn1w[i]:
			return []
	c = [1,2]
	for i in range(len(lm),len(cn1w)):
		c = mulperm(cycle(i-len(lm)+1,cn1w[i]),c)
	c = permtrim(c)
	res = [[[],v]]
	for i in range(len(lm)):
		res2 = []	
		for vlist,vplist in res:
			vp = vplist
			# cycle c_{1,lm[i]}^{-1}c_{1,cn1w[i]}
			vpl = divdiffable(vp,cycle(lm[i]+1,cn1w[i]-lm[i]))
			if vpl == []:
				continue
			# pull out var lm[i]+1
			vl = pull_out_var(lm[i]+1,vpl)
			for pw,vpl2 in vl:
				res2 += [[vlist+[pw],vpl2]]
		res = res2
	if len(lm) == len(cn1w):
		return res
	else:
		res2 = []
		for vlist,vplist in res:
			vp = vplist
			vpl = divdiffable(vp,c)
			if vpl == []:
				continue
			res2 += [[vlist,vpl]]
		return res2
		


dimen = 0
monom_to_vec = {}


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

def poly_to_vec(poly,vec0=None):	
	global dimen, monom_to_vec, n1, n2
	poly = expand(poly.subs({var3[1]: 0}))
	dc = poly.as_coefficients_dict()
	mn = list(dc.keys())
		
	
	if vec0 is None:
		init_basevec(dc)
	vec = np.zeros(dimen)
	
	for i in range(len(mn)-1,-1,-1):
		cf = dc[mn[i]]
		if cf == 0:
			continue
		cf = abs(cf)
		try:
			index = monom_to_vec[mn[i]]
		except Exception:
			return None
		if vec0 is not None and abs(vec0[index])<abs(cf):
			return None
		vec[index] = cf	
	return vec

def shiftsub(pol):
	subs_dict = dict([(var2[i],var2[i+1]) for i in range(99)])
	return sympify(pol).subs(subs_dict)

def shiftsubz(pol):
	subs_dict = dict([(var3[i],var3[i+1]) for i in range(99)])
	return sympify(pol).subs(subs_dict)	

def init_basevec(dc):
	global dimen, monom_to_vec
	monom_to_vec = {}
	index = 0
	for mn in dc:
		if dc[mn] == 0:
			continue
		monom_to_vec[mn] = index
		index+=1
	dimen = index
	

def split_flat_term(arg):
	arg = expand(arg)
	ys = []
	zs = []
	for arg2 in arg.args:
		if str(arg2).find("y") != -1:
			if isinstance(arg2,Mul):
				for i in range(int(arg2.args[0])):
					ys += [arg2.args[1]]
			else:
				ys += [arg2]
		else:
			if isinstance(arg2,Mul):
				for i in range(abs(int(arg2.args[0]))):
					zs += [-arg2.args[1]]
			else:
				zs += [arg2]
	return ys, zs

def is_flat_term(term):
	if isinstance(term,Integer) or isinstance(term,int):
		return True
	dc = expand(term).as_coefficients_dict()
	for t in dc:
		if str(t).count("y")+str(t).count("z")>1 or str(t).find("**")!=-1:
			return False
	return True

def flatten_factors(term,var2=var3,var3=var3):
	found_one = False
	if is_flat_term(term):
		return term, False
	elif isinstance(term,Pow):
		if is_flat_term(term.args[0]) and len(term.args[0].args)>2:
			ys, zs = split_flat_term(term.args[0])
			terms = [1]
			for i in range(len(ys)):
				terms2 = []
				for j in range(len(term.args[1])):
					for t in terms:
						terms2 += [t*(ys[i]+zs[i])]						
				terms = terms2			
			return Add(*terms)
		elif is_flat_term(term.args[0]):
			return term, False
		else:
			return flatten_factors(term.args[0])**term.args[1], True
	elif isinstance(term,Mul):
		terms = [1]
		for arg in term.args:	
			terms2 = []
			if isinstance(arg,Add) and not is_flat_term(expand(arg)):
				found_one = True				
				for term3 in terms:
					for arg2 in arg.args:
						flat, found = flatten_factors(arg2)
						terms2 += [term3*flat]				
			elif isinstance(arg,Add) and is_flat_term(arg) and len(arg.args)>2:
				found_one = True				
				ys, zs = split_flat_term(arg)
				for term3 in terms:
					for i in range(len(ys)):
						terms2 += [term3*(ys[i] + zs[i])]
			else:
				flat, found = flatten_factors(arg)
				if found:
					found_one = True
				for term3 in terms:					
					terms2+=[term3*flat]
			terms = terms2
		if len(terms) == 1:
			term = terms[0]
		else:
			term = Add(*terms)
		return term, found_one
	elif isinstance(term,Add):
		res = 0
		for arg in term.args:
			flat, found = flatten_factors(arg)
			if found:
				found_one = True
			res += flat
		return res, found_one

var2s = {var2[i]: i for i in range(len(var2))}
var3s = {var3[i]: i for i in range(len(var3))}

def split_mul(arg0):
	monoms = SortedList()
	for arg in arg0.args:
		if is_flat_term(arg):
			if isinstance(arg,Integer) or isinstance(arg,int):
				continue
			#print(type(arg))
			yval = arg.args[0]
			zval = arg.args[1]
			if str(yval).find("z")!=-1:
				yval, zval = zval, yval
			if str(zval).find("-") != -1:
				zval = -zval
			if str(yval).find("-") != -1:
				yval = -yval
			monoms += [(var2s[yval],var3s[zval])]
		elif isinstance(arg,Pow):
			arg2 = arg.args[0]
			yval = arg2.args[0]
			zval = arg2.args[1]
			if str(yval).find("z")!=-1:
				yval, zval = zval, yval
			if str(zval).find("-") != -1:
				zval = -zval
			if str(yval).find("-") != -1:
				yval = -yval
			tup = (var2s[yval],var3s[zval])
			for i in range(int(arg.args[1])):
				monoms += [tup]
	return monoms

def split_monoms(pos_part):
	arrs = SortedList()	
	if isinstance(pos_part,Add):		
		for arg0 in pos_part.args:
			monoms = split_mul(arg0)			
			arrs += [monoms]
	elif isinstance(pos_part,Mul) or isinstance(pos_part,Pow):
		arrs += [split_mul(pos_part)]
	else:
		return [pos_part]
	return arrs
	
def is_negative(term):
	sign = 1
	if isinstance(term,Integer) or isinstance(term,int):
		return term < 0
	elif isinstance(term,Mul):
		for arg in term.args:
			if isinstance(arg,Integer):
				sign *= arg
			elif isinstance(arg,Add):
				if str(arg).find("-y") != -1:
					sign *= -1
			elif isinstance(arg,Pow):
				mulsign = 1
				if str(arg.args[0]).find("-y") != -1:
					mulsign = -1
				sign *= mulsign**term.args[1]
	elif isinstance(term,Pow):
		mulsign = 1
		if str(term.args[0]).find("-y") != -1:
			mulsign = -1
		sign *= mulsign**term.args[1]
	return sign < 0


def compute_positive_rep(val,var2=var2,var3=var3,msg=False):
	notint = False	
	try:
		int(val)
		val2 = val
	except Exception:
		notint = True
	if notint:
		sval = str(val)
		regex = re.compile("y_([0-9]+)")
		regex2 = re.compile("z_([0-9]+)")					
		n1 = max([int(index) for index in regex.findall(sval)])+1
		n2 = max([int(index) for index in regex2.findall(sval)])+1

		var22 = sympy.symarray('y',n1)
		var33 = sympy.symarray('z',n2)			
		
		vec = poly_to_vec(val)

		base = []
		base_vectors = []
		base_monoms = []						
		
		val_poly = sympy.poly(sval,*var22[1:],*var33[1:])
		mn = val_poly.monoms()
		
		L1 = [0 for i in range(1,n1)]
		mn1L = []
		for mm0 in mn:
			if list(mm0[:n1-1]) == L1:
				mn1L += [mm0]			
	
		for mn1 in mn1L:				
			comblistmn1 = [1]
			for i in range(n2-1):
				if mn1[n1-1+i]!=0:
					arr = np.array(comblistmn1)
					comblistmn12 = []
					mn1_2 = tuple(list(mn1[:n1-1+i])+[0]+list(mn1[n1+i:]))
					newset	= []
					for mm0 in mn:
						if mm0[n1-1:] == mn1_2[n1-1:]:
							toadd = 1
							bad = False
							st = set([mm0[k] for k in range(n1-1)])
							if len(st.intersection(set([0,1])))<len(st):
								bad = True
							else:
								for k in range(n1-1):
									if mm0[k] == 1:
										toadd *= var2[k+1] - var3[i+1]
								if not bad:
									newset += [toadd]		
					for term in newset:
						comblistmn12 += (arr*term).tolist()
					comblistmn1 = comblistmn12
			for b1 in comblistmn1:
				vec0 = poly_to_vec(b1,vec)
				if vec0 is not None:
					base_vectors += [vec0]
					base_monoms += [b1]	
		vrs = [pu.LpVariable(name=f"a{i}",lowBound=0,cat='Integer') for i in range(len(base_vectors))]
		lp_prob = pu.LpProblem('Problem',pu.LpMinimize)		
		lp_prob += int(0)
		last = 0
		indices = [i for i in range(len(base_vectors[0]))]
		
		for i in indices:
			eq = 0
			for j in range(len(base_vectors)):
				if base_vectors[j][i]!=0:
					eq += int((base_vectors[j][i]))*vrs[j]
			lp_prob += eq == int(vec[i])			

		try:
			status = lp_prob.solve(pu.XPRESS_PY(msg=msg))
		except pu.PulpSolverError:
			try:
				solver = pu.PULP_CBC_CMD(msg=msg)
				status = lp_prob.solve(solver)
			except KeyboardInterrupt:
				import psutil
				current_process = psutil.Process()
				children = current_process.children(recursive=True)
				for child in children:
					child_process = psutil.Process(child.pid)
					child_process.terminate()
					child_process.kill()
				raise KeyboardInterrupt()

		val2 = 0
		for k in range(len(base_vectors)):
			x = vrs[k].value()
			b1 = base_monoms[k]
			if x!=0 and x is not None:
				val2 += int(x)*b1
	return val2

def posify(val,u2,v2,w2,var2=var2,var3=var3,msg=False):
	if inv(u2)+inv(v2) - inv(w2)<=1:
		return expand(val)
	if expand(val) == 0:
		return 0
	
	u, v, w = try_reduce_v(u2, v2, w2)
	if not will_formula_work(u,v) and not will_formula_work(v,u) and not one_dominates(u,w) and not is_reducible(v) and inv(w) - inv(u)!=1:
		u, v, w = try_reduce_u(u2, v2, w2)
		if not will_formula_work(u,v) and not will_formula_work(v,u) and not one_dominates(u,w) and not is_reducible(v) and inv(w) - inv(u)!=1:
			u, v, w = [*u2], [*v2], [*w2]
			if not will_formula_work(u,v) and not will_formula_work(v,u) and not one_dominates(u,w) and not is_reducible(v) and inv(w) - inv(u)!=1:
				w0 = list(w)
				u, v, w = reduce_descents(u,v,w)
				if not will_formula_work(u,v) and not will_formula_work(v,u) and not one_dominates(u,w)  and not is_reducible(v) and inv(w) - inv(u)!=1:
					u, v, w = reduce_coeff(u,v,w)
					if not will_formula_work(u,v) and not will_formula_work(v,u) and not one_dominates(u,w)  and not is_reducible(v) and inv(w) - inv(u)!=1:
						while not will_formula_work(u,v) and not will_formula_work(v,u) and tuple(permtrim(w0))!=tuple(permtrim(list(w))) and not one_dominates(u,w) and inv(w) - inv(u)!=1 and not is_reducible(v):
							w0 = w
							u, v, w = reduce_descents(u,v,w)						
							if not will_formula_work(u,v) and not will_formula_work(v,u) and not one_dominates(u,w) and not is_reducible(v) and inv(w) - inv(u)!=1:	
								u, v, w = reduce_coeff(u,v,w)
	u = tuple(u)
	v = tuple(v)	
	w = tuple(w)
			
	
	u3, v3, w3 = try_reduce_v(u, v, w)
	if is_reducible(v3):
		u, v, w = u3, v3, w3
	else:
		u3, v3, w3 = try_reduce_u(u, v, w)
		if one_dominates(u3,w3):
			u, v, w = u3, v3, w3
	if (will_formula_work(v,u) or dominates(u,w)):
		val = dualcoeff(u,v,w,var2,var3)
	elif inv(w) - inv(u) == 1:	
		a, b = -1, -1
		for i in range(len(w)):
			if a == -1 and u[i] != w[i]:
				a = i
			elif i>=len(u) and w[i] != i+1:
				b = i
			elif b == -1 and u[i] != w[i]:
				b = i
		arr = [[[],v]]
		d = -1
		for i in range(len(v)-1):
			if v[i]>v[i+1]:
				d = i + 1
		for i in range(d):
			arr2 = []
			if i in [a,b]:
				continue
			i2 = 1
			if i>b:
				i2 += 2
			elif i>a:
				i2 += 1
			for vr, v2 in arr:
				dpret = pull_out_var(i2,[*v2])
				for v3r, v3 in dpret:
					arr2 += [[vr + [v3r],v3]]
			arr = arr2
		val = 0
		for L in arr:
			v3 = [*L[-1]]
			if v3[0]<v3[1]:
				continue
			else:
				v3[0], v3[1] = v3[1], v3[0]
			toadd = 1
			for i in range(d):
				if i in [a,b]:
					continue
				i2 = i
				if i>b:
					i2 = i-2
				elif i>a:
					i2 = i-1
				oaf = L[0][i2]
				if i>=len(w):
					yv = i + 1
				else:
					yv = w[i]
				for j in range(len(oaf)):				
					toadd*= var2[yv] - var3[oaf[j]]
			toadd *= schubpoly(v3,[0,var2[w[a]],var2[w[b]]],var3)
			val += toadd
	elif expand(val) != 0:	
		c01 = code(u)
		c02 = code(w)
		c03 = code(v)
		
		c1 = code(inverse(u))
		c2 = code(inverse(w))
		
		if one_dominates(u,w):
			while c1[0] != c2[0]:				
				w = list(w)
				v = list(v)
				w[c2[0]-1], w[c2[0]] = w[c2[0]], w[c2[0]-1]
				v[c2[0]-1], v[c2[0]] = v[c2[0]], v[c2[0]-1]
				w = tuple(w)
				v = tuple(v)
				c2 = code(inverse(w))
				c03 = code(v)
				c01 = code(u)
				c02 = code(w)				
				
		
		if is_reducible(v):
			newc = []
			elemc = []
			for i in range(len(c03)):
				if c03[i]>0:
					newc += [c03[i]-1]
					elemc += [1]
				else:
					break
			v3 = uncode(newc)
			coeff_dict = schubmult({u: 1},tuple(permtrim(uncode(elemc))),var2,var3)
			val = 0
			for new_w in coeff_dict:
				tomul = coeff_dict[new_w]
				newval = schubmult({new_w: 1},tuple(permtrim(uncode(newc))),var2,var3).get(tuple(permtrim(list(w))),0)
				newval = posify(newval,new_w,permtrim(uncode(newc)),w,var2,var3,msg)
				val += tomul*shiftsubz(newval)
		elif c01[0] == c02[0] and c01[0] != 0:
			varl = c01[0]
			u3 = uncode([0] + c01[1:])
			w3 = uncode([0] + c02[1:])
			val = 0
			val = schubmult({tuple(permtrim(u3)): 1},tuple(permtrim(list(v))),var2,var3).get(tuple(permtrim(w3)),0)
			val = posify(val,u3,v,w3,var2,var3,msg)
			for i in range(varl):
				val = permy(val,i+1)			
		elif c1[0] == c2[0]:
			vp = pull_out_var(c1[0]+1,list(v))
			u3 = phi1(u)
			w3 = phi1(w)
			c3 = code(inverse(u3))
			c4 = code(inverse(w3))
			val = 0
			for arr, v3 in vp:
				tomul = 1
				for i in range(len(arr)):
					tomul*=var2[1] - var3[arr[i]]
				
				val2 = schubmult({tuple(permtrim(u3)): 1},tuple(permtrim(v3)),var2,var3).get(tuple(permtrim(w3)),0)
				val2 = posify(val2,u3,v3,w3,var2,var3,msg)
				val += tomul*shiftsub(val2)
		else:
			val2 = compute_positive_rep(val,var2,var3,msg)
			if val2 is not None:
				val = val2
	else:
		val = 0	
	return val

def split_perms(perms):
	perms2 = [perms[0]]
	for perm in perms[1:]:
		cd = code(perm)
		index = -1
		not_zero = False
		did = False
		for i in range(len(cd)):
			if cd[i] != 0:
				not_zero = True
			elif not_zero and cd[i] == 0:
				not_zero = False
				index = i
				num_zeros_to_miss = 0
				for j in range(index):
					if cd[j]!=0:
						num_zeros_to_miss = max(num_zeros_to_miss,cd[j]-(index-1-j))
				num_zeros = 0
				for j in range(index,len(cd)):
					if cd[j]!=0:
						break
					else:
						num_zeros+=1
				if num_zeros >= num_zeros_to_miss:
					cd1 = cd[:index]
					cd2 = [0 for i in range(index)]+cd[index:]
					perms2 += [tuple(permtrim(uncode(cd1))),tuple(permtrim(uncode(cd2)))]
					did = True
					break
		if not did:
			perms2 += [perm]
	return perms2

def schubpoly(v,var2=var2,var3=var3,start_var = 1):
	n = 0
	for j in range(len(v)-2,-1,-1):
		if v[j]>v[j+1]:
			n = j+1
			break
	if n == 0:
		return 1
	lst = pull_out_var(n,v)	
	ret = 0	
	for pw, vp in lst:
		tomul = 1
		for p in pw:
			tomul *= var2[start_var+n-1]-var3[p]
		ret += tomul * schubpoly(vp,var2,var3,start_var)	
	return ret

def permy(val,i):
	subsdict = {var2[i]: var2[i+1], var2[i+1]: var2[i]}
	return sympify(val).subs(subsdict)

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
		try:
			for s in sys.argv[1:]:
				if s == "-np" or s == "--no-print":
					pr = False
					continue
				if s == "-coprod":
					coprod = True
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
					print("Usage: schubmult_yz <-np|--no-print> <-code> <--display-positive> <--optimizer-message> perm1 - perm2 < - perm3 .. >")
					print("Alternative usage: schubmult_yz <-code> <--display-positive> -coprod perm - indexlist")
					exit(0)
				if s == "-":
					perms += [curperm]
					curperm = []
					continue
				curperm += [int(s)]
		except Exception:
			print("Usage: schubmult_yz <-np|--no-print> <-code> <--display-positive> <--optimizer-message> perm1 - perm2 < - perm3 .. >")
			print("Alternative usage: schubmult_yz <-code> <--display-positive> <--optimizer-message> -coprod perm - indexlist")
			exit(1)
				
		perms += [curperm]
		posified = False
		if coprod:
			subs_dict = {}
			if ascode:
				mperm = tuple(permtrim(uncode(perms[0])))
			else:
				mperm = tuple(permtrim(perms[0]))
	
			cd = code(mperm)
			perms[0] = mperm
			pos = perms[1]
	
			while cd[-1] == 0:
				cd.pop()
			k = len(pos)
			n = len(perms[0])
			kcd = [pos[i]-i-1 for i in range(len(pos))] + [n+1-k for i in range(k,n)]
			N = len(kcd)
			
			var2neg = [-var2[i] for i in range(len(var2))]
			var3neg = [-var3[i] for i in range(len(var3))]
			
			for i in range(1,100):
				if i<=N:
					subs_dict[var[i]] = var2[i]
				else:
					subs_dict[var[i]] = var3[i-N]
			
			kperm = inverse(uncode(kcd))
			coeff_dict = {tuple(kperm): 1}
			coeff_dict = schubmult(coeff_dict,perms[0],var,var2)
			
			inv_perm0 = inv(perms[0])
			inv_kperm = inv(kperm)
			inverse_kperm = inverse(kperm)
			if pr:
				coeff_perms = list(coeff_dict.keys())
				coeff_perms.sort(key=lambda x: (inv(x),*x))
			
				perm_pairs = []
				
				for perm in coeff_perms:
					downperm = mulperm(list(perm),inverse_kperm)
					if inv(downperm) == inv(perm) - inv_kperm:
						flag = True			
						for i in range(N):
							if downperm[i] > N:
								flag = False
								break
						if not flag:
							continue
						firstperm = downperm[0:N]
						secondperm = [downperm[i]-N for i in range(N,len(downperm))]
						perm_pairs += [[permtrim(firstperm),permtrim(secondperm)]]
			
				if ascode:
					width = max([len(str(trimcode(perm[0]))+" "+str(trimcode(perm[1]))) for perm in perm_pairs])
				else:
					width = max([len(str(perm[0])+" "+str(perm[1])) for perm in perm_pairs])
							
				for perm in coeff_perms:
					val = coeff_dict[perm]
					downperm = mulperm(list(perm),inverse_kperm)
					if inv(downperm) == inv(perm) - inv_kperm:
						flag = True			
						for i in range(N):
							if downperm[i] > N:
								flag = False
								break
						if not flag:
							continue
						firstperm = downperm[0:N]
						secondperm = [downperm[i]-N for i in range(N,len(downperm))]
						val = sympify(val).subs(subs_dict)
											
						if val != 0:
							if display_positive:												
								if expand(val) != 0:
									val2 = compute_positive_rep(val,var2neg,var3neg,msg)
									if val2 is not None:
										val = val2
								else:
									val = 0
							if val != 0:
								if not ascode:
									width2 = width - len(str(permtrim(firstperm))) - len(str(permtrim(secondperm)))
									print(f"{tuple(permtrim(firstperm))}{' ':>{width2}}{tuple(permtrim(secondperm))}  {str(val).replace('**','^').replace('*',' ')}")
								else:
									width2 = width - len(str(trimcode(firstperm))) - len(str(trimcode(secondperm)))
									print(f"{trimcode(firstperm)}{' ':>{width2}}{trimcode(secondperm)}  {str(val).replace('**','^').replace('*',' ')}")
		else:
			if ascode:
				for i in range(len(perms)):
					perms[i] = tuple(permtrim(uncode(perms[i])))
			else:
				for i in range(len(perms)):
					perms[i] = tuple(permtrim([*perms[i]]))
			
			size = 0
			L = len(perms)
			orig_perms = [*perms]
			while len(perms) != size:
				size = len(perms)
				perms = split_perms(perms)
			
			coeff_dict = {perms[0]: 1}
			check_coeff_dict = {perms[0]: 1}
			for perm in orig_perms[1:]:
				check_coeff_dict = schubmult(check_coeff_dict,perm)
			if display_positive and len(perms)==2 and will_formula_work(perms[0],perms[1]) and perms[0]:
				coeff_dict = {}
				th = theta(perms[1])
				muv = uncode(th)
				muvn1v = mulperm(inverse(muv),perms[1])
				coeff_dict2 = {perms[0]: 1}
				coeff_dict2 = schubmult(coeff_dict2,muv)
				
				for perm, val in coeff_dict2.items():		
					w = mulperm(list(perm),muvn1v)
					if inv(w)+inv(muvn1v) == inv(perm):
						coeff_dict[tuple(permtrim(w))] = val
				posified = True
				
			if display_positive and len(perms)>2:
				coeff_dict2 = dict(coeff_dict)
				for perm in perms[1:]:					
					coeff_dict3 = {}
					for u in coeff_dict2:
						coeff_dict4 = {u: 1}
						coeff_dict4 = schubmult(coeff_dict4,perm)
						for w in coeff_dict4:
							coeff_dict4[w] = coeff_dict2[u]*posify(coeff_dict4[w],u,perm,w,var2,var3,msg)
						coeff_dict3 = add_perm_dict(coeff_dict4,coeff_dict3)					
					coeff_dict2 = coeff_dict3
				coeff_dict = coeff_dict2
				posified = True
			else:				
				coeff_dict = check_coeff_dict
					
			if pr:
				if ascode:
					width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys()])
				else:
					width = max([len(str(perm)) for perm in coeff_dict.keys()])
				
				coeff_perms = list(coeff_dict.keys())
				coeff_perms.sort(key=lambda x: (inv(x),*x))
				
				for perm in coeff_perms:
					val = coeff_dict[perm]
					if val != 0:
						notint = False
						try:
							int(val)
						except Exception:
							notint = True
						if notint and display_positive:
							valu = val
							if len(perms) == 2 and not posified:
								val = posify(val,perms[0],perms[1],perm,var2,var3,msg)							
							elif not posified:								
								val2 = compute_positive_rep(val,var2,var3,msg)
								if val2 is not None:
									val2 = val
							if check and expand(val - check_coeff_dict.get(perm,0))!=0:
								print(f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {check_coeff_dict.get(perm,0)=}")
								print(f"{perms=},{perm=}")
								exit(1)
						if val!=0:
							if ascode:
								print(f"{str(trimcode(perm)):>{width}}  {str(val).replace('**','^').replace('*',' ')}")	
							else:
								print(f"{str(perm):>{width}}  {str(val).replace('**','^').replace('*',' ')}")	
	except BrokenPipeError:
		pass
		
if __name__ == "__main__":
	main()
