from symengine import *
import sys
from bisect import bisect_left
from functools import cache
from itertools import chain
from cachetools import cached
from cachetools.keys import hashkey
from schubmult.perm_lib import *
import numpy as np
import pulp as pu
import itertools as it
import sympy
import psutil
from sortedcontainers import SortedList 

def count_sorted(mn,tp):
	index = bisect_left(mn,tp)
	ct = 0
	if mn[index] == tp:
		while index<len(mn) and mn[index] == tp:
			ct += 1
	return ct

n = 100

var = tuple(symarray('x',n).tolist())
var2 = tuple(symarray('y',n).tolist())
var3 = tuple(symarray('z',n).tolist())


def forwardcoeff(u,v,perm,var2=var2,var3=var3):
	th = theta(v)
	muv = uncode(th)
	vmun1 = mulperm(inverse([*v]),muv)
	
	w = mulperm([*perm],vmun1)	
	if inv(w) == inv(vmun1) + inv(perm):
		coeff_dict = schubmult_one(tuple(permtrim([*u])),tuple(muv),var2,var3)
		return coeff_dict.get(tuple(permtrim(w)),0)
	return 0

def dualcoeff(u,v,perm,var2=var2,var3=var3):	
	if u == (1,2):
		vp = mulperm([*v],inverse(perm))
		if inv(vp) == inv(v) - inv(perm):
			val = schubpoly(vp,var2,var3)
		else:
			val = 0
	else:
		dpret = []
		if dominates(u,perm):
			dpret = dualpieri([*u],[*v],[*perm])		
		else:
			th = theta(u)
			muu = uncode(th)
			umun1 = mulperm(inverse([*u]),muu)
			w = mulperm([*perm],umun1)
			if inv(w) == inv(umun1) + inv(perm):
				dpret = dualpieri(muu,[*v],w)
		ret = 0
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
			vpl = divdiffable(vp,cycle(lm[i]+1,cn1w[i]-lm[i]))
			if vpl == []:
				continue
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

@cache
def schubmult_one(perm1,perm2,var2=var2,var3=var3):
	return schubmult({perm1: 1},perm2,var2,var3)

def schubmult(perm_dict,v,var2=var2,var3=var3):
	vn1 = inverse(v)
	th = theta(vn1)
	if th[0]==0:
		return perm_dict		
	mu = permtrim(uncode(th))
	vmu = permtrim(mulperm([*v],mu))
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
						sumval = vpathsums[up].get(v,zero)
						if sumval == 0:
							continue
						for v2,vdiff,s in vpathdicts[index][v]:
							newpathsums[up2][v2] = newpathsums[up2].get(v2,zero)+s*sumval*elem_sym_func(th[index],index+1,up,up2,v,v2,udiff,vdiff,var2,var3)
			vpathsums = newpathsums
		toget = tuple(vmu)
		ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget,0) for ep in vpathsums},ret_dict)
	return ret_dict

fvar = 0

def poly_to_vec(poly,vec0=None):	
	global dimen, monom_to_vec, base_vec
	poly = expand(poly.xreplace({var3[1]: 0}))
	
	dc = poly.as_coefficients_dict()
	
	if vec0 is None:
		init_basevec(dc)
	
	vec = {}
	for mn in dc:
		cf = dc[mn]		
		if cf == 0:
			continue
		cf = abs(int(cf))
		try:
			index = monom_to_vec[mn]
		except KeyError:
			return None		
		if vec0 is not None and vec0[index]<cf:
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
	global dimen, monom_to_vec, base_vec
	monom_to_vec = {}
	index = 0
	for mn in dc:
		if dc[mn] == 0:
			continue
		monom_to_vec[mn] = index
		index+=1	
	dimen = index
	base_vec = [0 for i in range(dimen)]
	

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


def fres(v):
	for s in v.free_symbols:
		return s

def split_mul(arg0,var2=var2,var3=var3):
	monoms = SortedList()
	
	
	
	var2s = {fres(var2[i]): i for i in range(len(var2))}
	var3s = {fres(var3[i]): i for i in range(len(var3))}
	#print(f"{type(arg0)=} {arg0=}")
	if isinstance(arg0,Pow):
		arg = arg0
		arg2 = expand(arg.args[0])
		yval = arg2.args[0]
		zval = arg2.args[1]
		if str(yval).find("z")!=-1:
			yval, zval = zval, yval
		if str(zval).find("-") != -1:
			zval = -zval
		if str(yval).find("-") != -1:
			yval = -yval
		tup = (var2s[fres(yval)],var3s[fres(zval)])
		for i in range(int(arg0.args[1])):
			monoms += [tup]
	else:
		for arg in arg0.args:
			if is_flat_term(arg):
				if isinstance(arg,Integer) or isinstance(arg,int):
					continue
				arg = expand(arg)
				if arg == 0:
					break
				yval = arg.args[0]
				zval = arg.args[1]
				if str(yval).find("z")!=-1:
					yval, zval = zval, yval
				if str(zval).find("-") != -1:
					zval = -zval
				if str(yval).find("-") != -1:
					yval = -yval
				monoms += [(var2s[fres(yval)],var3s[fres(zval)])]
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
				tup = (var2s[fres(yval)],var3s[fres(zval)])
				for i in range(int(arg.args[1])):
					monoms += [tup]
	return monoms

def split_monoms(pos_part,var2,var3):
	arrs = SortedList()	
	if isinstance(pos_part,Add):		
		for arg0 in pos_part.args:
			monoms = split_mul(arg0,var2,var3)			
			arrs += [monoms]
	elif isinstance(pos_part,Mul) or isinstance(pos_part,Pow):
		monoms = split_mul(pos_part,var2,var3)
		arrs += [monoms]		
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

def find_base_vectors(monom_list,monom_list_neg,var2,var3,depth):
	size = 0
	mn_fullcount = {}
	#pairs_checked = set()
	monom_list = set([tuple(mn) for mn in monom_list])
	additional_set = set(monom_list)
	ct = 0
	while ct < depth and size!=len(monom_list):
		size = len(monom_list)
		#found = False
		#for mn in mons2:
		#	if mn not in monom_list:
		#		found = True
		#		break
		#if not found:
		#	print("Breaking")
		#	break
		
		monom_list2 = set(monom_list)
		additional_set2 = set()
		for mn in monom_list:
		
			#res = 1
			#for tp in mn:
			#	res *= var2[tp[0]] - var3[tp[1]]
			#if poly_to_vec(res,vec) is None:
			#	continue
			
			mncount = mn_fullcount.get(mn,{})
			if mncount == {}:
				for tp in mn:
					mncount[tp] = mncount.get(tp,0) + 1
				mn_fullcount[mn] = mncount
			for mn2 in monom_list:
				#if (mn,mn2) in pairs_checked:
				#	continue
				mn2count = mn_fullcount.get(mn2,{})
				if mn2count == {}:
					for tp in mn2:
						mn2count[tp] = mn2count.get(tp,0) + 1
					mn_fullcount[mn2] = mn2count
				num_diff = 0				
				for tp in mncount:
					pt = mn2count.get(tp,0) - mncount[tp]
					num_diff += abs(pt)
					if num_diff > 1:						
						break											
				if num_diff == 1:
					diff_term1 = None
					diff_term2 = None
					for tp in mn2count:
						if mn2count[tp] > mncount.get(tp,0):
							diff_term2 = tp
							break
					for tp2 in mncount:
						if mncount[tp2]>mn2count.get(tp2,0):
							diff_term1 = tp2
							break
					#print(f"{mn,mn2}")
					if diff_term1 is None or diff_term2 is None:
						print(f"{mn=} {mn2=}")
						exit(1)
					if diff_term2[1] == diff_term1[1]:
						continue
					new_term1 = (diff_term1[0],diff_term2[1])
					new_term2 = (diff_term2[0],diff_term1[1])
					#mn3 = [*mn]
					#mn4 = list(mn2)
					index = bisect_left(mn,diff_term1)
					mn3 = list(mn[:index])+list(mn[index+1:])
					index = bisect_left(mn3,new_term1)
					mn3_t = tuple(mn3[:index]+[new_term1]+mn3[index:])
					index2 = bisect_left(mn2,diff_term2)
					mn4 = list(mn2[:index2])+list(mn2[index2+1:])
					index2 = bisect_left(mn4,new_term2)
					mn4_t = tuple(mn4[:index2]+[new_term2]+mn4[index2:])
					#res = 1
					#for tp in mn3_t:
					#	res *= var2[tp[0]] - var3[tp[1]]
					#if poly_to_vec(res,vec) is not None:
					if mn3_t not in monom_list2:
						additional_set2.add(mn3_t)
					monom_list2.add(mn3_t)
					#res = 1
					#for tp in mn4_t:
					#	res *= var2[tp[0]] - var3[tp[1]]
					##
					##	additional_set2.add(mn3_t)
					#if poly_to_vec(res,vec) is not None:
					if mn4_t not in monom_list2:
						additional_set2.add(mn4_t)
					monom_list2.add(mn4_t)
		additional_set = additional_set2
		#print(f"{len(additional_set)=} {len(monom_list2)=}")
		monom_list = monom_list2		
		ct+=1
	#print(f"{monom_list=}")
	ret = []
	#print("done monom")
	for mn in monom_list:
		if len(mn)!=len(set(mn)):
			continue
		res = 1
		for tp in mn:
			res *= var2[tp[0]] - var3[tp[1]]
		ret+=[res]
	return ret, monom_list

def compute_positive_rep(val,var2=var2,var3=var3,msg=False,do_pos_neg=True):
	notint = False	
	try:
		int(expand(val))
		val2 = expand(val)
	except Exception:
		notint = True
	if notint:
	
		
		frees = val.free_symbols
		var2list = [*var2]
		var3list = [*var3]
		
		for i in range(len(var2list)):
			symset = var2list[i].free_symbols
			for sym in symset:
				var2list[i] = sym
				
		for i in range(len(var3list)):
			symset = var3list[i].free_symbols
			for sym in symset:
				var3list[i] = sym
		
		varsimp2 = [m for m in frees if m in var2list]		
		varsimp3 = [m for m in frees if m in var3list]
		varsimp2.sort(key=lambda k: var2list.index(k))
		varsimp3.sort(key=lambda k: var3list.index(k))
				
		
		var22 = [sympy.sympify(m) for m in varsimp2]
		var33 = [sympy.sympify(m) for m in varsimp3]
		n1 = len(varsimp2)
		n2 = len(varsimp3)		
		
		for i in range(len(varsimp2)):
			varsimp2[i] = var2[var2list.index(varsimp2[i])]
		for i in range(len(varsimp3)):
			varsimp3[i] = var3[var3list.index(varsimp3[i])]
		
		base_vectors = []
		base_monoms = []						
		vec = poly_to_vec(val,None)
		
		if do_pos_neg:
			smp = val
			flat, found_one = flatten_factors(smp)
			while found_one:
				flat, found_one = flatten_factors(flat,varsimp2,varsimp3)
			pos_part = 0
			neg_part = 0
			if isinstance(flat,Add) and not is_flat_term(flat):	
				for arg in flat.args:
					if expand(arg) == 0:
						continue
					if not is_negative(arg):
						pos_part += arg
					else:
						neg_part -= arg
				if neg_part == 0:
					#print("no neg")
					return pos_part
			depth = 1
					
			mons = split_monoms(pos_part,varsimp2,varsimp3)
			mons = set([tuple(mn) for mn in mons])
			mons2 = split_monoms(neg_part,varsimp2,varsimp3)
			mons2 = set([tuple(mn2) for mn2 in mons2])
			
			#mons2 = split_monoms(neg_part)
			#for mn in mons2:
			#	if mn not in mons:
			#		mons.add(mn)
			#print(mons)
			status = 0
			size = len(mons)
			while status != 1:
				base_monoms, mons = find_base_vectors(mons,mons2,varsimp2,varsimp3,depth)
				if len(mons) == size:
					print(f"Found counterexample")
					return -1
				size = len(mons)
				base_vectors = []
				bad = False
				bad_vectors = []
				for i in range(len(base_monoms)):
					vec0 = poly_to_vec(base_monoms[i],vec)
					if vec0 is not None:
						base_vectors += [vec0]
					else:
						bad_vectors += [i]
				for j in range(len(bad_vectors)-1,-1,-1):
					base_monoms.pop(bad_vectors[j])
			
				vrs = [pu.LpVariable(name=f"a{i}",lowBound=0,cat='Integer') for i in range(len(base_vectors))]
				lp_prob = pu.LpProblem('Problem',pu.LpMinimize)		
				lp_prob += int(0)
				eqs = [*base_vec]
				for j in range(len(base_vectors)):
					for i in base_vectors[j]:
						bvi = base_vectors[j][i]
						if bvi == 1:
							eqs[i] += vrs[j]
						else:
							eqs[i] += bvi*vrs[j]
				for i in range(dimen):
					try:
						lp_prob += eqs[i] == vec[i]
					except TypeError:
						bad = True
						break
				if bad:
					continue
				try:
					solver = pu.PULP_CBC_CMD(msg=msg)
					status = lp_prob.solve(solver)
				except KeyboardInterrupt:			
					current_process = psutil.Process()
					children = current_process.children(recursive=True)
					for child in children:
						child_process = psutil.Process(child.pid)
						child_process.terminate()
						child_process.kill()
					raise KeyboardInterrupt()
				status = lp_prob.status
		else:
			val_poly = sympy.poly(expand(val),*var22,*var33)
			vec = poly_to_vec(val)
			mn = val_poly.monoms()
			L1 = tuple([0 for i in range(n1)])
			mn1L = []
			lookup = {}
			for mm0 in mn:
				key = mm0[n1:]
				if key not in lookup:
					lookup[key] = []
				mm0n1 = mm0[:n1]
				st = set(mm0n1)
				if len(st.intersection(set([0,1])))==len(st) and 1 in st:
					lookup[key]+=[mm0]
				if mm0n1 == L1:
					mn1L += [mm0]
			for mn1 in mn1L:
				comblistmn1 = [1]
				for i in range(n1,len(mn1)):
					if mn1[i]!=0:
						arr = np.array(comblistmn1)
						comblistmn12 = []
						mn1_2 = tuple([*mn1[n1:i]]+[0]+[*mn1[i+1:]])
						for mm0 in lookup[mn1_2]:												
							comblistmn12 += (arr*np.prod([varsimp2[k] - varsimp3[i - n1] for k in range(n1) if mm0[k]==1])).tolist()						
						comblistmn1 = comblistmn12
				for i in range(len(comblistmn1)):
					b1 = comblistmn1[i]
					vec0 = poly_to_vec(b1,vec)
					if vec0 is not None:
						base_vectors += [vec0]
						base_monoms += [b1]
			vrs = [pu.LpVariable(name=f"a{i}",lowBound=0,cat='Integer') for i in range(len(base_vectors))]
			lp_prob = pu.LpProblem('Problem',pu.LpMinimize)		
			lp_prob += int(0)
			eqs = [*base_vec]
			for j in range(len(base_vectors)):
				for i in base_vectors[j]:
					bvi = base_vectors[j][i]
					if bvi == 1:
						eqs[i] += vrs[j]
					else:
						eqs[i] += bvi*vrs[j]
			for i in range(dimen):
				lp_prob += eqs[i] == vec[i]
			try:
				solver = pu.PULP_CBC_CMD(msg=msg)
				status = lp_prob.solve(solver)
			except KeyboardInterrupt:			
				current_process = psutil.Process()
				children = current_process.children(recursive=True)
				for child in children:
					child_process = psutil.Process(child.pid)
					child_process.terminate()
					child_process.kill()
				raise KeyboardInterrupt()
		#print(f"{pos_part=}")
		#print(f"{neg_part=}")
		#else:
		#	print(f"No dice {flat=}")
		#	exit(1)
		#	#val = pos_part - neg_part
		
			#depth+=1
		val2 = 0
		for k in range(len(base_vectors)):
			x = vrs[k].value()
			b1 = base_monoms[k]
			if x!=0 and x is not None:
				val2 += int(x)*b1
	return val2

def is_split_two(u,v,w):
	if inv(w)-inv(u)!=2:
		return False, []
	diff_perm = mulperm(inverse([*u]),[*w])
	identity = [i+1 for i in range(len(diff_perm))]
	cycles = []
	for i in range(len(identity)):
		if diff_perm[i] != identity[i]:
			cycle0 = set()
			cycle = set([i+1])
			last = i
			while len(cycle0) != len(cycle):
				cycle0 = cycle
				last = diff_perm[last] - 1
				cycle.add(last+1)
			if len(cycle)>1 and cycle not in cycles:
				cycles += [cycle]
			if len(cycles)>2:
				break
	if len(cycles) == 2:
		return True, cycles
	else:
		return False, []

def is_coeff_irreducible(u,v,w):
	return not will_formula_work(u,v) and not will_formula_work(v,u) and not one_dominates(u,w) and not is_reducible(v) and inv(w) - inv(u)>1 and not is_split_two(u,v,w)[0] and len([i for i in code(v) if i !=0]) > 1

def is_hook(cd):
	started = False
	done = False
	found_zero_after = False
	for i in range(len(cd)):
		if (done or found_zero_after) and cd[i]!=0:
			return False
		if cd[i] == 1 and not started:
			started = True
		if cd[i] > 1:
			done = True
		if started and cd[i] == 0:
			found_zero_after = True
	if started or done:
		return True
	return False

def div_diff(i,poly):
	return sympify(sympy.div(sympy.sympify(poly - permy(poly,i)),sympy.sympify(var2[i]-var2[i+1]))[0])

def skew_div_diff(u,w,poly):
	d = -1
	for i in range(len(w)-1):
		if w[i]>w[i+1]:
			d = i
			break
	d2 = -1
	for i in range(len(u)-1):
		if u[i]>u[i+1]:
			d2 = i
			break
	if d == -1:
		if d2 == -1:
			return poly
		return 0
	w2 = [*w]
	w2[d], w2[d+1] = w2[d+1], w2[d]
	if d<len(u)-1 and u[d]>u[d+1]:		
		u2 = [*u]		
		u2[d], u2[d+1] = u2[d+1], u2[d]
		return skew_div_diff(u2,w2,permy(poly,d+1))
	else:
		return skew_div_diff(u,w2,div_diff(d+1,poly))

@cached(cache={}, key=lambda val, u2,v2,w2,var2=var2,var3=var3,msg=False,do_pos_neg=True,sign_only=False: hashkey(u2,v2,w2,var2,var3,msg,do_pos_neg,sign_only))
def posify(val,u2,v2,w2,var2=var2,var3=var3,msg=False,do_pos_neg=True,sign_only=False):
	if inv(u2)+inv(v2) - inv(w2) == 0:
		return val
	cdv = code(v2)	
	if set(cdv) == set([0,1]) and do_pos_neg:	
		return val
	#if is_hook(cdv):
	#	print(f"Could've {cdv}")
	if not sign_only and expand(val) == 0:
		return 0
	
	u, v, w = try_reduce_v(u2, v2, w2)
	if is_coeff_irreducible(u,v,w):
		u, v, w = try_reduce_u(u2, v2, w2)
		if is_coeff_irreducible(u,v,w):
			u, v, w = [*u2], [*v2], [*w2]
			if is_coeff_irreducible(u,v,w):
				w0 = [*w]
				u, v, w = reduce_descents(u,v,w)
				if is_coeff_irreducible(u,v,w):
					u, v, w = reduce_coeff(u,v,w)
					if is_coeff_irreducible(u,v,w):
						while is_coeff_irreducible(u,v,w) and tuple(permtrim(w0))!=tuple(permtrim([*w])):
							w0 = w
							u, v, w = reduce_descents(u,v,w)						
							if is_coeff_irreducible(u,v,w):	
								u, v, w = reduce_coeff(u,v,w)
	u = tuple(u)
	v = tuple(v)	
	w = tuple(w)
	
	if w != w2 and sign_only:
		return 0
			
	if is_coeff_irreducible(u,v,w):
		u3, v3, w3 = try_reduce_v(u, v, w)
		if not is_coeff_irreducible(u3,v3,w3):
			u, v, w = u3, v3, w3
		else:
			u3, v3, w3 = try_reduce_u(u, v, w)		
			if not is_coeff_irreducible(u3,v3,w3):
				u, v, w = u3, v3, w3
	split_two_b, split_two = is_split_two(u,v,w)
	
	if len([i for i in code(v) if i !=0]) == 1:
		if sign_only:
			return 0
		cv = code(v)
		for i in range(len(cv)):
			if cv[i]!=0:
				k = i+1
				p = cv[i]
				break					
		inv_u = inv(u)
		r = inv(w) - inv_u
		val = 0
		w2 = w
		hvarset = [w2[i] for i in range(min(len(w2),k))]+[i+1 for i in range(len(w2),k)] + [w2[b] for b in range(k,len(u)) if u[b]!=w2[b]]+[w2[b] for b in range(len(u),len(w2))]
		val = elem_sym_poly(p-r,k+p-1,[-var3[i] for i in range(1,n)],[-var2[i] for i in hvarset])	
	elif (will_formula_work(v,u) or dominates(u,w)):
		if sign_only:
			return 0
		val = dualcoeff(u,v,w,var2,var3)
	elif inv(w) - inv(u) == 1:	
		if sign_only:
			return 0
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
	elif split_two_b:
		if sign_only:
			return 0
		cycles = split_two
		a1, b1 = cycles[0]
		a2, b2 = cycles[1]
		a1-=1
		b1-=1
		a2-=1
		b2-=1
		spo = sorted([a1,b1,a2,b2])
		real_a1 = min(spo.index(a1),spo.index(b1))
		real_a2 = min(spo.index(a2),spo.index(b2))
		real_b1 = max(spo.index(a1),spo.index(b1))
		real_b2 = max(spo.index(a2),spo.index(b2))
		
		good1 = False
		good2 = False
		if real_b1 - real_a1 == 1:
			good1 = True
		if real_b2 - real_a2 == 1:
			good2 = True
		a, b = -1, -1
		if good1 and not good2:
			a, b = min(a2,b2), max(a2,b2)
		if good2 and not good1:
			a, b = min(a1,b1), max(a1,b1)
		arr = [[[],v]]
		d = -1
		for i in range(len(v)-1):
			if v[i]>v[i+1]:
				d = i + 1
		for i in range(d):
			arr2 = []
						
			if i in [a1,b1,a2,b2]:
				continue
			i2 = 1
			i2 += len([aa for aa in [a1,b1,a2,b2] if i>aa])
			for vr, v2 in arr:
				dpret = pull_out_var(i2,[*v2])
				for v3r, v3 in dpret:
					arr2 += [[vr + [(v3r,i+1)],v3]]
			arr = arr2
		val = 0
		
		if good1:
			arr2 = []
			for L in arr:
				v3 = [*L[-1]]
				if v3[real_a1]<v3[real_b1]:
					continue
				else:
					v3[real_a1], v3[real_b1] = v3[real_b1], v3[real_a1]
				arr2 += [[L[0],v3]]
			arr = arr2
			if not good2:
				for i in range(4):
					arr2 = []
					
					if i in [real_a2,real_b2]:
						continue
					if i == real_a1:
						var_index = min(a1,b1) + 1
					elif i == real_b1:
						var_index = max(a1,b1) + 1
					i2 = 1
					i2 += len([aa for aa in [real_a2,real_b2] if i>aa])
					for vr, v2 in arr:
						dpret = pull_out_var(i2,[*v2])
						for v3r, v3 in dpret:
							arr2 += [[vr + [(v3r,var_index)],v3]]
					arr = arr2
		if good2:
			arr2 = []
			for L in arr:
				v3 = [*L[-1]]
				try:
					if v3[real_a2]<v3[real_b2]:
						continue
					else:
						v3[real_a2], v3[real_b2] = v3[real_b2], v3[real_a2]
				except IndexError:
					continue
				arr2 += [[L[0],v3]]
			arr = arr2
			if not good1:
				for i in range(4):
					arr2 = []
								
					if i in [real_a1,real_b1]:
						continue
					i2 = 1
					i2 += len([aa for aa in [real_a1,real_b1] if i>aa])
					if i == real_a2:
						var_index = min(a2,b2) + 1
					elif i == real_b2:
						var_index = max(a2,b2) + 1
					for vr, v2 in arr:
						dpret = pull_out_var(i2,[*v2])
						for v3r, v3 in dpret:
							arr2 += [[vr + [(v3r,var_index)],v3]]
					arr = arr2
		
		for L in arr:
			v3 = [*L[-1]]
			tomul = 1
			doschubpoly = True
			if (not good1 or not good2) and v3[0]<v3[1] and (good1 or good2):
				continue
			elif (good1 or good2) and (not good1 or not good2):
				v3[0], v3[1] = v3[1], v3[0]
			elif not good1 and not good2:
				doschubpoly = False
				if v3[0] < v3[1]:
					dual_u = uncode([2,0])
					dual_w = [4,2,1,3]
					coeff = permy(dualcoeff(dual_u,v3,dual_w,var2,var3),2)
					
				elif len(v3)<3 or v3[1] < v3[2]:
					if len(v3)<=3 or v3[2]<v3[3]:
						coeff = 0
						continue
					else:
						v3[0], v3[1] = v3[1], v3[0]
						v3[2], v3[3] = v3[3], v3[2]
						coeff = permy(schubpoly(v3,var2,var3),2)
				elif len(v3)<=3 or v3[2] < v3[3]:
					if len(v3)<=3:
						v3 += [4]
					v3[2], v3[3] = v3[3], v3[2]
					coeff = permy(posify(schubmult_one((1,3,2),tuple(permtrim([*v3])),var2,var3).get((2,4,3,1),0),(1,3,2),tuple(permtrim([*v3])),(2,4,3,1),var2,var3,msg,do_pos_neg),2)					
				else:
					coeff = permy(schubmult_one((1,3,2),tuple(permtrim([*v3])),var2,var3).get((2,4,1,3),0),2)
				tomul = sympify(coeff)
			toadd = 1
			for i in range(len(L[0])):
				var_index = L[0][i][1]
				oaf = L[0][i][0]
				if var_index-1>=len(w):
					yv = var_index
				else:
					yv = w[var_index-1]
				for j in range(len(oaf)):				
					toadd*= var2[yv] - var3[oaf[j]]
			if (not good1 or not good2) and (good1 or good2):
				varo = [0,var2[w[a]],var2[w[b]]]
			else:
				varo = [0,*[var2[w[spo[k]]] for k in range(4)]]
			if doschubpoly:
				toadd *= schubpoly(v3,varo,var3)
			else:
				subs_dict3 = {var2[i]: varo[i] for i in range(len(varo))}
				toadd*=tomul.subs(subs_dict3)
			val += toadd
	elif will_formula_work(u,v):
		if sign_only:
			return 0
		val = forwardcoeff(u,v,w,var2,var3)
	#elif inv(w) - inv(u) == 2:
	#	indices = []
	#	for i in range(len(w)):
	#		if i>=len(u) or u[i]!=w[i]:
	#			indices += [i+1]
	#	arr = [[[],v]]
	#	d = -1
	#	for i in range(len(v)-1):
	#		if v[i]>v[i+1]:
	#			d = i + 1
	#	for i in range(d):
	#		arr2 = []
	#					
	#		if i+1 in indices:
	#			continue
	#		i2 = 1
	#		i2 += len([aa for aa in indices if i+1>aa])
	#		for vr, v2 in arr:
	#			dpret = pull_out_var(i2,[*v2])
	#			for v3r, v3 in dpret:
	#				arr2 += [[vr + [(v3r,i+1)],v3]]
	#		arr = arr2
	#	val = 0
	#	
	#	for L in arr:
	#		v3 = [*L[-1]]
	#		tomul = 1
	#		pooly = skew_div_diff(u,w,schubpoly(v3,[0,*[var2[a] for a in indices]],var3))
	#		coeff = compute_positive_rep(pooly,var2,var3,msg,False)
	#		if coeff == -1:
	#			return -1
	#		tomul = sympify(coeff)
	#		toadd = 1
	#		for i in range(len(L[0])):
	#			var_index = L[0][i][1]
	#			oaf = L[0][i][0]
	#			if var_index-1>=len(w):
	#				yv = var_index
	#			else:
	#				yv = w[var_index-1]
	#			for j in range(len(oaf)):				
	#				toadd*= var2[yv] - var3[oaf[j]]
	#		toadd*=tomul#.subs(subs_dict3)
	#		val += toadd		
	else:
		c01 = code(u)
		c02 = code(w)
		c03 = code(v)
		
		c1 = code(inverse(u))
		c2 = code(inverse(w))
		
		if one_dominates(u,w):
			if sign_only:
				return 0
			while c1[0] != c2[0]:				
				w = [*w]
				v = [*v]
				w[c2[0]-1], w[c2[0]] = w[c2[0]], w[c2[0]-1]
				v[c2[0]-1], v[c2[0]] = v[c2[0]], v[c2[0]-1]
				w = tuple(w)
				v = tuple(v)
				c2 = code(inverse(w))
				c03 = code(v)
				c01 = code(u)
				c02 = code(w)				
				
		
		if is_reducible(v):
			if sign_only:
				return 0
			newc = []
			elemc = []
			for i in range(len(c03)):
				if c03[i]>0:
					newc += [c03[i]-1]
					elemc += [1]
				else:
					break
			v3 = uncode(newc)
			coeff_dict = schubmult_one(tuple(permtrim([*u])),tuple(permtrim(uncode(elemc))),var2,var3)
			val = 0
			for new_w in coeff_dict:
				tomul = coeff_dict[new_w]
				newval = schubmult_one(new_w,tuple(permtrim(uncode(newc))),var2,var3).get(tuple(permtrim([*w])),0)
				newval = posify(newval,new_w,tuple(permtrim(uncode(newc))),w,var2,var3,msg,do_pos_neg)
				val += tomul*shiftsubz(newval)
		elif c01[0] == c02[0] and c01[0] != 0:
			if sign_only:
				return 0
			varl = c01[0]
			u3 = uncode([0] + c01[1:])
			w3 = uncode([0] + c02[1:])
			val = 0
			val = schubmult_one(tuple(permtrim(u3)),tuple(permtrim([*v])),var2,var3).get(tuple(permtrim(w3)),0)
			val = posify(val,tuple(permtrim(u3)),tuple(permtrim([*v])),tuple(permtrim(w3)),var2,var3,msg,do_pos_neg)
			for i in range(varl):
				val = permy(val,i+1)			
		elif c1[0] == c2[0]:
			if sign_only:
				return 0
			vp = pull_out_var(c1[0]+1,[*v])
			u3 = tuple(permtrim(phi1(u)))
			w3 = tuple(permtrim(phi1(w)))
			c3 = code(inverse(u3))
			c4 = code(inverse(w3))
			val = 0
			for arr, v3 in vp:
				tomul = 1
				for i in range(len(arr)):
					tomul*=var2[1] - var3[arr[i]]
				
				val2 = schubmult_one(tuple(permtrim(u3)),tuple(permtrim(v3)),var2,var3).get(tuple(permtrim(w3)),0)
				val2 = posify(val2,u3,tuple(permtrim(v3)),w3,var2,var3,msg,do_pos_neg)
				val += tomul*shiftsub(val2)
		elif inv(w)-inv(u)==2:
			indices = []
			for i in range(len(w)):
				if i>=len(u) or u[i]!=w[i]:
					indices += [i+1]
			arr = [[[],v]]
			d = -1
			for i in range(len(v)-1):
				if v[i]>v[i+1]:
					d = i + 1
			for i in range(d):
				arr2 = []
							
				if i+1 in indices:
					continue
				i2 = 1
				i2 += len([aa for aa in indices if i+1>aa])
				for vr, v2 in arr:
					dpret = pull_out_var(i2,[*v2])
					for v3r, v3 in dpret:
						arr2 += [[vr + [(v3r,i+1)],v3]]
				arr = arr2
			val = 0
			
			for L in arr:
				v3 = [*L[-1]]
				tomul = 1						
				toadd = 1
				for i in range(len(L[0])):
					var_index = L[0][i][1]
					oaf = L[0][i][0]
					if var_index-1>=len(w):
						yv = var_index
					else:
						yv = w[var_index-1]
					for j in range(len(oaf)):				
						toadd*= var2[yv] - var3[oaf[j]]
				pooly = skew_div_diff(u,w,schubpoly(v3,[0,*[var2[a] for a in indices]],var3))
				if toadd == 0:
					continue
				if pooly !=0:
					coeff = compute_positive_rep(pooly,var2,var3,msg,False)
				else:
					coeff = 0
				if coeff == -1:
					return -1
				tomul = sympify(coeff)
				toadd*=tomul#.subs(subs_dict3)
				val += toadd						
		else:
			if not sign_only:
				if inv(u)+inv(v)-inv(w)==1:
					val2 = compute_positive_rep(val,var2,var3,msg,False)
				else:
					val2 = compute_positive_rep(val,var2,var3,msg,do_pos_neg)					
				if val2 is not None:				
					val = val2							
			else:
				#st = str(expand(val))
				#if st.find("-")!=-1:
				#	return -1
				#else:
				#	return val
				d = expand(val).as_coefficients_dict()
				for v in d.values():
					if v<0:
						return -1
				return 1
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
			print("**** schubmult_yz ****")
			print("Purpose: Compute products (and coproducts) of double Schubert polynomials in different sets of variables")
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
			pos2 = []
			last_descent = -1
			poso = []
			for i in range(len(perms[0])-1):
				if perms[0][i]>perms[0][i+1]:
					last_descent = i+1
			for i in range(1,last_descent+1):
				if i not in pos:
					pos2 += [i-1]
				else:
					poso += [i-1]
			
			mu_W = uncode(theta(inverse(perms[0])))				
	
			the_top_perm = tuple(permtrim(mulperm(list(perms[0]),mu_W)))
	
			muA = uncode(mu_A(code(mu_W),poso))
			muB = uncode(mu_A(code(mu_W),pos2))
			#print(f"{code(mu_W)=},{pos=},{poso=} {pos2=},{p_trans(code(mu_W))=},{p_trans(code(muA))=},{p_trans(code(muB))=},{code(muA)=},{code(muB)=}")
			while cd[-1] == 0:
				cd.pop()
			k = len(pos)
			n = len(perms[0])
			kcd = [pos[i]-i-1 for i in range(len(pos))] + [n+1-k for i in range(k,n)]
			N = len(kcd)
			
			var2neg = np.array([-var2[i] for i in range(len(var2))])
			var3neg = np.array([-var3[i] for i in range(len(var3))])
			
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
								if val != 0:
									val2 = posify(val,tuple(permtrim(mulperm(firstperm,muA))),tuple(permtrim(mulperm(secondperm,muB))),the_top_perm,tuple(var2neg.tolist()),tuple(var3neg.tolist()),msg,False)
									if expand(val - val2) != 0:
										print(f"error; write to schubmult@gmail.com with the case {perms=}\n{code(firstperm)=} {code(secondperm)=}\n{val2=}\n{val=}")
										print(f"{code(tuple(permtrim(mulperm(firstperm,muA))))=},{code(tuple(permtrim(mulperm(secondperm,muB))))=},{code(the_top_perm)=}\n{expand(val-val2)=}")
										exit(1)
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
			if display_positive and len(perms)==2 and will_formula_work(perms[0],perms[1]):
				coeff_dict = {}
				th = theta(perms[1])
				muv = uncode(th)
				muvn1v = mulperm(inverse(muv),perms[1])
				coeff_dict2 = {perms[0]: 1}
				coeff_dict2 = schubmult(coeff_dict2,muv)
				for perm, val in coeff_dict2.items():					
					w = mulperm([*perm],muvn1v)
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
			elif not posified:	
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
							try:
								if len(perms) == 2 and not posified:
									val = posify(val,perms[0],perms[1],perm,var2,var3,msg)							
								elif not posified:								
									val = compute_positive_rep(val,var2,var3,msg)
							except TypeError:
								print(f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {check_coeff_dict.get(perm,0)=}")
								exit(1)
							if check and expand(val - check_coeff_dict.get(perm,0))!=0:
								print(f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {check_coeff_dict.get(perm,0)=}")
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
