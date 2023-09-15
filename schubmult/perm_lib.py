from symengine import *
from functools import cache
from itertools import chain

def inv(perm):
	return inv_cache(tuple(perm))

@cache
def inv_cache(perm):
	res = 0	
	for i in range(len(perm)):
		for j in range(i+1,len(perm)):
			if perm[i]>perm[j]:
				res+=1
	return res

def mulperm(perm1,perm2):
	if len(perm1)<len(perm2):
		return [perm1[perm2[i]-1] if perm2[i]<=len(perm1) else perm2[i] for i in range(len(perm2))]
	else:
		return [perm1[perm2[i]-1] for i in range(len(perm2))]+perm1[len(perm2):]

def uncode(cd):
	cd2 = list(cd)
	if cd2 == []:
		return [1,2]
	max_required = max([cd2[i]+i for i in range(len(cd2))])
	cd2 += [0 for i in range(len(cd2),max_required)]
	fullperm = [i+1 for i in range(len(cd2)+1)]
	perm = []
	for i in range(len(cd2)):
		perm += [fullperm[cd2[i]]]
		fullperm.pop(cd2[i])
	perm += [fullperm[0]]
	return perm


def code(perm):
	ret = []
	for i in range(len(perm)-1):
		ret += [0]
		for j in range(i+1,len(perm)):
			if perm[i]>perm[j]:
				ret[-1] += 1
	return ret

def inverse(perm):
	retperm = [0 for i in range(len(perm))]
	for i in range(len(perm)):
		retperm[perm[i]-1] = i+1
	return retperm

def permtrim(perm):
	L = len(perm)
	while L > 2 and perm[-1] == L:
		L = perm.pop() - 1
	return perm
	
def has_bruhat_descent(perm,i,j):
	if perm[i]<perm[j]:
		return False
	for p in range(i+1,j):
		if perm[i]>perm[p] and perm[p]>perm[j]:
			return False
	return True
	
def has_bruhat_ascent(perm,i,j):
	if perm[i]>perm[j]:
		return False
	for p in range(i+1,j):
		if perm[i]<perm[p] and perm[p]<perm[j]:
			return False
	return True	

def elem_sym_perms(orig_perm,p,k):	
	total_list = [(orig_perm,0)]
	up_perm_list = [(orig_perm,1000)]
	for pp in range(p):
		perm_list = []
		for up_perm, last in up_perm_list:	
			up_perm2 = list(up_perm) + [len(up_perm)+1]
			if len(up_perm2) < k + 1:
				up_perm2 += [i+1 for i in range(len(up_perm2),k+2)]
			pos_list = [i for i in range(k) if (((i>=len(orig_perm) and up_perm2[i] == i+1) or (i<len(orig_perm) and orig_perm[i] == up_perm[i])))]
			for j in range(k,len(up_perm2)):
				if up_perm2[j]>=last:
					continue
				for i in pos_list:			
					if has_bruhat_ascent(up_perm2,i,j):
						up_perm2[i],up_perm2[j] = up_perm2[j],up_perm2[i]
						new_perm = permtrim(list(up_perm2))
						up_perm2[i],up_perm2[j] = up_perm2[j],up_perm2[i]
						new_perm_add = tuple(new_perm)
						perm_list += [(new_perm_add,up_perm2[j])]
						total_list+=[(new_perm_add,pp+1)]
		up_perm_list = perm_list
	return total_list

# perms and inversion diff
def kdown_perms(perm,monoperm,p,k):	
	inv_m = inv(monoperm)
	inv_p = inv(perm)
	full_perm_list = []
	
	if inv(mulperm(list(perm),monoperm)) == inv_m - inv_p:
		full_perm_list+= [(tuple(perm),0,1)]
	
	down_perm_list = [(perm,1)]
	if len(perm)<k:
		return full_perm_list
	a2 = k-1
	for pp in range(1,p+1):
		down_perm_list2 = []
		for perm2, s in down_perm_list:
			L = len(perm2)
			if L<k:
				continue	
			s2 = -s
			for b in chain(range(k-1),range(k,L)):
				if perm2[b] != perm[b]:
					continue
				if b < a2:
					i,j = b,a2
				else:
					i,j,s2 = a2,b,s
				if has_bruhat_descent(perm2,i,j):
					new_perm = [*perm2]
					new_perm[a2],new_perm[b] = new_perm[b],new_perm[a2]
					permtrim(new_perm)
					down_perm_list2 += [(new_perm,s2)]
					if inv(mulperm(new_perm,monoperm)) == inv_m - inv_p + pp:
						full_perm_list += [(tuple(new_perm),pp,s2)]
		down_perm_list = down_perm_list2
	return full_perm_list
	
def compute_vpathdicts(th,vmu):
	vpathdicts = [{} for index in range(len(th))]
	vpathdicts[-1][tuple(vmu)] = None
	S = sum(th)
	thL = len(th)
	th2 = [*th]
	top = code(inverse(uncode(th)))
	for i in range(thL-1,-1,-1):
		top2 = code(inverse(uncode(top)))
		while top2[-1]==0:
			top2.pop()
		top2.pop()		
		top = code(inverse(uncode(top2)))
		monoperm = uncode(top)
		if len(monoperm) < 2:
			monoperm = [1,2]
		k = i+1
		for last_perm in vpathdicts[i]:
			newperms = kdown_perms(last_perm,monoperm,th[i],k)
			vpathdicts[i][last_perm]=newperms
			if i>0:
				for trip in newperms:
					vpathdicts[i-1][trip[0]] = None
	vpathdicts2 = [{} for i in range(len(th))]
	for i in range(len(th)):
		for key,valueset in vpathdicts[i].items():
			for value in valueset:
				key2 = value[0]
				if key2 not in vpathdicts2[i]:
					vpathdicts2[i][key2]=set()
				vpathdicts2[i][key2].add((key,value[1],value[2]))
	#print(vpathdicts2)
	return vpathdicts2


def theta(perm):
	cd = code(perm)
	for i in range(len(cd)-1,0,-1):
		for j in range(i-1,-1,-1):
			if cd[j] < cd[i]:
				cd[i]+=1
	cd.sort(reverse=True)
	return cd

def add_perm_dict(d1,d2):
	for k,v in d2.items():
		d1[k] = d1.get(k,0)+v
	return d1

def elem_sym_poly(p,k,varl1,varl2,xstart=0,ystart=0):
	if p>k:
		return 0
	if p == 0:
		return 1
	if p == 1:
		res = varl1[xstart] - varl2[ystart]
		for i in range(1,k):
			res += varl1[xstart+i] - varl2[ystart+i]
		return res
	if p == k:
		res = (varl1[xstart] - varl2[ystart]) * (varl1[xstart+1] - varl2[ystart])
		for i in range(2,k):
			res *= (varl1[i+xstart] - varl2[ystart])
		return res
	mid = k//2	
	xsm = xstart + mid
	ysm = ystart + mid	
	kmm = k - mid
	res = elem_sym_poly(p,mid,varl1,varl2,xstart,ystart)+elem_sym_poly(p,kmm,varl1,varl2,xsm,ysm)
	for p2 in range(max(1,p-kmm),min(p,mid+1)):
		res += elem_sym_poly(p2,mid,varl1,varl2,xstart,ystart)*elem_sym_poly(p-p2,kmm,varl1,varl2,xsm,ysm-p2)
	return res

def elem_sym_func(k,i,u1,u2,v1,v2,udiff,vdiff,varl1,varl2):
	newk = k - udiff
	p = k-udiff-vdiff
	if newk<0 or p<0:
		return 0
	if p == 0:
		return 1
	yvars = [varl1[j+1] for j in range(len(u2),k)] + \
			[varl1[u2[j]] for j in range(len(u1),len(u2)) if u2[j] == j+1] + \
			[varl1[u2[j]] for j in range(len(u1)) if u1[j] == u2[j]]
	v3 = list(v2)+[j+1 for j in range(len(v2),i)]	
	zvars = [varl2[v3[i-1]]] + \
			[varl2[v3[j]] for j in range(len(v1),len(v3)) if v3[j]!=j+1 and j!=i-1] + \
			[varl2[v3[j]] for j in range(len(v1)) if v1[j]!=v3[j] and j!=i-1]
	return elem_sym_poly(p,newk,yvars,zvars)
