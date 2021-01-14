from random import random,randrange,sample,shuffle
from numpy import arange,array,cos,exp,inf,isnan,linspace,log,ma,matrix,mean,nan_to_num,percentile,pi,sin,sqrt,where,zeros
from copy import deepcopy
from scipy.spatial.distance import euclidean
from igraph import Graph
from numpy.random import lognormal,normal
from difflib import SequenceMatcher
import csv
import string
import rasterio
import os
import sys
import itertools
from time import sleep

def wgs84(lat_lon):
	return [90-lat_lon[0],lat_lon[1]-180]



def rand_tre(tre):
	r = random()
	if r<tre:
		return tre
	else:
		return r


##########################################################################################
###STRING MATCHING FUNCTIONS TO COMPARE PHENOTYPES
##########################################################################################
'''
Function to evaluate a consumer's ability to use a resource, based on the
resource and consumer phenotypes (i.e. strings including a random sample of letters of size varying
between 1 and 10), the trait adjacency matrix, and the minimum/maximum compatibility values
estimated empirically in the previous step. As explained above, compatibility is obtained by
summing up all the ij entries in the  adjacency matrix for each i trait of the consumer and each j
trait of the potential resource.
'''

def compute_tc(resource,consumer,ft,ft_mat, min_tc, max_tc):
	sc = 0
	for c_trait in consumer:
		for r_trait in resource:
			sc += ft_mat[ft.index(c_trait)][ft.index(r_trait)]
	val = 1-(max_tc-sc)/(max_tc-min_tc)
	if val>1:
		val = 1
	elif val<0:
		val = 0
	return val


##########################################################################################
### GENERATE FOOD WEBS FROM SPECIES POOL
##########################################################################################
#####threshold affects connectance; 0.5 or 0.75?


def get_all_pair(sp12,ft,ft_mat,min_tc,max_tc):
	bs_dict = dict([['A_A', [1.0, 9.0]],	['A_B', [16.0, 122.0]],['A_M', [0, 4553.0]],['A_R', [213.0, 726.0]],['B_A', [3.0, 11.0]],['B_B', [11.0, 14.0]],['B_M', [132.0, 1850.0]],['B_R', [178.0, 1171.0]],['M_A', [22.0, 94.0]],['M_B', [15.0, 20.0]],['M_M', [109.0, 198.0]],['M_R', [22.0, 94.0]],['R_A', [0.0, 2.0]],['R_B', [29.0, 71.0]],['R_M', [44.0, 103.0]],['R_R', [179.0, 1082.0]]])
	sp1,sp2 = sp12
	if sp1[2][0]>sp2[2][0] and sp1[2][1]<=sp2[2][0]+1:
		tc = sp2[6]+'_'+sp1[6]
		ll,ul = bs_dict[tc]
		if ll<=(sp1[3]/sp2[3])<=ul:
			fc = compute_tc(sp1[-2],sp2[-2],ft,ft_mat,min_tc,max_tc)
			if fc>=0.55:
				fc = 1-(1-fc)/(1-0.55)
				return ([sp2[0],sp1[0],fc])
	elif	 sp2[2][0]>sp1[2][0] and sp2[2][1]<=sp1[2][0]+1:
		tc = sp1[6]+'_'+sp2[6]
		ll,ul = bs_dict[tc]
		if ll<=(sp2[3]/sp1[3])<=ul:
			fc = compute_tc(sp2[-2],sp1[-2],ft,ft_mat,min_tc,max_tc)
			if fc>=0.55:
				fc = 1-(1-fc)/(1-0.55)
				return ([sp1[0],sp2[0],fc])
	else:
		pass


def net_from_pool(pool,glob_net,prune='no'): #removed randomness
	g = glob_net.subgraph(set([i[0] for i in pool]))
	bas = set([g.vs['name'].index(i[0]) for i in pool if round(i[2][1])==1.0])
	if bas == set([]):
		g = Graph(directed = True)
		g.es['weight'] = []
		g.vs['name'] = []
		return g
	ddd = g.shortest_paths_dijkstra(target = bas, mode="IN")
	#check the case where ddd = []
	to_del = [i for i in range(len(ddd)) if min(ddd[i]) == inf]
	if prune == 'yes':
		to_del += where(array(g.degree())==0)[0].tolist()
	g.delete_vertices(to_del)
	return g


def write_pool(pools,space,file):
	out = open(file,'w')
	for i in range(len(space)):
		out.write(','.join(map(str,space[i]+[j[0] for j in pools[i]]))+'\n')
	out.close()


def write_list(lll,space,file):
	out = open(file,'w')
	for i in range(len(space)):
		out.write(','.join(map(str,space[i]+[j for j in lll[i]]))+'\n')
	out.close()


def write_map_eq_file(pools,space,file):
	out = open(file,'w')
	for i in range(len(space)):
		lat = 90-space[i][0]
		lon = space[i][1]-180
		for sp in pools[i]:
			out.write(','.join(map(str,[sp[0],lat,lon]))+'\n')
	out.close()

###################################################################################################
### SIMULATE CO-EXTINCTION CASCADES
###################################################################################################

def ext_casc(g,to_del,spp_dict,tre_cons,tre_res):
	g.delete_vertices(to_del)
	gw = array(list(g.get_adjacency(attribute='weight')))
	names = g.vs['name']
	bas = [i for i in g.vs['name'] if round(spp_dict[i][2][1])==1.0] ##connected to plants, including omnivorous
	bas_id = [names.index(i) for i in g.vs['name'] if round(spp_dict[i][2][0])==1.0] ##obligate herbivores
	ddd = g.shortest_paths_dijkstra(target = bas, mode="IN")
	done = []#names.index(i) for i in bas]
	if bas == []:
		g = Graph(directed = True)
		g.es['weight'] = []
		g.vs['name'] = []
		return g
	to_del_g = [names.index(g.vs['name'][i]) for i in range(len(ddd)) if min(ddd[i]) == inf]
	wmat = ma.masked_invalid(gw*(gw.T/gw.sum(1)).T)
	to_del_m_cons = [i for i in where(wmat.sum(0)<tre_cons)[0] if i not in done+bas_id]
	to_del_m_res = [i for i in where(wmat.sum(1)>tre_res)[0] if i not in done]
	to_del_m = list(set(to_del_g+to_del_m_cons+to_del_m_res))
	while to_del_m!=[]:
		g.delete_vertices([names[i] for i in to_del_m])
		gw[:,to_del_m] = 0.0
		gw[to_del_m,:] = 0.0
		done+=to_del_m
		bas = [i for i in g.vs['name'] if round(spp_dict[i][2][1])==1.0]
		if bas == []:
			g = Graph(directed = True)
			g.es['weight'] = []
			g.vs['name'] = []
			return g
		ddd = g.shortest_paths_dijkstra(target = bas, mode="IN")
		to_del_g = [names.index(g.vs['name'][i]) for i in range(len(ddd)) if min(ddd[i]) == inf]
		wmat = ma.masked_invalid(gw*(gw.T/gw.sum(1)).T)
		to_del_m_cons = [i for i in where(wmat.sum(0)<tre_cons)[0] if i not in done+bas_id]
		to_del_m_res = [i for i in where(wmat.sum(1)>tre_res)[0] if i not in done]
		to_del_m = list(set(to_del_g+to_del_m_cons+to_del_m_res))
	return g



#####get neighboring nodes
def find_nei_p(y,x,d):
	angle = random()*2*pi
	xOff = cos(angle)*d
	yOff = sin(angle)*d
	new_x = int(round(x + xOff))
	new_y = int(round(y + yOff))
	if new_x>359:
		new_x-=360
	if new_x<0:
		new_x+=360
	if new_y>179:
		new_y = 179
	if new_y<0:
		new_y = 0
	return (new_y,new_x)



###################################################################################################
### SPECIES' BIDIMENSIONAL NICHE & ADAPTATION
###################################################################################################
def fsigmoid(x, a, b):
    return 1.0 / (1.0 + exp(-a*(x-b)))


def sigmoid(vmin,vmax,vmean,pmm=0.95):
	a = 0.0001
	b = (log(pmm/(1-pmm))+vmax*a)/a
	while fsigmoid(vmean,a,b)>0.0001:
		a+=0.01
		b = (log(pmm/(1-pmm))+vmax*a)/a
	c = 0.0001
	d = (log((1-pmm)/pmm)+vmin*c)/c
	while fsigmoid(vmean,c,d)<0.9998:
		c+=0.01
		d = (log((1-pmm)/pmm)+vmin*c)/c
	return a,b,c,d


def niche_sigma(sp,v,pot):
	tmean,tmax,tmin,ta,tb,tc,td,prmean,prmax,prmin,pra,prb,prc,prd = sp[1]
	if pot == 't':
		t = v
		if t<=tmean:
			pt = 1-fsigmoid(t,tc,td)
		else:
			pt = fsigmoid(t,ta,tb)
		return pt
	else:
		pr = v
		if pr<=prmean:
			ppr = 1-fsigmoid(pr,prc,prd)
		else:
			ppr = fsigmoid(pr,pra,prb)
		return ppr



def adapt(sp,loc_clim,p_adapt_factor):
	pars = sp[1][:]
	t = sp[4]*p_adapt_factor
	x0,y0 = array([pars[0],pars[7]])
	x1,y1 = loc_clim
	xmean,ymean = ((1-t)*x0+t*x1),((1-t)*y0+t*y1)
	sx,sy = xmean-x0,ymean-y0
	xmin = pars[2]+sx
	ymin = pars[9]+sy
	xmax = pars[1]+sx
	ymax = pars[8]+sy
	xa,xb,xc,xd = sigmoid(xmin,xmax,xmean)
	ya,yb,yc,yd = sigmoid(ymin,ymax,ymean)
	return [xmean,xmax,xmin,xa,xb,xc,xd,ymean,ymax,ymin,ya,yb,yc,yd]
