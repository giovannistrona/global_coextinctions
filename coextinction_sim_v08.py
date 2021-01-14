from coextinction_functions_v08 import*

try:
	rcp = sys.argv[1]
except:
	rcp = sample(['45','6','85'],1)[0]


###create output folder
sleep(randrange(60))
res_fold_n = 0
while True:
	if not os.path.exists('./results/'+str(res_fold_n)+'_'+rcp):
		try:
			os.makedirs('./results/'+str(res_fold_n)+'_'+rcp)
			break
		except:
			res_fold_n+=1
			pass
	else:
		res_fold_n+=1


res_fold_n = str(res_fold_n)+'_'+rcp
###LOAD REFERENCE RASTERS
rst_fn = 'earth.tif'
rst = rasterio.open(rst_fn)
earth_mat = rst.read(1)
meta = rst.meta.copy()
meta.update(dtype='float64')
meta.update(count = 22) #layers to be saved
sss=earth_mat.shape
land_use = rasterio.open('./land_use_layers/'+rcp+'/2019.tif').read(1)
eco_csv = [i for i in csv.reader(open('ecoregions.csv','r',encoding = "ISO-8859-1"))]
eco_dict = dict([[float(i[0]),i[1:]] for i in eco_csv[1:]])
eco_dict[0.0] = ['na','na','na','na']
eco_mat = rasterio.open('biomes.tif').read(1)
##########################################################################################
###GENERATE GLOBAL SPECIES POOL
#### functional traits
ft = string.ascii_lowercase
ft_mat = zeros([26,26])
for i in range(26):
	for j in range(26):
		if i!=j:
			sign = sample([-1,1],1)[0]
			ft_mat[i][j] = sign*random()
		else:
			ft_mat[i][j] = 0


res = []
for rep in range(1000000):
	resource = sample(ft,randrange(1,11))
	consumer = sample(ft,randrange(1,11))
	sc = 0
	for c_trait in consumer:
		for r_trait in resource:
			sc += ft_mat[ft.index(c_trait)][ft.index(r_trait)]#
	res.append(sc)


min_tc, max_tc  = min(res),max(res)
niche_data = [i for i in csv.reader(open('species_niches_'+rcp+'.csv','r'))]
tl_vs_m = [i for i in csv.reader(open('tl_vs_m.csv','r'))]
mass =  [i for i in csv.reader(open('mass_db.csv','r'))]
for i in range(len(mass)):
	if float(mass[i][1]) == 0:
		mass[i][1] = 0.1


mass_dict = dict([[i[0],[float(i[1]),i[2]]] for i in mass])
niche_mass = []
for i in niche_data:
	try:
		sp_mass = mass_dict[i[1].replace('_',' ')]
		niche_mass.append(sp_mass+i[2:])
	except:
		pass


###adjust by proportion of richness x taxon (according to IUCN)
tax_rich = [['M',5513],['B',10425],['R',10038],['A',7302]]
tot_r = float(sum([i[1] for i in tax_rich]))
tax_p = dict([[i[0],i[1]/tot_r] for i in tax_rich])
tot_rich = sum([i[1] for i in tax_rich])
spp = []
for taxon in ['A','B','M','R']:
	niche_mass_t = [i for i in niche_mass if i[1]==taxon]
	for i in range(int(round(tot_rich*tax_p[taxon]))):
		tl_vs_m = sample(tl_vs_m,len(tl_vs_m))
		rsp = sample(niche_mass_t,1)[0]
		t = rsp[1]
		tl = ''
		tre = 0.01
		while tl == '':
			for j in tl_vs_m:
				if j[2] == t and (max(float(j[1]),rsp[0])-min(float(j[1]),rsp[0]))/max(float(j[1]),rsp[0])<tre:
					tl = [round(float(j[5]),1)-1,round(float(j[4]),1)-1] #max trophic level; mean is j[3]
					if round(tl[1]) == 2 and j[6] == 'yes': #connected to invertebrate?
						tl[1] -= 1
			tre+=0.01
		phenotype = ''.join(sample(ft,randrange(1,11)))
		adapt_val = random()
		spp.append([str(len(spp)),list(map(float,rsp[2:])),tl,rsp[0],adapt_val,phenotype,t])


comb = itertools.combinations(spp,2)
map_f = map(get_all_pair,comb, itertools.repeat(ft),itertools.repeat(ft_mat),itertools.repeat(min_tc),itertools.repeat(max_tc))
res = [val for val in map_f if val is not None]
out = open('./results/'+res_fold_n+'/glob_net.csv','w')
for i in res:
	out.write(','.join(map(str,i))+'\n')


out.close()
glob_net = Graph.TupleList(res,directed=True,weights=True)
spp = [i for i in spp if i[0] in glob_net.vs['name']]
spp_dict = dict([[i[0],i] for i in spp]) #create a dictionary linking species identifiers to species features

tas_2015,pr_2015 = [],[]
for mon in range(60):
	tas_2015.append(rasterio.open('./clim_layers/'+rcp+'/tas'+str(mon)+'.tif').read(1))
	pr_2015.append(rasterio.open('./clim_layers/'+rcp+'/pr'+str(mon)+'.tif').read(1)*10**6)


tas_0_min,pr_0_min =  [[[min([lay[mon][j][i] for mon in range(60)]) for i in range(360)] for j in range(180)] for lay in [tas_2015,pr_2015]]
tas_0_max,pr_0_max =  [[[max([lay[mon][j][i] for mon in range(60)]) for i in range(360)] for j in range(180)] for lay in [tas_2015,pr_2015]]
all_locs = [[i,j] for i in range(180) for j in range(360) if land_use[i][j]>0]
glob_div = 4500
pools = []	#empty list of species pools; each pool is the list of species found in a given locality
nets = []	#empty list of food webs
space = []	#empty list of localities' position and climate (lat,lon,minT,maxT)
for i in all_locs:
	att = 0
	lat,lon = i
	luf = land_use[lat][lon]#**0.5
	pool_size = int(round(glob_div*luf))
	if pool_size<100:
		pool_size = 100
	while att<10:
		if pool_size>len(spp):
			pool_size = len(spp)
		min_t,max_t,min_p,max_p = tas_0_min[lat][lon],tas_0_max[lat][lon],pr_0_min[lat][lon],pr_0_max[lat][lon]
		pool=[]		#generate empty pool
		for sp in sample(spp,pool_size):	#drop the selected number of random species in the new locality, one after another
			if 0.05>=max([niche_sigma(sp,v,'t') for v in [min_t,max_t]]+[niche_sigma(sp,v,'pr') for v in [min_p,max_p]]):
				pool.append(sp)
		net = net_from_pool(pool,glob_net,prune='yes')
		if len(net.es)>0:
			g = net.copy()
			gw = array(list(g.get_adjacency(attribute='weight')))
			www = ma.masked_invalid(gw*(gw.T/gw.sum(1)).T)#[not_bas]
			tre_cons = percentile([val for val in www.sum(0) if val>0],5)
			tre_res = percentile(www.sum(1),90)
			net = ext_casc(net,[],spp_dict,tre_cons,tre_res)
		if len(net.vs)>0: 	#if the attempt to build a food web succeeds:
			nets.append(net)	#add the web to the list of networks;
			in_net=set(net.vs['name'])	#check which species are in the network
			pool = [j for j in pool if j[0] in in_net]	#include in pool only species that have links in the network
			pools.append(pool)	#add local pool to pools list
			space.append([lat,lon])	#add the new locality to the locality list
			print (len(net.vs),len(net.es))
			att = 10
		else:
			pool_size+=100
			att+=1


tre_vals_res,tre_vals_cons = [],[]
for i in range(len(nets)):
	if len(nets[i].es)>0:
		g = nets[i].copy()
		gw = array(list(g.get_adjacency(attribute='weight')))
		www = ma.masked_invalid(gw*(gw.T/gw.sum(1)).T)
		tre_vals_cons.append(min([val for val in www.sum(0) if val>0]))
		tre_vals_res.append(www.sum(1).max())
	else:
		tre_vals_res.append('null')
		tre_vals_cons.append('null')


min_tre_cons = min([i for i in tre_vals_cons if i!='null'])
max_tre_res = max([i for i in tre_vals_res if i!='null'])
tre_vals_cons = [i if i!='null' else min_tre_cons for i in tre_vals_cons]
tre_vals_res = [i if i!='null' else max_tre_res for i in tre_vals_res]
max_bas = [len([j for j in pool if j[2][1]==1]) for pool in pools]
get_loc_id = dict([[tuple(space[i]),i] for i in range(len(space))])
pre_col_steps = 100
p_adapt_factor = 0.01
adapt_freq = 0.001
col_count = 0
############################
#pre_col stage
############################
pools_ = deepcopy(pools)
for step in range(pre_col_steps): #set number of step
		for substep in range(len(pools)):#randrange(len(pools),len(pools)*10)): ###i is the source pool, j the recipient
			i = randrange(len(pools))
			if len(pools[i])>0:
				lat,lon = space[i]
				d = 1+lognormal(0, 1)
				nei = list(find_nei_p(lat,lon,d))
				nei_att = 0
				while nei not in space and nei_att<10:
					d = 1+lognormal(0, 0.5)
					nei = find_nei_p(lat,lon,d)
					nei_att+=1
				if nei in space:
					j = get_loc_id[tuple(nei)]
					lat,lon = space[j]
					min_t,max_t,min_p,max_p = tas_0_min[lat][lon],tas_0_max[lat][lon],pr_0_min[lat][lon],pr_0_max[lat][lon]
					#invasability of net j is not considered in this preliminary phase.
					col_sp=sample(pools[i],1)[0][:]	#sample a potential colonizer from locality i; note that a copy of the colonizer is created (by adding '[:]' at the end of the line); this is important to avoid that time of permanence is reset for the colonizer in both the source and the target locality (it needs to be reset only in the latter)
					if col_sp not in pools_[j]:
						if 0.05>= max([niche_sigma(col_sp,v,'t') for v in [min_t,max_t]]+[niche_sigma(col_sp,v,'pr') for v in [min_p,max_p]]):
							pools_[j].append(col_sp)
							col_count+=1#print ('colonization!')
		if step>0 and step%10 == 0:
			for i in range(len(pools_)):
				bas = [j for j in pools_[i] if j[2][1]==1]
				if len(bas)>max_bas[i]:
					lat,lon = space[i]
					min_t,max_t,min_p,max_p = tas_0_min[lat][lon],tas_0_max[lat][lon],pr_0_min[lat][lon],pr_0_max[lat][lon]
					to_del_bas = sorted([[max([niche_sigma(bsp,v,'t') for v in [min_t,max_t]]+[niche_sigma(bsp,v,'pr') for v in [min_p,max_p]]),bsp] for bsp in bas])[max_bas[i]:]
					pools_[i] = [j for j in pools_[i] if j not in [k[1] for k in to_del_bas]]
			nets = [net for net in map(net_from_pool,pools_,itertools.repeat(glob_net),itertools.repeat('yes'))]
			nets_ = []
			for net in range(len(nets)):
				net_ = ext_casc(nets[net],[],spp_dict,tre_vals_cons[net],tre_vals_res[net])
				nets_.append(net_)
			nets = nets_
			pools = [[sp for sp in pools_[j] if sp[0] in nets[j].vs['name']] for j in range(len(pools_))]
			pools_ = deepcopy(pools)
		print (step,col_count)


for i in range(len(pools_)):
	bas = [j for j in pools_[i] if j[2][1]==1]
	if len(bas)>max_bas[i]:
		lat,lon = space[i]
		min_t,max_t,min_p,max_p = tas_0_min[lat][lon],tas_0_max[lat][lon],pr_0_min[lat][lon],pr_0_max[lat][lon]
		to_del_bas = sorted([[max([niche_sigma(bsp,v,'t') for v in [min_t,max_t]]+[niche_sigma(col_sp,v,'pr') for v in [min_p,max_p]]),bsp] for bsp in bas])[max_bas[i]:]
		pools_[i] = [j for j in pools_[i] if j not in [k[1] for k in to_del_bas]]


nets = [net for net in map(net_from_pool,pools_,itertools.repeat(glob_net),itertools.repeat('yes'))]
nets_ = []
for net in range(len(nets)):
	net_ = ext_casc(nets[net],[],spp_dict,tre_vals_cons[net],tre_vals_res[net])
	nets_.append(net_)


nets = nets_
pools = [[sp for sp in pools_[j] if sp[0] in nets[j].vs['name']] for j in range(len(pools_))]
pools_,nets_,space_,max_bas_ = [],[],[],[]
for i in range(len(space)):
	if len(pools[i])>0:
		pools_.append(pools[i])
		nets_.append(nets[i])
		space_.append(space[i])
		max_bas_.append(max_bas[i])


pools = pools_
nets = nets_
space = space_
max_bas = max_bas_
get_loc_id = dict([[tuple(space[i]),i] for i in range(len(space))])
header = 'year,mean_div,mean_div_control,comp_ext,clim_ext,lu_ext,sec_ext,adapted\n'
out = open('./results/'+res_fold_n+'/summary.csv','w')
out.write(header)
out.close()
############################
###simulation phase
############################
tre_vals_res,tre_vals_cons = [],[]
for i in range(len(nets)):
	if len(nets[i].es)>0:
		g = nets[i].copy()
		gw = array(list(g.get_adjacency(attribute='weight')))
		www = ma.masked_invalid(gw*(gw.T/gw.sum(1)).T)
		tre_vals_cons.append(min([val for val in www.sum(0) if val>0]))
		tre_vals_res.append(www.sum(1).max())
	else:
		tre_vals_res.append('null')
		tre_vals_cons.append('null')


min_tre_cons = min([i for i in tre_vals_cons if i!='null'])
max_tre_res = max([i for i in tre_vals_res if i!='null'])
tre_vals_cons = [i if i!='null' else min_tre_cons for i in tre_vals_cons]
tre_vals_res = [i if i!='null' else max_tre_res for i in tre_vals_res]
###prepare results files
eco_reg = [eco_dict[eco_mat[i[0]][i[1]]][1:-1] for i in space]
var_names = 'year,lat,lon,eco_reg,biome,div1,div2,div_tl2_1,div_tl2_2,div_tl3_1,div_tl3_2,div_tl4_1,div_tl4_2,tl1_mean,tl2_mean,tl1_max,tl2_max,bs1_mean,bs2_mean,bs1_max,bs2_max,node_n,edge_n,conn,diam,giant,clust_n\n'
out_rel,out_raw = open('./results/'+res_fold_n+'/res_rel.csv','w'),open('./results/'+res_fold_n+'/res_raw.csv','w')
out_rel.write(var_names)
for loc in range(len(space)):
	out_rel.write(','.join(map(str,[2015]+space[loc]+eco_reg[loc]+[1.0 for i in range(len(var_names)-5)]))+'\n')


out_rel.close()
out_raw.write(var_names)
out_raw.close()
land_use_0 = rasterio.open('./land_use_layers/'+rcp+'/2019.tif').read(1)
land_use = rasterio.open('./land_use_layers/'+rcp+'/2019.tif').read(1)
pools_1=deepcopy(pools)	#create a copy of the pools list for the
pools_2=deepcopy(pools)
nets_1=deepcopy(nets)	#igraph weighted networks
time_step = 0
for year in range(2015,2101): #set number of step
	coloniz,comp_ext,clim_ext,sec_ext,lu_ext,adapted = [[] for i in range(len(pools))],[[] for i in range(len(pools))],[[] for i in range(len(pools))],[[] for i in range(len(pools))],[[] for i in range(len(pools))],[[] for i in range(len(pools))]
	tas_y,pr_y = [],[] #save yearly climate for adaptation
	if year>2019:
		land_use = rasterio.open('./land_use_layers/'+rcp+'/'+str(year)+'.tif').read(1)
		land_use_0 = rasterio.open('./land_use_layers/'+rcp+'/'+str(year-1)+'.tif').read(1)
	for month in range(12):
		tas = rasterio.open('./clim_layers/'+rcp+'/tas'+str(time_step)+'.tif').read(1)
		pr = rasterio.open('./clim_layers/'+rcp+'/pr'+str(time_step)+'.tif').read(1)*10**6
		tas_y.append(tas)
		pr_y.append(pr)
		time_step+=1
	pools_2_ = deepcopy(pools_2)
	for i in range(len(pools)): ###i is the source pool, j the recipient
		lat,lon = space[i]
		if land_use_0[lat][lon]>0:
			d = 1+lognormal(0, 1)
			nei = list(find_nei_p(lat,lon,d))
			nei_att = 0
			while nei not in space and nei_att<10:
				d = 1+lognormal(0, 0.5)
				nei = find_nei_p(lat,lon,d)
				nei_att+=1
			if nei in space:
				j = get_loc_id[tuple(nei)]
				lat_j,lon_j = space[j]
				loc_tas,loc_pr = tas[lat_j][lon_j],pr[lat_j][lon_j]
				if len(pools_2[i])>0:
					col_sp = sample(pools_2[i],1)[0][:]	#sample a potential colonizer from locality i; note that a copy of the colonizer is created (by adding '[:]' at the end of the line); this is important to avoid that time of permanence is reset for the colonizer in both the source and the target locality (it needs to be reset only in the latter)
					if col_sp not in pools_2_[j]:
						pools_2_[j].append(col_sp)
						coloniz[j].append(col_sp[0])
				if len(pools_1[i])>0:
					col_sp=sample(pools_1[i],1)[0][:] #select at random a potential colonizer; note that a copy of the colonizer is created (by adding '[:]' at the end of the line); this is important to avoid that time of permanence is reset for the colonizer in both the source and the target locality (it needs to be reset only in the latter)
					if col_sp not in pools_1[j]:
						pools_1[j].append(col_sp)
			loc_clim = array([mean([tm[lat][lon] for tm in tas_y]),mean([prm[lat][lon] for prm in pr_y])])
			loc_clim_min = array([min([tm[lat][lon] for tm in tas_y]),min([prm[lat][lon] for prm in pr_y])])
			loc_clim_max = array([max([tm[lat][lon] for tm in tas_y]),max([prm[lat][lon] for prm in pr_y])])
			for sp in range(len(pools_1[i])):
				if random()<adapt_freq:
					pools_1[i][sp][1] = adapt(pools_1[i][sp],loc_clim,p_adapt_factor)
			for sp in range(len(pools_2_[i])):
				if random()<adapt_freq:
					pools_2_[i][sp][1] = adapt(pools_2_[i][sp],loc_clim,p_adapt_factor)
					adapted[i].append(pools_2_[i][sp][0])
			bas_2 = [j for j in pools_2_[i] if j[2][1]==1]
			if len(pools_1[i])>len(pools[i]):
				to_del_pool = sorted([[max([niche_sigma(psp,loc_clim_min[0],'t'),niche_sigma(psp,loc_clim_max[0],'t'),niche_sigma(psp,loc_clim_min[1],'pr'),niche_sigma(psp,loc_clim_max[1],'pr')]),psp] for psp in pools_1[i]])[len(pools[i]):]
				pools_1[i] = [j for j in pools_1[i] if j not in [k[1] for k in to_del_pool]]
			if len(bas_2)>max_bas[i]:
				to_del_bas = sorted([[max([niche_sigma(bsp,loc_clim_min[0],'t'),niche_sigma(bsp,loc_clim_max[0],'t'),niche_sigma(bsp,loc_clim_min[1],'pr'),niche_sigma(bsp,loc_clim_max[1],'pr')]),bsp] for bsp in bas_2])[max_bas[i]:]
				pools_2_[i] = [j for j in pools_2_[i] if j not in [k[1] for k in to_del_bas]]
				comp_ext[i] = set([k[1][0] for k in to_del_bas])-set(coloniz[i])
	nets_1 = [net for net in map(net_from_pool,pools_2_,itertools.repeat(glob_net),itertools.repeat('no'))]
	pools_2 = [[sp for sp in pools_2_[j] if sp[0] in nets_1[j].vs['name']] for j in range(len(pools))]
	for i in range(len(pools)): #iterate over localities/communities
		lat,lon = space[i]
		if land_use_0[lat][lon]>0:
			ext_1,ext_2=[],[]	#identify primary extinctions due to T change
			if land_use_0[lat][lon]>land_use[lat][lon]:
				luf = (1-(land_use_0[lat][lon]-land_use[lat][lon])/land_use_0[lat][lon])#**0.5
				ext_1 = [j[0] for j in sample(pools_1[i],len(pools_1[i])-int(round(luf*len(pools_1[i]))))]
				ext_2 = [j[0] for j in sample(pools_2[i],len(pools_2[i])-int(round(luf*len(pools_2[i]))))]
				lu_ext[i] = set(ext_2[:])-set(coloniz[i])
			for month in range(12):
				for j in pools_1[i]:	# evaluate the compatibility between all species in the pool with the new climate in the tolerance scenario
					#if random()<max(niche_sigma(j,tas_y[month][lat][lon],'t'),niche_sigma(j,pr_y[month][lat][lon],'pr'))**1.5:
					if rand_tre(0.05)<max(niche_sigma(j,tas_y[month][lat][lon],'t'),niche_sigma(j,pr_y[month][lat][lon],'pr')):
						ext_1.append(j[0])
				for j in pools_2[i]:	# do the same for the co-extinction scenario; evaluation is made separately, because species pools will clearly become different one from another throughout the simulation
					if rand_tre(0.05)<max(niche_sigma(j,tas_y[month][lat][lon],'t'),niche_sigma(j,pr_y[month][lat][lon],'pr')):
						ext_2.append(j[0])
			clim_ext[i] = set(ext_2)-(set(lu_ext[i])|set(coloniz[i]))
			if len(nets_1[i].vs)>0:	#check if local food web has extant links
				names=nets_1[i].vs['name']
				new_g = ext_casc(nets_1[i],ext_2,spp_dict,tre_vals_cons[i],tre_vals_res[i])
				nets_1[i]=new_g #replace the starting network with the one obtained after the co-extinction cascade
			sec_ext[i] = set([j[0] for j in pools_2[i] if j[0] not in nets_1[i].vs['name']])-set(ext_2)
			pools_1[i]=[j for j in pools_1[i] if j[0] not in ext_1] #eliminate species that went extinct (primary extinctions)
			pools_2[i] = [j for j in pools_2[i] if j[0] in nets_1[i].vs['name']] #reduce species pool to species still in the network after the co-extinction cascade
	ext_fract = [sum([len(ds[row]) for row in range(len(pools))]) for ds in [comp_ext,clim_ext,lu_ext,sec_ext,adapted]]
	res = [year,mean([len(i) for i in pools_2]),mean([len(i) for i in pools_1])]+ext_fract
	print (res)
	out = open('./results/'+res_fold_n+'/summary.csv','a')
	row_l = out.write(','.join(list(map(str,res)))+'\n')
	out.close()
	div1 = array([len(set([i[0] for i in pool])) for pool in pools_2],dtype=float)
	div2 = array([len(set([i[0] for i in pool])) for pool in pools_1],dtype=float)
	div_tl2_1 = array([len([i for i in pool if i[2][0]>1]) for pool in pools_2],dtype=float)
	div_tl2_2 = array([len([i for i in pool if i[2][0]>1]) for pool in pools_1],dtype=float)
	div_tl3_1 = array([len([i for i in pool if i[2][0]>2]) for pool in pools_2],dtype=float)
	div_tl3_2 = array([len([i for i in pool if i[2][0]>2]) for pool in pools_1],dtype=float)
	div_tl4_1 = array([len([i for i in pool if i[2][0]>3]) for pool in pools_2],dtype=float)
	div_tl4_2 = array([len([i for i in pool if i[2][0]>3]) for pool in pools_1],dtype=float)
	tl1_mean = array([mean([i[2][0] for i in pool]) if pool!=[] else 0 for pool in pools_2],dtype=float)
	tl2_mean = array([mean([i[2][0] for i in pool]) if pool!=[] else 0 for pool in pools_1],dtype=float)
	tl1_max = array([max([i[2][0] for i in pool]) if pool!=[] else 0 for pool in pools_2],dtype=float)
	tl2_max = array([max([i[2][0] for i in pool]) if pool!=[] else 0 for pool in pools_1],dtype=float)
	bs1_mean = array([mean([i[3] for i in pool]) if pool!=[] else 0 for pool in pools_2],dtype=float)
	bs2_mean = array([mean([i[3] for i in pool]) if pool!=[] else 0 for pool in pools_1],dtype=float)
	bs1_max = array([max([i[3] for i in pool]) if pool!=[] else 0 for pool in pools_2],dtype=float)
	bs2_max = array([max([i[3] for i in pool]) if pool!=[] else 0 for pool in pools_1],dtype=float)
	node_n = array([float(len(net.vs)) for net in nets_1])
	edge_n = array([float(len(net.es)) for net in nets_1])
	conn = array([edge_n[i]/(node_n[i])**2 if node_n[i]>0 else 0.0 for i in range(len(node_n))])
	diam = array([net.diameter(directed=True) for net in nets_1])
	clusts = [net.clusters(mode='weak') for net in nets_1]
	giant = array([max([len(i) for i in clusts[clust]])/node_n[clust] if node_n[clust]!=0  else 0.0 for clust in range(len(clusts))])
	clust_n = array([len(clusts[clust])/node_n[clust] if node_n[clust]!=0  else 0.0 for clust in range(len(clusts))])
	res_all = [div1,div2,div_tl2_1,div_tl2_2,div_tl3_1,div_tl3_2,div_tl4_1,div_tl4_2,tl1_mean,tl2_mean,tl1_max,tl2_max,bs1_mean,bs2_mean,bs1_max,bs2_max,node_n,edge_n,conn,diam,giant,clust_n]
	if year == 2015:
		res_0 = deepcopy(res_all)
	out_rel,out_raw = open('./results/'+res_fold_n+'/res_rel.csv','a'),open('./results/'+res_fold_n+'/res_raw.csv','a')
	for loc in range(len(space)):
		out_raw.write(','.join(map(str,[year]+space[loc]+eco_reg[loc]+[i[loc] for i in res_all]))+'\n')
		if year>2015:
			out_rel.write(','.join(map(str,[year]+space[loc]+eco_reg[loc]+[(res_all[i]/res_0[i])[loc] for i in range(len(res_all))]))+'\n')
	out_rel.close()
	out_raw.close()
	mmm = [zeros(sss) for lay in range(22)]
	lays = [div1,div2,div_tl2_1,div_tl2_2,div_tl3_1,div_tl3_2,div_tl4_1,div_tl4_2,tl1_mean,tl2_mean,tl1_max,tl2_max,bs1_mean,bs2_mean,bs1_max,bs2_max,node_n,edge_n,conn,diam,giant,clust_n]
	for i in range(len(pools)):
		lat,lon = space[i]
		for band in range(22):
			mmm[band][lat][lon] = lays[band][i]
	out=rasterio.open('./results/'+str(res_fold_n)+'/'+str(year)+'.tif', 'w', **meta)
	for band in range(22):
		out.write_band(band+1,mmm[band].astype(rasterio.float64))
	out.close()
