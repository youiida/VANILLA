from UDFManager import *
from sys import *
import os
import math
import gnuplot
from numpy import *
from collections import defaultdict
from operator import itemgetter
from math import fsum
import MakeGraphSheetLib

def less(a,b):
	if a<b:
		return a
	else:
		return b
def more(a,b):
	if a>b:
		return a
	else:
		return b
def fh(chi, x, n1, n2):
	if x>=1.0:
		print("change mesh or dx!",x)
		x = 1.0-1.0e-9
	elif x<=0.0:
		print("change mesh or dx!",x)
		x = 1.0e-9
	return chi*x*(1.0-x) + x*math.log(x)/n1 + (1.0-x)*math.log(1.0-x)/n2
def cp(chi, x, n1, n2):
	if x>=1.0:
		print("change mesh or dx!",x)
		x = 1.0-1.0e-9
	elif x<=0.0:
		print("change mesh or dx!",x)
		x = 1.0e-9
	return chi*(1.0-2.0*x) + (1.0+math.log(x))/n1 - (1.0+math.log(1.0-x))/n2
def cp1d(chi, x, n1, n2):
	return -2.0*chi + 1.0/(x*n1) + 1.0/ ((1.0-x)*n2)
def phi_sp(chi, n1, n2):
	val  = 2.0*n1*n2*chi
	coef = n1-n2-val
	root = math.sqrt(coef*coef-4.0*val*n2)
	return -0.5*(coef+root)/val, -0.5*(coef-root)/val
def critical(n1, n2):
	phic  = math.sqrt(n2)/(math.sqrt(n1)+math.sqrt(n2))
	chic  = 0.5*(math.sqrt(n1)+math.sqrt(n2))*(math.sqrt(n1)+math.sqrt(n2))/(n1*n2)
	return chic, phic
def intercept(x, y, grad):
	return y-grad*x
def icDiff(chi, fix, opt, n1, n2, th):
	error   = 1.0
	target  = cp(chi, fix, n1, n2)
	old_phi = opt
	counter = 0
	dx = th
	while error > th:
		new_phi = old_phi - (cp(chi, old_phi, n1, n2)-target) / cp1d(chi, old_phi, n1, n2)
		if new_phi >= 1.0:
			new_phi = 1.0-dx
			dx = dx*0.5
		if new_phi <= 0.0:
			new_phi = dx
			dx = dx*0.5
		error   = abs(target - cp(chi, new_phi, n1, n2))
		old_phi = new_phi
		counter = counter + 1
		if counter>10000:
			print("can not converge(icDiff)!")
			break
	
	opt_phi = old_phi
	opt_mu  = cp(chi, new_phi, n1, n2)
	opt_fe  = fh(chi, opt_phi, n1, n2)
	opt_ic  = intercept(opt_phi, opt_fe, opt_mu)
	fix_phi = fix
	fix_mu  = target
	fix_fe  = fh(chi, fix_phi, n1, n2)
	fix_ic  = intercept(fix_phi, fix_fe, fix_mu)
	return  fix_ic-opt_ic, opt_phi, opt_mu

def gradIcDiff(chi, stag, fix, n1, n2, torelance):
	minus = stag - 0.5*torelance
	plus  = stag + 0.5*torelance
	
	if minus < 0.0:
		minus = abs(minus)
	elif minus == 0.0:
		minus = 0.1*torelance
	if plus > 1.0:
		plus  = 2.0 - plus
	elif plus == 0.0:
		plus = 1.0-0.1*torelance
		
	temp_m = less(minus, plus)
	temp_p = more(minus, plus)
	minus, plus = temp_m, temp_p
	delta  = plus - minus
	m_diff, m_opt_phi, m_opt_cp = icDiff(chi, minus, fix, n1, n2, torelance)
	p_diff, p_opt_phi, p_opt_cp = icDiff(chi, plus, fix, n1, n2, torelance)
	
	if m_diff == p_diff:
		return 1.0e-10
	
	return (p_diff-m_diff)/delta

def phase_diagram(_udf_, index1, index2, dev, max_plus, type, csv_file):
	raw_n1 = _udf_.get("component_properties.components[].degree_of_polymerization",[index1])
	raw_n2 = _udf_.get("component_properties.components[].degree_of_polymerization",[index2])
	name1 = _udf_.get("component_properties.components[].name",[index1])
	name2 = _udf_.get("component_properties.components[].name",[index2])
	
	scale = 1.0/raw_n1
	n1    = 1.0
	n2    = raw_n2 * scale
	rsc   = raw_n1
	
	dataList = []#binodal, spinodal phi and chi
	chic, phic = critical(n1, n2)
	dataList.append([phic, chic/rsc])
	
	
	if max_plus<0.0:
		max = chic + 1.0
	else:
		max = chic + max_plus*rsc
	chi = chic

	binordal  = []
	spinordal = []

	spinordal.append([phic, chic/rsc])
	binordal.append([phic, chic/rsc])

	delta = (max-chic)/dev
	
	while chi < max:
		chi = chi + delta
		sp_l, sp_r = phi_sp(chi, n1, n2)
		if type=="spinodal":
			dataList.append([sp_l, chi/rsc])
			dataList.append([sp_r, chi/rsc])
		spinordal.append([sp_l, chi/rsc])
		spinordal.append([sp_r, chi/rsc])

	chi = chic

	if type=="binodal" or csv_file=="yes":
		while chi < max:
			flag = "true"
			chi = chi + delta
			dx  = 1.0/float(10000)
			th  = dx
		
			sp_l, sp_r = phi_sp(chi, n1, n2)
		
			mu = []
			phi = dx
			while phi < 1.0:
				mu.append([phi, cp(chi, phi, n1, n2)])
				phi = phi + dx
			lm = []
			prev = -1
			for i in range(len(mu)):
				if mu[i][1] > 0.0:
					now = 1
				else:
					now = -1
				if prev*now < 0.0 and now > 0.0:
					lm.append([mu[i][0], fh(chi, mu[i][0], n1, n2)])
				prev = now
			
			min = 0.0
			for i in lm:
				if fh(chi, i[0], n1, n2)<min:
					phi_min = i[0]
					min = fh(chi, i[0], n1, n2)
			
			if phi_min < sp_l:
				phi_b2 = sp_r
				sign = 1.0
			elif phi_min > sp_r:
				phi_b2 = sp_l
				sign = -1.0
			else:
				flag = "false"
				print("aho")
				break
		
			delta_fe = 1.0
			phi_b1   = phi_min
			counter  = 0
			while phi_b1 > 0.0:
				phi_var = phi_b2
				while phi_var > 0.0 and phi_var < 1.0 and delta_fe > 0.0:
					delta_fe = fh(chi, phi_var, n1, n2) - (cp(chi, phi_b1, n1, n2)*(phi_var-phi_b1) + fh(chi, phi_b1, n1, n2))
					phi_var  = phi_var + sign*dx
				if delta_fe < 0.0:
					break
				counter = counter + 1
				if counter > 10000:
					print("can not converge(delta_fe)")
					flag = "false"
					break
				phi_b1 = phi_b1 + sign*dx
				
			error   = 1.0
			old_phi = phi_var
			counter = 0
			while error > th:
				diff, opt_phi, opt_cp = icDiff(chi, old_phi, phi_b1, n1, n2, th)
				new_phi = old_phi - diff / gradIcDiff(chi, old_phi, phi_b1, n1, n2, th)
				if new_phi >= 1.0:
					new_phi = 1.0-dx
					dx = dx*0.5
				if new_phi <= 0.0:
					new_phi = dx
					dx = dx*0.5
				error = abs(cp(chi, new_phi, n1, n2) - opt_cp)
				old_phi = new_phi
				counter = counter + 1
				if counter > 10000:
					flag = "false"
					print("can not converge(opt_phi)")
					break
			
			if flag == "true":
				dataList.append([new_phi, chi/rsc])
				dataList.append([opt_phi, chi/rsc])
				binordal.append([new_phi, chi/rsc])
				binordal.append([opt_phi, chi/rsc])
	
	dataList.sort(key=lambda x:x[0])

	binordal.sort(key=lambda x:x[0])
	spinordal.sort(key=lambda x:x[0])

	if csv_file=="yes":
		file = open("phase_diagram_"+name1+"_"+name2+".csv", 'w')

	
	result = [[],[]]
	for i in range(len(dataList)):
		result[0].append(dataList[i][0])
		result[1].append(dataList[i][1])
	
	title_string = "phase diagram (" + name1 + " and " + name2 + ")"
	tag_list = [name1 + " volume fraction", type + ": " +name1+" "+str(raw_n1)+"  and  "+name2+" "+str(raw_n2)]
	
	gnuplot.plot(data=result, title=title_string, labels=tag_list)


	if csv_file=="yes":
		file.write("#binordal_phi,binordal_chi,spinordal_phi,spinordal_chi,n1="+str(raw_n1)+",n2="+str(raw_n2)+"\n")
		for i in range(len(binordal)):
			file.write(str(binordal[i][0]) + "," + str(binordal[i][1]) + "," 
					 + str(spinordal[i][0]) + "," + str(spinordal[i][1]) + "\n")
		file.close()

def shimazaki_histogram(data, start, end):
	_max = max(data)
	_min = min(data)

	results = []
	
	for N in range(start, end):
		width = float(_max - _min) / N

		hist = defaultdict(int)
		for x in data:
			i = int((x - _min) / width)
			if i >= N:       # Mimicking the behavior of matlab.histc(), and
				i = N - 1    # matplotlib.hist() and numpy.histogram().
			y = _min + width * i
			hist[y] += 1

		# Compute the mean and var.
		k = fsum(hist[x] for x in hist) / N
		v = fsum(hist[x]**2 for x in hist) / N - k**2

		C = (2 * k - v) / (width**2)

		results += [(hist, C, N, width)]

	optimal = min(results, key=itemgetter(1))

	if 0: # if true, print bin-widths and C-values, the cost function.
		for (hist, C, N, width) in results:
			print(width, C)

	return optimal


def free_energy(_udf_, time):
	title_string	=	"free energy"
	if time=="record":
		locationList	=	[
			'output.steps',
			'output.free_energy'
		]
		tag_list	= [
			'record',
			'free energy'
		]
	elif time=="time":
		locationList	=	[
			'output.time',
			'output.free_energy'
		]
		tag_list	= [
			'time',
			'free energy'
		]
	
	dataList=MakeGraphSheetLib.makeGraphSheet(_udf_,locationList,tag_list)
	gnuplot.plot(data=dataList,title=title_string,labels=tag_list)

	#console window, ascii file
	name = _udf_.udfFilename().split(".")[0] + "_free_energy.dat"
	file = open(name,'w')
	file.write("#"+time+"\tenergy\n")
	
	print(time+"\t",_udf_.udfFilename().split(".")[0])

	if time=="record":
		for rec in range(0,_udf_.totalRecord()):
			_udf_.jump(rec)
			step=_udf_.get("output.steps")
			fe = _udf_.get("output.free_energy")
			print(str(step)+"\t"+str(fe))
			file.write(str(step)+"\t"+str(fe)+"\n")
	elif time=="time":
		for rec in range(0,_udf_.totalRecord()):
			_udf_.jump(rec)
			time=_udf_.get("output.time")
			fe = _udf_.get("output.free_energy")
			print(str(time)+"\t"+str(fe))
			file.write(str(time)+"\t"+str(fe)+"\n")

def volume_fraction(_udf_):
	title_string = "volume fractions"
	tag_list = ["record"]
	dataList = [[]]
	total = _udf_.totalRecord()
	comps = _udf_.size("output.components[]")
	print(total,comps)
	
	for c in range(comps):
		name = _udf_.get("output.components[].name",[c])
		tag_list.append(name)
		dataList.append([])
		
	for rec in range(total):
		_udf_.jump(rec)
		dataList[0].append(rec)
		for c in range(comps):
			vol = _udf_.get("output.components[].volume_fraction",[c])
			dataList[c+1].append(vol)
			
	gnuplot.plot(data=dataList,title=title_string,labels=tag_list)

def concentration_distribution(_udf_, interval, index):
	all_data = []
	tag_list = ["volume fraction"]
	
	records = int(_udf_.totalRecord()/interval) + 1
	
	for i in range(records):
		rec = i*interval
		_udf_.jump(rec)
		name = _udf_.get("output.components[].name",[index])
		tag_list.append(name+":record "+str(rec))
		
		size = _udf_.size("output.components[].density[]",[index])
		for ix in range(size):
			vol = _udf_.get("output.components[].density[]",[index,ix])
			all_data.append(vol)
	
	
	optimal = shimazaki_histogram(all_data, 1, 100)
	(hist, C, N, width) = optimal
	key = []	
	for x in hist:
		key.append(x)
	key.sort()
	
	title_string="concentration distribution"

	dataList = []
	dataList.append(key)
	
	for i in range(records):
		rec = i*interval
		_udf_.jump(rec)
		data = []
		count= []
		for x in range(len(key)):
			count.append(0)
			
		size = _udf_.size("output.components[].density[]",[index])
		for ix in range(size):
			vol = _udf_.get("output.components[].density[]",[index,ix])
			data.append(vol)
			
		for j in data:
			for k in range(len(dataList[0]))[::-1]:
				if j >= dataList[0][k]:
					count[k] = count[k] + 1
					break
		for y in range(len(count)):
			count[y] = float(count[y])/float(size)
		dataList.append(count)


	gnuplot.plot(data=dataList,title=title_string,labels=tag_list)

def rich_phase_rate(_udf_, index, threshold):
	title_string = "rich phase rate"
	name = _udf_.get("output.components[].name",[index])
	tag_list = ["record", name + "-rich phase rate"]
	data_list = [[],[]]
	
	if threshold < 0:
		threshold = _udf_.get("dynamics.visco_elastic.critical_phi")
	size = _udf_.size("output.components[].density[]",[index])
	for rec in range(_udf_.totalRecord()):
		_udf_.jump(rec)
		counter = 0
		data_list[0].append(rec)
		for ix in range(size):
			val = _udf_.get("output.components[].density[]", [index, ix])
			if val >= threshold:
				counter = counter + 1
		data_list[1].append(float(counter)/float(size))

	gnuplot.plot(data=data_list,title=title_string,labels=tag_list)

def px(i,nx):
	if (i+1)%nx!=0:
		return i+1
	else:
		return i -(nx-1)
def py(i,nx,ny):
	if i%(nx*ny) < nx*(ny-1):
		return i+nx
	else:
		return i-nx*(ny-1)
def pz(i,nx,ny,nz):
	if i<nx*ny*(nz-1):
		return i+nx*ny
	else:
		return i-nx*ny*(nz-1)

def surface_area(_udf_, index, time):
	title_string = "surface area"
	name = _udf_.get("output.components[].name",[index])
	dx = _udf_.get("simulation_conditions.mesh_conditions.mesh_size.dx")
	dy = _udf_.get("simulation_conditions.mesh_conditions.mesh_size.dy")
	dz = _udf_.get("simulation_conditions.mesh_conditions.mesh_size.dz")
	nx = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nx")
	ny = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Ny")
	nz = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nz")

	tag_list = []
	if time=="record":
		tag_list.append("record")
	else:
		tag_list.append("time")

	
	dataList = [[],[]]
	t=0
	for rec in range(_udf_.totalRecord()):
		_udf_.jump(rec)
		
		if time=="record":
			t = _udf_.get("output.steps")
		else:
			t = _udf_.get("output.time")
		dataList[0].append(t)
		
		th = _udf_.get("output.components[].volume_fraction",[index])
		
		
		data = []
		for i in range(_udf_.size("output.components[].density[]",[index])):
			val = _udf_.get("output.components[].density[]",[index,i])
			data.append(val)
		area = 0.0
		for i in range(len(data)):
			if (data[i]-th)*(data[px(i,nx)]-th)<0.0:
				area = area + dy*dz
			if (data[i]-th)*(data[py(i,nx,ny)]-th)<0.0:
				area = area * dz*dx
			if (data[i]-th)*(data[pz(i,nx,ny,nz)]-th)<0.0:
				area = area + dx*dy
		dataList[1].append(area)

	gnuplot.plot(data=dataList,title=title_string,labels=tag_list)

def xyz(index, nx,ny,nz):
	x = index%nx
	z = index/(nx*ny)
	y = (index-x-nx*ny*z)/nx
	return x,y,z
	
def scattering1d(_udf_, index):
	dx = _udf_.get("simulation_conditions.mesh_conditions.mesh_size.dx")
	dy = _udf_.get("simulation_conditions.mesh_conditions.mesh_size.dy")
	dz = _udf_.get("simulation_conditions.mesh_conditions.mesh_size.dz")
	nx = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nx")
	ny = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Ny")
	nz = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nz")
	name = _udf_.get("output.components[].name",[index])
	field = zeros((nx,ny,nz))
	volume = nx*ny*nz

	for i in range(_udf_.size("output.components[].density[]",[index])):
		val = _udf_.get("output.components[].density[]",[index,i])
		x,y,z = xyz(i,nx,ny,nz)
		field[x][y][z] = val
	cfield = fft.fftn(field)

	rdata = []
	for i in range(nx/2):
		for j in range(ny/2):
			for k in range(nz/2):
				r = math.sqrt( (i*dx)*(i*dx) + (j*dy)*(j*dy) + (k*dz)*(k*dz) )
				rdata.append(r)
	max=-1
	if max<nx:
		max = nx
	if max<ny:
		max = ny
	if max<nz:
		max = nz
	
	optimal = shimazaki_histogram(rdata, 1, max)
	(hist, C, N, width) = optimal
	key = []	
	for x in hist:
		key.append(x)
	key.sort()
	xvalues = []
	counter = []
	yvalues = []
	
	for i in range(len(key)-1):
		xvalues.append((key[i]+key[i+1])*0.5)
		counter.append(0)
		yvalues.append(0.0)
	
	for i in range(-nx/2, nx/2):
		for j in range(-ny/2, ny/2):
			for k in range(-nz/2, nz/2):
				r = math.sqrt( (i*dx)*(i*dx) + (j*dy)*(j*dy) + (k*dz)*(k*dz) )
				val = abs(cfield[i+nx/2][j+ny/2][k+nz/2])/volume
				for s in range(len(key)-1):
					if r>key[s] and r<key[s+1]:
						counter[s] = counter[s] + 1
						yvalues[s] = yvalues[s] + val
	
	for s in range(len(key)-1):
		yvalues[s] = yvalues[s]/float(counter[s])
	dataList = [xvalues,yvalues]
	tag_list = ["q-range","intensity"]
	title_string = "scattering(1d, " + name + ")"
	gnuplot.plot(data=dataList,title=title_string,labels=tag_list)
	

def get_time_obj(_udf_):
	model = _udf_.get("dynamics.model")
	if model=="time_depend_Ginzburg_Landau":
		dt = _udf_.get("dynamics.time_depend_Ginzburg_Landau.delta_t")
	elif model=="visco_elastic":
		dt = _udf_.get("dynamics.visco_elastic.delta_t")
	rec = _udf_.get("simulation_conditions.total_records")
	inter = _udf_.get("simulation_conditions.output_interval_steps")

	return dt,rec,inter

def get_volume(_udf_):
	nx = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nx")
	ny = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Ny")
	nz = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nz")
	dx = _udf_.get("simulation_conditions.mesh_conditions.mesh_size.dx")
	dy = _udf_.get("simulation_conditions.mesh_conditions.mesh_size.dy")
	dz = _udf_.get("simulation_conditions.mesh_conditions.mesh_size.dz")
	return nx*ny*nz*dx*dy*dz

def reaction_estimation(_udf_, reaction_id):
	model = _udf_.get("reaction.reaction[].model",[reaction_id])
	
	dt,total,inter = get_time_obj(_udf_)
	
	if model=="first_order":
		r = _udf_.get("reaction.reaction[].first_order.reactant", [reaction_id] )
		p = _udf_.get("reaction.reaction[].first_order.product", [reaction_id] )
		rate = _udf_.get("reaction.reaction[].first_order.reaction_rate", [reaction_id] )
		
		r_phi = _udf_.get("component_properties.components[].volume_fraction",[r])
		p_phi = _udf_.get("component_properties.components[].volume_fraction",[p])
		r_name =  _udf_.get("component_properties.components[].name",[r])
		p_name =  _udf_.get("component_properties.components[].name",[p])

		title_string = "reaction estimation : " + str(reaction_id)
		tag_list = ["record", "reactant:"+r_name, "product:"+p_name]
		data_list = [[0],[r_phi],[p_phi]]
		k = rate*dt*inter
		for rec in range(1,total):
			dp = r_phi * ( 1.0 - exp(-k*rec) )
			data_list[0].append(rec)
			data_list[1].append(r_phi-dp)
			data_list[2].append(p_phi+dp)

		gnuplot.plot(data=data_list,title=title_string,labels=tag_list)
		
		for i in range(len(data_list[0])):
			print(data_list[0][i],data_list[1][i],data_list[2][i])
		
	elif model=="second_order":
		print("under constraction")
	elif model=="chain_growth":
		print("under constraction")
	else:
		print("invalid reaction type")
'''
def OLD_volume_fraction(_udf_):
	title_string = "volume fractions"
	tag_list = ["record"]
	dataList = [[]]
	total = _udf_.totalRecord()
	comps = _udf_.size("output.components[]")
	size  = _udf_.size("output.components[].density[]",[0])
	print total,comps,size
	
	for c in range(comps):
		tag_list.append("#"+str(c))
		dataList.append([])
		
	for rec in range(total):
		_udf_.jump(rec)
		dataList[0].append(rec)
		for c in range(comps):
			sum = 0.0
			for s in range(size):
				print s, sum, _udf_.get("output.components[].density[]",[c,s])
				sum = sum + _udf_.get("output.components[].density[]",[c,s])
			dataList[c+1].append(sum/size)
			
	gnuplot.plot(data=dataList,title=title_string,labels=tag_list)
'''
