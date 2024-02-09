from UDFManager import *
from sys import *
import os
import math
import numpy.random as ran

def getIndex(_udf_, show_type, list):
	index = []
	if len(list)==0 :
		index.append(0)
		return index
	if(show_type=="segment_type"):
		polys = _udf_.size("component_properties.components[]")
		for p in range(polys):
			name = _udf_.get("component_properties.components[].segment",[p])
			for i in list:
				if i==name:
					index.append(p)
	return index

def make_cognac_density_file(_udf_, filename, lx, ly, lz):
	cogfile = os.path.join(_udf_.udfDirectory(),filename+".udf")
	oobj = open(cogfile, 'w')
	oobj.write("\\include{\"cognac90.udf\"}")
	oobj.close()
	cudf = UDFManager(cogfile)
	cudf.newRecord("dummy",0,1)
	cudf.newRecord("density_grid")
	
	x = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nx")
	y = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Ny")
	z = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nz")
	#put data grid density
	cudf.put(_udf_.udfFilename().split(".")[0], "Grid_Density.mesh.name")
	cudf.put("REGULAR","Grid_Density.mesh.type")
	
	cudf.put(0.0, "Grid_Density.mesh.axes[0].values[0]")
	cudf.put(lx, "Grid_Density.mesh.axes[0].values[1]")
	cudf.put(x, "Grid_Density.mesh.axes[0].values[2]")
	cudf.put(0.0, "Grid_Density.mesh.axes[1].values[0]")
	cudf.put(ly, "Grid_Density.mesh.axes[1].values[1]")
	cudf.put(y, "Grid_Density.mesh.axes[1].values[2]")
	cudf.put(0.0, "Grid_Density.mesh.axes[2].values[0]")
	cudf.put(lz, "Grid_Density.mesh.axes[2].values[1]")
	cudf.put(z, "Grid_Density.mesh.axes[2].values[2]")
	cudf.put(0, "Grid_Density.mesh.index_rule[0]")
	cudf.put(1, "Grid_Density.mesh.index_rule[1]")
	cudf.put(2, "Grid_Density.mesh.index_rule[2]")
	
	#periodic only
	cudf.put("PERIODIC", "Grid_Density.boundary_condition.conditions[0].axis_conditions[0]")
	cudf.put("PERIODIC", "Grid_Density.boundary_condition.conditions[1].axis_conditions[0]")
	cudf.put("PERIODIC", "Grid_Density.boundary_condition.conditions[2].axis_conditions[0]")
	
	cudf.put(1.0,"Statistics_Data.Density.Instantaneous")
	cudf.put(lx, "Structure.Unit_Cell.Cell_Size.a")
	cudf.put(ly, "Structure.Unit_Cell.Cell_Size.b")
	cudf.put(lz, "Structure.Unit_Cell.Cell_Size.c")
	cudf.put(90, "Structure.Unit_Cell.Cell_Size.alpha")
	cudf.put(90, "Structure.Unit_Cell.Cell_Size.beta")
	cudf.put(90, "Structure.Unit_Cell.Cell_Size.gamma")

	cudf.write()
	return cudf

def export_cognac(_udf_, filename, lx, ly, lz):
	cudf = make_cognac_density_file(_udf_, filename, lx, ly, lz)
	
	N = _udf_.size("component_properties.components[]")
	for i in range(N):
		name = _udf_.get("component_properties.components[].name", [i])
		cudf.put(name, "Grid_Density.atom_name[]", [i])
	
	cudf.put("Phi", "Grid_Density.phi.name")
	cudf.put(N, "Grid_Density.phi.num_of_component")
	
	v = []
	for i in range(N):
		val = _udf_.get("output.components[].density[]", [i])
		v.append(val)
	
	size = len(val)
	for i in range(size):
		for j in range(N):
			cudf.put(v[j][i], "Grid_Density.phi.value[].comp[]", [i,j])
	
	cudf.jump(-1)
	for i in range(N):
		name = _udf_.get("component_properties.components[].name", [i])
		cudf.put(name, "Set_of_Molecules.molecule[].atom[].Atom_Name",[0,i])
		cudf.put(name, "Set_of_Molecules.molecule[].atom[].Atom_Type_Name",[0,i])
	
	cudf.write()


def pdist2(p1,p2,lx,ly,lz):
	dist = [abs(p1[0]-p2[0]), abs(p1[1]-p2[1]), abs(p1[2]-p2[2])]
	if dist[0]>lx*0.5:
		dist[0] = lx-dist[0]
	if dist[1]>ly*0.5:
		dist[1] = ly-dist[1]
	if dist[2]>lz*0.5:
		dist[2] = lz-dist[2]
	return dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]


def export_filler(_udf_, filename, flist, beta, lx, ly, lz):
	cudf = make_cognac_density_file(_udf_, filename, lx, ly, lz)
	#filler_list = [[name,vol,size,bias]]
	#flist = [["Filler-A", 0.3, 10.0, 0],["Filler-B", 0.1, 20.0, 0]]

	N = _udf_.size("component_properties.components[]")
	fnames = []
	for i in range(len(flist)):
		name = flist[i][0]
		flag = True
		for j in fnames:
			if j==name:
				flag=False
		if flag==True:
			fnames.append(name)
	M = len(fnames)

	for i in range(N):
		name = _udf_.get("component_properties.components[].name", [i])
		cudf.put(name, "Grid_Density.atom_name[]", [i])
	for i in range(M):
		name = fnames[i]
		cudf.put(name, "Grid_Density.atom_name[]", [N+i])
		
	cudf.put("Phi", "Grid_Density.phi.name")
	cudf.put(N+M, "Grid_Density.phi.num_of_component")

	v = []
	for i in range(N):
		val = _udf_.get("output.components[].density[]", [i])
		v.append(val)
	
	size = len(val)
	for i in range(size):
		for j in range(N):
			cudf.put(v[j][i], "Grid_Density.phi.value[].comp[]", [i,j])
	
	cudf.jump(-1)
	for i in range(N):
		name = _udf_.get("component_properties.components[].name", [i])
		cudf.put(name, "Set_of_Molecules.molecule[].atom[].Atom_Name",[0,i])
		cudf.put(name, "Set_of_Molecules.molecule[].atom[].Atom_Type_Name",[0,i])
	for i in range(M):
		name = fnames[i]
		cudf.put(name, "Set_of_Molecules.molecule[].atom[].Atom_Name",[0,N+i])
		cudf.put(name, "Set_of_Molecules.molecule[].atom[].Atom_Type_Name",[0,N+i])


	list = []
	Nx = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nx")
	Ny = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Ny")
	Nz = _udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nz")


	#generate filler
	#flist = [["Filler-A", 0.3, 10.0, 0],["Filler-B", 0.1, 20.0, 0]]
	for stat in flist:
		name = stat[0]
		vol  = stat[1]
		size = stat[2]
		bias = stat[3]
		onev = 4*math.pi*size*size*size/3.0
		num  = int(vol*lx*ly*lz/onev)
		c    = 0
		#count= 0
		print(stat, num*onev/(lx*ly*lz), num)
		while(c<num):
			flag = False
			rx = ran.rand()
			ry = ran.rand()
			rz = ran.rand()
			
			ix = int( Nx*rx )
			iy = int( Ny*ry )
			iz = int( Nz*rz )
			val = 0.5 - _udf_.get("output.components[].density[]", [bias, Nz*Ny*ix+Ny*iy+iz])
			if val<0.0:
				flag = True
			else:
				prob = ran.rand()
				if prob < math.exp(-beta*val):
					flag = True
			if flag == True:
				px = lx*rx
				py = ly*ry
				pz = lz*rz
				for i in range(len(list)):
					pos = list[i][0]
					s   = list[i][1]
					d2  = pdist2([px,py,pz], pos, lx, ly, lz)
					if d2 < (size+s)*(size+s):
						flag = False
				if flag == True:
					list.append([[px,py,pz], size, name])
					c = c+1

	#write filler-phi
	phi = []
	dx  = lx / float(Nx)
	dy  = ly / float(Ny)
	dz  = lz / float(Nz)
	
	for f in range(len(fnames)):
		v = Nx*Ny*Nz*[0]
		for l in list:
			if l[2] == fnames[f]:
				fp    = l[0]
				fs    = l[1]
				size2 = fs*fs
				total = 0
				for i in range(Nx):
					for j in range(Ny):
						for k in range(Nz):
							d2 = pdist2([i*dx, j*dy, k*dz], fp, lx, ly, lz)
							if d2<size2:
								v[Nz*Ny*i+Nz*j+k] = 1.0
								total = total + 1
		for i in range(len(v)):
			cudf.put(v[i], "Grid_Density.phi.value[].comp[]", [i,N+f])
	cudf.write()

