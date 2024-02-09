import sys, os
import numpy as np
from UDFManager import *

def findfiles( directory, pattern ):
	objects = os.listdir(directory)
	ret = []
	for o in objects:
		if pattern in o:
			ret.append(o)
	return ret

files = findfiles(sys.argv[1], sys.argv[2]) #get SUSHI files, 

for sushi_file in files:
	sobj = UDFManager(sushi_file)
	name = "VANILLA_"+sushi_file
	with open(name,'w') as f:
		f.write("\\include{\"VANILLA_def10.udf\"}")
	vobj = UDFManager(name)
	
	#monomers
	mono = sobj.get("SUSHIInput.monomers[]")
	for i in range(len(mono)):
		name = mono[i][0]
		vobj.put(name,"component_properties.segment[].name",[i])

	#components
	#only support polymer-block
	comps = sobj.get("SUSHIInput.components.polymers[]")
	ids = []
	vfs = sobj.get("SUSHIInput.volume_fractions.polymer_volume_fractions[]")
	for i in vfs:
		id = i[0]
		v  = i[1]
		if v>0.0:
			ids.append((id,v))
	c = 0
	c0= 0
	for id,v in ids:
		comp = comps[id]
		totalN = 0.0
		Nlist = []
		for b in range(len(comp[1])):
			name = "polymer"+str(id)+"block"+str(b)
			seg  = comp[1][b][0]
			N    = comp[1][b][1]
			vobj.put(name,"component_properties.components[].name",[c])
			vobj.put(v,"component_properties.components[].volume_fraction",[c])
			vobj.put(seg,"component_properties.components[].segment",[c])
			vobj.put(N,"component_properties.components[].degree_of_polymerization",[c])
			vobj.put(name,"output.components[].name",[c])
			totalN = totalN + N
			Nlist.append(N)
			c = c+1

		for b in range(len(comp[1])):
			vol = v*Nlist[b]/totalN
			vobj.put(vol,"output.components[].volume_fraction",[c0])
			c0 = c0+1

	#chi_parameter
	#chi = sobj.get("SUSHIInput.chi_parameters[]")
	#i = 0
	#for c in chi:
	#	vobj.put(c[0],"component_properties.chi_parameter[].segment_name1",[i])
	#	vobj.put(c[1],"component_properties.chi_parameter[].segment_name2",[i])
	#	vobj.put(c[2],"component_properties.chi_parameter[].parameter",[i])
	#	i = i+1
	
	#mesh
	#only support REGULAR-3D
	nx = int(sobj.get("SUSHIInput.mesh.axes[0].values[2]"))
	ny = int(sobj.get("SUSHIInput.mesh.axes[1].values[2]"))
	nz = int(sobj.get("SUSHIInput.mesh.axes[2].values[2]"))
	sobj.jump(sobj.totalRecord()-1)
	dx = sobj.get("MeshData.position[].x",[1])
	dy = sobj.get("MeshData.position[].y",[nx])
	dz = sobj.get("MeshData.position[].z",[nx*ny])
	vobj.put(nx,"simulation_conditions.mesh_conditions.number_of_mesh.Nx")
	vobj.put(ny,"simulation_conditions.mesh_conditions.number_of_mesh.Ny")
	vobj.put(nz,"simulation_conditions.mesh_conditions.number_of_mesh.Nz")
	vobj.put(dx,"simulation_conditions.mesh_conditions.mesh_size.dx")
	vobj.put(dy,"simulation_conditions.mesh_conditions.mesh_size.dy")
	vobj.put(dz,"simulation_conditions.mesh_conditions.mesh_size.dz")

	#density
	phi = np.array(sobj.get("SUSHIOutput.phi.value[]")).T
	#print(phi.shape)
	for p in range(len(phi)):
		#print("phip",phi[0][p])
		for k in range(nx*ny*nz):
			#print("phipk",phi[0][p][k])
			vobj.put(phi[p][0][k],"output.components[].density[]",[p,k])

	vobj.write()
