from UDFManager import *
from sys import *
import os
import math


def getList(_udf_):
	ni		=	_udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nx")
	nj		=	_udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Ny")
	nk		=	_udf_.get("simulation_conditions.mesh_conditions.number_of_mesh.Nz")
	dx		=	_udf_.get("simulation_conditions.mesh_conditions.mesh_size.dx")
	dy		=	_udf_.get("simulation_conditions.mesh_conditions.mesh_size.dy")
	dz		=	_udf_.get("simulation_conditions.mesh_conditions.mesh_size.dz")
	xlist	=	[0.0, ni*dx]
	ylist	=	[0.0, nj*dy]
	zlist	=	[0.0, nk*dz]
	divList	=	[ni, nj, nk]
	return xlist, ylist, zlist, divList

def setDraw(_udf_, meshf, contour_color, value, max, min):
	if contour_color == "blue_red":
		meshf.isovalue(value, _udf_.ccolor(1))
	elif contour_color == "rainbow":
		meshf.isovalue(value, [0,0,1,1,-240])
	else:
		meshf.isovalue(value, [0,0,0,1,1,1,1,1])

	if max > min:
		_udf_.crange(min, max)
	else:
		_udf_.crange(-1.0, 1.0)

def setSurface(mf, divList):
	for i in range(divList[0]):
		mf.set( [i, divList[1], divList[2]], mf.get([i,0,0]) )
		mf.set( [i, divList[1], 0], mf.get([i,0,0]) )
		for j in range(divList[1]):
			mf.set( [i,j,divList[2]], mf.get([i,j,0]) )

	for j in range(divList[1]):
		mf.set( [divList[0], j, divList[2]], mf.get([0,j,0]) )
		mf.set( [0, j, divList[2]], mf.get([0,j,0]) )
		for k in range(divList[2]):
			mf.set( [divList[0],j,k], mf.get([0,j,k]) )

	for k in range(divList[2]):
		mf.set( [divList[0], divList[1], k], mf.get([0,0,k]) )
		mf.set( [divList[0], 0, k], mf.get([0,0,k]) )
		for i in range(divList[0]):
			mf.set( [i,divList[1],k], mf.get([i,0,k]) )

	mf.set( [divList[0],divList[1],divList[2]], mf.get([0,0,0]) )
	


def getMaxMin(max, min, value):
	if max < value:
		max = value
	if min > value:
		min = value
	return max, min
	
def mean(_udf_, phiList):
	vs	=	0.0
	for p in range(len(phiList)):
		vs	=	vs + _udf_.get("output.components[].volume_fraction", [phiList[p],p])
	return vs
				
def density(_udf_, phiList, frame_attr, contour_color, max, min):
	xlist, ylist, zlist, divList = getList(_udf_)
	coord_list_list = [xlist, ylist, zlist]
	list = [divList[0], divList[1], divList[2]]
	mf			=	_udf_.meshfield("regular", coord_list_list, list)
	temp_max	=	0.0
	temp_min	=	1.0
	r = 0

	for k in range(divList[2]):
		for j in range(divList[1]):
			for i in range(divList[0]):
				v	=	0.0
				for p in range(len(phiList)):
					v	=	v + _udf_.get("output.components[].density[]",[phiList[p],r])
				temp_max, temp_min	=	getMaxMin(temp_max, temp_min, v)
				mf.set([i,j,k], v)
				r = r+1

	setSurface(mf, divList)

	if(max<0.0):
		max = temp_max
	if(min<0.0):
		min = temp_min

	setDraw(_udf_, mf, contour_color, mean(_udf_, phiList), max, min)
	mf.draw(frame=frame_attr, subdivision=2,iso_side_surface=1)


def mises(sxx, sxy, sxz, syy, syz, szz):
	xx = sxx
	yy = syy
	zz = szz
	j2 = -(xx*yy + yy*zz + zz*xx) + sxy*sxy + sxz*sxz + syz*syz
	return 3.0*math.sqrt(j2)

def stress(_udf_, frame_attr, contour_color, th, type, index):
	xlist, ylist, zlist, divList = getList(_udf_)
	coord_list_list = [xlist, ylist, zlist]
	list = [divList[0], divList[1], divList[2]]
	mf			=	_udf_.meshfield("regular", coord_list_list, list)

	max	=	-1e9
	min	=	1.0e9
	r = 0
	for k in range(divList[2]):
		for j in range(divList[1]):
			for i in range(divList[0]):
				if type=="mises":
					sxx	=	_udf_.get("output.components[].shear_stress[].xx",[index,r])
					sxy	=	_udf_.get("output.components[].shear_stress[].xy",[index,r])
					sxz	=	_udf_.get("output.components[].shear_stress[].xz",[index,r])
					syy	=	_udf_.get("output.components[].shear_stress[].yy",[index,r])
					syz	=	_udf_.get("output.components[].shear_stress[].yz",[index,r])
					szz	=	_udf_.get("output.components[].shear_stress[].zz",[index,r])
					v	=	mises(sxx, sxy, sxz, syy, syz, szz)
				elif type=="xx":
					v	=	_udf_.get("output.components[].shear_stress[].xx",[index,r])
				elif type=="xy":
					v	=	_udf_.get("output.components[].shear_stress[].xy",[index,r])
				elif type=="xz":
					v	=	_udf_.get("output.components[].shear_stress[].xz",[index,r])
				elif type=="yy":
					v	=	_udf_.get("output.components[].shear_stress[].yy",[index,r])
				elif type=="yz":
					v	=	_udf_.get("output.components[].shear_stress[].yz",[index,r])
				elif type=="zz":
					v	=	_udf_.get("output.components[].shear_stress[].zz",[index,r])
				elif type=="bulk":
					v	=  _udf_.get("output.components[].bulk_stress[]",[index,r])
				max, min	=	getMaxMin(max, min, v)
				mf.set([i,j, k], v)
				r = r+1

	setSurface(mf, divList)


	setDraw(_udf_, mf, contour_color, th, max, min)
	mf.draw(frame=frame_attr, subdivision=2,iso_side_surface=1)
