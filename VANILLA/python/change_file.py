# -*- coding: utf-8 -*-
import os,re, sys
from UDFManager import *



def findFile( dir, pat ):
	ret = []
	cm = re.compile( pat )
	for file in os.listdir(dir):
		if cm.search( file ):
			ret.append( file )
	return ret

## how to use
## $> python change_file.py (directory that you want to edit udfs in)
##

list = findFile(sys.argv[1], ".udf")
#make udfs-list to change paramters

for name in list:
	print(name)
	uobj = UDFManager(name)
	# set parameters
	uobj.put(1.0,"component_properties.segment[0].interfacial_tension")
	uobj.put(1.0,"component_properties.segment[1].interfacial_tension")
	uobj.write()
