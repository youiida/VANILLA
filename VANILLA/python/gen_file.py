# -*- coding: utf-8 -*-
import os, re, shutil, sys
from UDFManager import *

def findFile( dir, pat ):
	ret = []
	cm = re.compile( pat )
	for file in os.listdir(dir):
		if cm.search( file ):
			ret.append( file )
	return ret

def gen_files(parameter_name, para_path, values):
	#make udfs-list to change paramters
	files = findFile(sys.argv[1], ".udf")
	for file in files:
		for value in values:
			nfname = parameter_name+str(value)+"_"+file
			shutil.copy(file, nfname)
			uobj = UDFManager(nfname)
			uobj.put(value, para_path)
			uobj.write()
		
## how to use
## $> python gen_file.py (directory that you want to edit udfs in)
##

gen_files("tb", "component_properties.segment[0].bulk_relaxation_time", [1])
gen_files("ts", "component_properties.segment[0].shear_relaxation_time", [1])




