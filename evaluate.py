import subprocess
import sys
import map2densityFunc_single as mapfunc_single
import numpy as np
import pandas as pd 
import shutil

class Xbond:
	def __init__(self, datadir, format, structures):
		self.datadir = datadir
		self.format = format
		self.structures = structures
		self.path = str(format)+'/'+str(structures) 

	def get_data(self):
		pass	

	def get_poses(self): # extract docking poses from dlg file 
		subprocess.call(["bash", "extract_poses.sh", self.format, self.structures, self.datadir])


	def get_geom_props(self):  # print H-bond/X-bond distances and angles 
			subprocess.call(["bash", "print_geom_batch.sh", self.format, self.structures, self.datadir])


	def get_density(self): # get densities
		# load atoms properties from the reference system 
		don = pd.read_csv(self.datadir+'/H-bond-donors_property_reference.csv')		
		don['L1overL3'] = don['L1']/don['L3']
		return don

	def map_to_density_single(self, verbose=True):
		don = self.get_density()
		fin=open(self.datadir+'/fit_params.txt')
		params = fin.readlines()[0].rstrip('\n').split(',')
		nbins = params[1]
		x_bins, y_bins, don_l1overl3_avg, don_props_sum_avg = mapfunc_single.get_polar_coord_props(don, int(nbins), 'L1overL3')
		mapfunc_single.predict_prop(self.path, self.datadir, don, don_l1overl3_avg, don_props_sum_avg, y_bins, x_bins, float(params[0]), verbose)
	

if __name__ == "__main__":

	datadir = sys.argv[1] # directory with input data
	format = sys.argv[2] # pdbs or docking poses
	structures = sys.argv[3] # name of the folder with structures to parse in data/folder   

	xbond = Xbond(datadir, format, structures)
	if (format == "poses"):
		xbond.get_poses()   

	xbond.get_geom_props()
	xbond.map_to_density_single(verbose=True) 
	shutil.rmtree(format)


