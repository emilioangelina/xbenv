# import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats



def get_polar_coord_props(don, prec, prop): 
	np.seterr(divide='ignore', invalid='ignore')
	# Get the intervals (they should be the same used in the 2d binned scatter plots?)
	counts, xbins, ybins = np.histogram2d(don['Angle'], don['distance'], bins=(prec,prec)) 
	sums, _, _ = np.histogram2d(don['Angle'], don['distance'], weights=don[prop], bins=(xbins,ybins)) 

	# Average the charge density-based property values whitin 2d bins squares (angle bin x distance bin)
	prop_avg = np.nan_to_num(sums/counts)

	# Sum up togheter the charge density-derived properties for each complex in the reference system
	don_props_sum_avg = get_averaged_sum(don, prop_avg, ybins, xbins, 'L1overL3')
	return xbins, ybins, prop_avg ,don_props_sum_avg

def get_averaged_sum(all_ref_local_props, bin_avg_prop, dist_bins, angle_bins, property):
    all_pdbids = []
    all_global_props = []
    for group, frame in all_ref_local_props.groupby(['pdb id']):
        global_prop = 0
        for index, row in frame.iterrows():
            avg_local_prop = coord2bin(row, dist_bins, angle_bins, bin_avg_prop)
            global_prop += avg_local_prop
        all_pdbids.append(group)
        all_global_props.append(global_prop)
    return pd.DataFrame(all_global_props, index = all_pdbids, columns = [property])	


# find bin-averaged property for the matched point in the reference system 
def coord2bin(ref_props, dist_bins, angle_bins, bins_avg_props):
	val = 0.
	dist_condition = (ref_props['distance'] >= dist_bins[0]) and (ref_props['distance'] <= dist_bins[-1])
	ang_condition = (ref_props['Angle'] >= angle_bins[0]) and (ref_props['Angle'] <= angle_bins[-1])
	if (dist_condition & ang_condition):  
		i = np.where(ref_props['distance'] >= dist_bins[:-1])[0][-1] 
		j = np.where(ref_props['Angle'] >= angle_bins[:-1])[0][-1] 
		val = bins_avg_props[j,i]
		if (val == 0.):
			x = bins_avg_props[:,i]	
			val = np.mean(x[np.nonzero(x)])
		return val
	if (dist_condition & (not ang_condition)):
		i = np.where(ref_props['distance'] >= dist_bins[:-1])[0][-1] 
		x = bins_avg_props[:,i]	
		val = np.mean(x[np.nonzero(x)])
		return val
	if ((not dist_condition) & ang_condition):
		j = np.where(ref_props['Angle'] >= angle_bins[:-1])[0][-1] 
		x = bins_avg_props[j,:]	
		val = np.mean(x[np.nonzero(x)])
		return val
	if ((not dist_condition) & (not ang_condition)):
		return 0

def predict_bcp(datadir, bcp_props):
	with open(datadir+'/parameters.npy', 'rb') as f:
		# load model parameters 
		w0 = np.float((np.load(f)))
		w = np.load(f)
		mean = np.load(f)
		std = np.load(f)
		feature_names = np.load(f)
		# build handcrafted features 	

	pocket_column_names = bcp_props.index.values.tolist()
	bcp_props['ho hbond'] = ((bcp_props['ho dist'] <= 2.40) & (bcp_props['ho don_angle'] >= 100))*1
	bcp_props['hx acc_angle_center'] = np.abs(bcp_props['hx acc_angle'] - 90)	
	x = np.array(bcp_props[feature_names].values)

	# Scale data 
	x_scaled = (x - mean)/std 
	import math
	def sigmoid(x):
	    return 1 / (1 + pow(math.e, -x))
 
	result = 0
	result += w0
	for i in range(0, 5):
	    result += x_scaled[i] * w[i]
	result = sigmoid(result)
	return result		


def predict_prop(path, datadir, all_ref_local_props, bin_avg_prop, all_ref_global_props, dist_bins, angle_bins,prob_cutoff, verbose):
    import math 
    import glob
    path2 = path + "/geom_props"	
    qfiles = glob.glob(path2+'/*/dist_angles_H.txt')
    diff = []	
    diff_sq = 0
    predictions = {}
    for f in qfiles:	
        pdbid = f.split('/')[3]
        poses = pd.read_csv(f)	
        poses_O = pd.read_csv(path2+'/'+pdbid+'/dist_angles_O.txt') # to check that X-bond is actually formed
        best_global_prop = 0
        best_pose = None
        avg_global_prop = []
        for group, frame in poses.groupby('pose'):
            global_prop = 0
            for index_q, row_q in frame.iterrows():
                frame_O = poses_O.groupby('pose').get_group(row_q['pose'])[['o name','xb_angle']]
                XB = frame_O[(frame_O['o name'] == "O") & (frame_O['xb_angle'] >= 130)]
                if XB.empty:
#                        print('pose {} from {}, no X-bond found'.format(row_q['pose'], pdbid))
                        break
                prob = predict_bcp(datadir, row_q)
                if prob >= prob_cutoff:
                    row_q_tmp = pd.Series({'distance' : row_q['hx dist'], 'Angle' : row_q['hx acc_angle']}) 
                    avg_local_prop = coord2bin(row_q_tmp, dist_bins, angle_bins, bin_avg_prop)
                    global_prop += avg_local_prop
            if global_prop <= best_global_prop:
                best_global_prop = global_prop
                best_pose = (group,frame)
            predictions[pdbid] = best_global_prop
	print(' ')
	print('structure: {}, best_pose: {}'.format(pdbid, best_pose[0]))
        for index_qb, row_qb in best_pose[1].iterrows():	
            prob = predict_bcp(datadir, row_qb)
	    row_q_tmp = pd.Series({'distance' : row_qb['hx dist'], 'Angle' : row_qb['hx acc_angle']}) 
            avg_local_prop = coord2bin(row_q_tmp, dist_bins, angle_bins, bin_avg_prop)
	    if prob > prob_cutoff:
                print('H at_name: {}, H res_name: {}, H res_num: {}, BCP probability: {}, BCP lamba1/lambda3: {}'.format(row_qb['h name'],row_qb['h resName'],row_qb['h resNum'], prob,avg_local_prop))
    return 





