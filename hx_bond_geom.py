import glob
import sys
import Atoms 
import numpy as np
import math
import pandas as pd


# constants
HBOND_DIST_THRESHOLD=4.0
XBOND_DIST_THRESHOLD=3.7
XBOND_ANGLE_THRESHOLD=126.
HALOGEN = ['CL', 'Cl']

# functions 
def dist(a,b):
    return math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2)

def vector(a, b):
    return [b.x-a.x, b.y-a.y, b.z-a.z]

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

ext = sys.argv[1]
path = sys.argv[2]
lst = glob.glob(path+'/*.'+ext)
lst.sort()


ocolumns = ['pdb id', 'pose', 'x name', 'o name', 'resName', 'resNum', 'dist', 'xb_angle']
hcolumns = ['pdb id', 'pose', 'x name', 'h name', 'h resName', 'h resNum', 'hx dist', 'hx acc_angle','hx don_angle','o resNum', 'ho dist', 'ho don_angle', 'ox dist']

df_xbonds = open(path+'/dist_angles_O.txt', 'w')
df_xbonds.write(','.join(map(str,ocolumns))+'\n')

df_hbonds = open(path+'/dist_angles_H.txt', 'w')
df_hbonds.write(','.join(map(str,hcolumns))+'\n')


for f in lst:
	recName=path[-4:]
	pose=f.split('_')[-1].split('.')[0]
	mol=Atoms.pdb2mol(f) 

	# Find halogen atoms 	
	x_atoms = [atom for atom in mol if atom.element.strip() in HALOGEN if atom.resName.strip() not in ['CL', 'Cl']]
	num_x_atoms = len(x_atoms)
	halogens={}
	for x in x_atoms:
		hydrogens=[]
		oxygens=[]
		for atom in mol:  
			if (atom.element.strip() == 'H'):
				d = dist(x,atom)
				if d <= (HBOND_DIST_THRESHOLD):
					hydrogens.append(atom)
			if (atom.element.strip() == 'O') and (atom.resName != 'HOH'):
				d = dist(x,atom)
				if d <= (XBOND_DIST_THRESHOLD):
					oxygens.append(atom)
			if atom.resName == x.resName:
				if atom.element.strip() == 'C':
					if dist(x,atom) < 2.5:
						car = atom

		
		if car:	
			v1 = vector(car, x)				 
			OXs=[]
			for o in oxygens:
				v2 = vector(o, x)				 
				ang=math.degrees(angle(v1,v2))
				if (ang > XBOND_ANGLE_THRESHOLD) and (o.name.strip() == 'O'):
					ox = [o, dist(x,o), ang]
					OXs.append(ox)

			if OXs:			
				v1 = vector(x, car)				 
				HXs=[]		
				for h in hydrogens:
					v2 = vector(x, h)
					hb_acc_angle=math.degrees(angle(v1,v2))
					for atom in mol:
						if ((dist(h,atom) < 1.5) and (atom.resName == h.resName) and (atom.element.strip() != 'H')):
							donor = atom 				
							v3 = vector(h, donor)
							hb_don_angle=math.degrees(angle(vector(h,x),v3))
							dho=[] 
							best_ang = 0.
							for o in OXs:
								if o[2] > best_ang:
									best_ang = o[2]
									best_acc = o
							donoh_ang=math.degrees(angle(vector(h,best_acc[0]),v3))
							hx = [h,dist(x,h),hb_acc_angle,hb_don_angle, best_acc[0].resNum,dist(best_acc[0],h),donoh_ang, best_acc[1]] # fix 
							HXs.append(hx)
				halogens[x.name.strip()] = [OXs, HXs]				
				num_oxygens = len(np.asarray(halogens.values()).flatten()[0])
				num_hydrogens = len(np.asarray(halogens.values()).flatten()[1])
				xbonds = np.zeros((num_oxygens, len(ocolumns)),dtype=object)
				hbonds = np.zeros((num_hydrogens, len(hcolumns)),dtype=object)
				df_hbonds = pd.DataFrame(columns=hcolumns)
				for key, value in halogens.items():
					for i,o in enumerate(value[0]):
						xbonds[i,:] = np.array((recName, int(pose), key, o[0].name.strip(), o[0].resName.strip(), int(o[0].resNum), o[1], o[2]), dtype=object)	
#			
					for j,h in enumerate(value[1]):
						hbonds[j,:] = np.array((recName, int(pose), key, h[0].name.strip(), h[0].resName.strip(), h[0].resNum, h[1], h[2], h[3], h[4], h[5], h[6], h[7]), dtype=object)						

				df_xbonds = pd.DataFrame(xbonds, columns=ocolumns)						
				df_xbonds.to_csv(path+'/dist_angles_O.txt',mode='a', header=False)			

				df_hbonds = pd.DataFrame(hbonds, columns=hcolumns)			
				df_hbonds.to_csv(path+'/dist_angles_H.txt',mode='a', header=False)			
		



