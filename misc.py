import math
#from matplotlib import pyplot as plt

def gauss(x,sig):
    """x is arg, sig is spread"""
    return math.exp( -(x)**2 / (2*sig**2))

def smoothgauss(x, sig, smooth=0):
    x = (abs(x) > smooth) * (abs(x) - smooth)
    return math.exp( -(x)**2 / (2*sig**2))


def LJ(r, zero, repulse, attract):
    term = zero / r                 # sigma/r
    repulsive = term ** repulse     # (s/r)^12
    attractiv = term ** attract      # (s/r)^6
    return 4 * (repulsive - attractiv)

def MGL_LJ(r, rmin, n, m, smooth):
    """ As in AutoDock3.5 Manual """
    pass

#X = [i*0.1 for i in range(-50,50,1)]
#Y = []
#for x in X:
#    Y.append(gauss(x,3,2))
#    
#
#plt.plot(X,Y)
#plt.show()

class GPF():
    def __init__(self, filename):
        self.digest(filename)
        self.filename = filename
        self.set_odd_npts() # ZnGrid class depends on this

    def digest(self, filename):
        f = open(filename)
        self.maps = []
        for line in f:
            fields = line.split()
            if line.startswith('npts'):
                self.npts = tuple([int(x) for x in fields[1:4]])
            if line.startswith('gridcenter'):
                self.gridcenter = tuple([float(x) for x in fields[1:4]])
            if line.startswith('spacing'):
                self.spacing = float(fields[1])
            if line.startswith('receptor'):
                self.receptor = fields[1]
            if line.startswith('ligand_types'):
                self.ligand_types = []
                for field in fields:
                    if '#' in field:
                        break
                    self.ligand_types.append(field)
            if line.startswith('map'):
                self.maps.append(fields[1])
            if line.startswith('gridfld'):
                self.gridfld = fields[1]
        f.close()
        # missing: gridfld, elecmap, dmap, dieletric, smooth, receptor_types

    def set_odd_npts(self):
        self.npts = [int(2*math.floor(self.npts[i]/2)+1) for i in range(3)]

    def get_odd_npts(self):
        return [int(2*math.floor(self.npts[i]/2)+1) for i in range(3)]

class MapHeader():
    def __init__(self,filename):
        f = open(filename)
        for i in range(6):
            line = f.readline()
            fields = line.split()
            if fields[0] == 'GRID_PARAMETER_FILE':
                self.gpf = fields[1]
            if fields[0] == 'GRID_DATA_FILE':
                self.gridfld = fields[1]
            if fields[0] == 'MACROMOLECULE':
                self.receptor = fields[1]
            if fields[0] == 'NELEMENTS':
                self.npts = tuple([int(i) for i in fields[1:4]])
            if fields[0] == 'SPACING':
                self.spacing = float(fields[1])
            if fields[0] == 'CENTER':
                self.gridcenter = tuple([float(i) for i in fields[1:4]])
        f.close()
        self.set_odd_npts() # ZnGrid class depends on this

    def set_odd_npts(self):
        self.npts = [int(2*math.floor(self.npts[i]/2)+1) for i in range(3)]


def key_from_item(your_dict, target):
    out = []
    for key in your_dict:
        if your_dict[key] == target:
            out.append(key)
    if len(out) == 1:
        return out[0]
    elif len(out) > 1:
        return out
    else:
        return None

def ki2kcal(ki):
    R = 1.9858775e-3
    T = 298.15
    kcal = math.log(ki)*R*T
    return kcal

