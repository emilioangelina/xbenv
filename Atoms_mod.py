#!/usr/bin/python

import math
import misc
import Geom

class Atom():
    def __init__(self, x, y, z, atomnr_elemnt):
        # Add element and atomic nr
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.atomnr, self.element = atomnr_element_converter(
            atomnr_element)

    def dist(self, other_atom):
        return Geom.dist(self.getcoords(), other_atom.getcoords())
   
    def angle(self, atom_2, atom_3):
        """Calculate the angle between three atoms. This atom in the middle"""
        d12 = self.dist(atom_2)
        d13 = self.dist(atom_3)
        d23 = atom_2.dist(atom_3)
        #round. To avoid things like 1.000000001
        angle = math.acos(round((d12**2 + d13**2 - d23**2)/(2*d12*d13),7))
        return angle*180/math.pi

    def dihedral(self, b, c, d):
        """ Dihedral with this atom in the beggining """
        a = self.getcoords()
        b = b.getcoords()
        c = c.getcoords()
        d = d.getcoords()
        return Geom.dihedral(a,b,c,d)

    def getcoords(self):
        return (self.x, self.y, self.z)

    def setcoords(self, coords):
        self.x, self.y, self.z = coords

    def isbound(self, other_atom, cut_off_percent = .5):
        """ Depends on atomic number """
        threshold = cut_off_percent * (
            .5 * atomnr_vdw(self.atomnr) +
            .5 * atomnr_vdw(other_atom.atomnr))
        is_bound = self.dist(other_atom) < threshold
        return is_bound

        pass

class PDB(Atom):
    def __init__(self, line):
        self._parse_common(line)     # used here (PDB) and in PDBQT
        self._parse_specific(line)   # differs in PDBQT

    def getline(self):
        txt = self._print_common()    # no \n; PDB + PDBQT
        txt += self._print_specific() # differs in PDBQT
        return txt

    # blanks: [11:12], [20:21], [27:30], [66:76](pdb) or [66:68](pdbqt)

    def _parse_common(self, line):
        self.keyword     = line      [ 0: 6]     # ATOM or HETATM
        self.serial      = int(line  [ 6:11])    # atom id
        #                            [11:12]
        self.name        = line      [12:16]     # atom name
        self.altLoc      = line      [16:17]     # Alternate location
        self.resName     = line      [17:20]     # Residue name
        #                            [20:21] 
        self.chain       = line      [21:22]     # chain
        self.resNum      = int(line  [22:26])    # Residue number
        self.icode       = line      [26:27]     # ???
        #                            [27:30]
        self.x           = float(line[30:38])    # X
        self.y           = float(line[38:46])    # Y
        self.z           = float(line[46:54])    # Z
        self.occupancy   = float(line[54:60])    # Occupancy
        self.bfact       = float(line[60:66])    # Temperature factor

    def _parse_specific(self, line):
        self.element     = line      [76:79]     # Element mod
        self.frmlcharge  = line      [78:80]     # Charge as in PDB
        self.element = self.element.strip().upper()
        self.atomnr = atomnr_elemnt_converter(self.element)
        
    def _print_common(self):
        """ Characters [0:68]"""
        linestr = ''
        linestr += '%6s' % (self.keyword)
        linestr += '%5d' % (self.serial)
        linestr += ' ' 
        linestr += '%4s' % (self.name)
        linestr += '%1s' % (self.altLoc) 
        linestr += '%3s' % (self.resName)
        linestr += ' ' 
        linestr += '%1s' % (self.chain)
        linestr += '%4d' % (self.resNum)
        linestr += '%1s' % (self.icode)
        linestr += ' ' * 3 
        linestr += '%8.3f' % (self.x)
        linestr += '%8.3f' % (self.y)
        linestr += '%8.3f' % (self.z)
        linestr += '%6.2f' % (self.occupancy)
        linestr += '%6.2f' % (self.bfact)
        return linestr

    def _print_specific(self):
        """ PDB characters [68:80]"""
        linestr = ''
        linestr += ' ' * 10                 # [66:76]
        linestr += '%2s' % self.element     # [76:78]
        linestr += '%2s' % self.frmlcharge  # [78:80]
        linestr += '\n'
        return linestr

class PDBQT(PDB):
    def _parse_specific(self, line):
        """ PDBQT characters [68:79] """
        self.charge      = float(line[68:76])   # Charge
        self.atype       = line      [77:79]    # Atom type
        self.atype = self.atype.strip().upper()
        self.atomnr = AutoDockType_to_atomnr(self.atype)

    def _print_specific(self):
        """ PDBQT characters [68:79] """
        linestr =  ' ' * 2                      # [66:68]
        linestr += '%8.3f' % (self.charge)      # [68:76]
        linestr += ' ' * 1                      # [76:77]
        linestr += '%2s' % (self.atype)       # [77:79]
        linestr += '\n'
        return linestr

# Aux Functions

def centerMol(atom_list):
    center = [0.0, 0.0, 0.0]
    for atom in atom_list:
        center[0] += atom.x
        center[1] += atom.y
        center[2] += atom.z
    return [center[i]/len(atom_list) for i in [0,1,2]]

def pdbMissingRes(name):
    read1 = False
    read2 = False
    read3 = False
    f = open(name)
    output = []
    for line in f:
        # Part 1 - STOP conditions
        stop_reading = (
            len(line.replace(' ','')) < 10 or not (
            line.startswith('REMARK 465') or
            line.startswith('REMARK 570') or
            line.startswith('REMARK 610'))
            )
        if stop_reading:
            read1 = False
            read2 = False
            read3 = False
        if line.startswith('ATOM'):
            return output
        # Part 2 - READ
        if read1:
            resnum  = line[22:26]
        if read2:
            resnum  = line[20:24]
        if read3:
            resnum  = line[21:25]
        if read1 or read2 or read3:
            resname = line[15:18]
            chain   = line[19:20]
            output.append((int(resnum), resname.strip(' '), chain))
        # Part 3 - START conditions
        if "REMARK 465   M RES C SSSEQI" in line:
            read1 = True
        if "REMARK 470   M RES CSSEQI  ATOMS" in line:
            read2 = True
        if "REMARK 610   M RES C SSEQI" in line:
            read3 = True

def pdb2mol(name, inType = 'pdb'):
    atoms_list = []
    f = open(name)
    if inType == 'pdb':
        for line in f:
            if line.startswith('ATOM  ') or line.startswith('HETATM'):
                atoms_list.append(PDB(line))
    elif inType == 'pdbqt':
         for line in f:
            if line.startswith('ATOM  ') or line.startswith('HETATM'):
                atoms_list.append(PDBQT(line))
    else:
        raise RuntimeError('Unknown input type: Must be "pdb" or "pdbqt"')
    f.close()
    return atoms_list

def mol2pdb(name, atoms_list):
    f = open(name,'w')
    for atom in atoms_list:
        f.write(atom.getline())
    f.close()

def minMolDist(atom_listA, atom_listB):
    minDist = float('inf')
    for A in atom_listA:
        for B in atom_listB:
            dist = A.dist(B)
            minDist = min(dist, minDist)
    return minDist

def atomnr_elemnt_converter(x):
    D = {1:'H', 6:'C', 7:'N', 8:'O', 9:'F', 12:'MG', 15:'P', 16:'S', 30:'ZN'}
    
    if str(x).isdigit():
        return x, D[x] 
    else:
        return misc.key_from_item(D,x.strip().upper()), x

def AutoDockType_to_atomnr(atype):
    D = {'H':1, 'HD':1, 'HS':1, 'C':6, 'A':6, 'N':7, 'NA':7, 'NS':7, 'OA':8,
        'OS':8, 'F':9, 'MG':12, 'S':16, 'SA':16, 'CL':17, 'CA':20, 'MN':25,
        'FE':26, 'ZN':30, 'BR':35, 'I':53, 'G':6, 'J':6, 'P':15, 'HG':1}
    try:
        return D[atype.strip().upper()]
    except:
        print('STRANGE ATOM: %s' % atype)
        return -1 

def atomnr_vdw(atomnr):
    D = {1:2.0, 6:4.0, 7:3.5, 8:3.2, 9:3.1, 12:1.3, 15:4.2, 16:4.0, 17:4.1, 
        20:2.0, 25:1.3, 26:1.3, 30:1.5, 35:4.3, 53:4.7} 
    return D[atomnr]
