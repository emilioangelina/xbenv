import math
def dist(a,b):
    return math.sqrt(sum([(a[i]-b[i])**2 for i in range(min(len(a),len(b)))]))

def dihedral(a,b,c,d):
    """ Calculate dihedral considering a in the beggining"""
    v1 = [b[i] - a[i] for i in range(3)]
    v2 = [c[i] - b[i] for i in range(3)]
    v3 = [d[i] - c[i] for i in range(3)]
    temp = [dist((0,0,0),v2) * v1[i] for i in range(3)]
    y = dotProd(temp ,crossProd(v2,v3))
    x = dotProd(crossProd(v1,v2),crossProd(v2,v3))
    rad = math.atan2(y,x)
    return rad*(180/math.pi) 

def dotProd(a,b):
    N = min(len(a),len(b))
    return sum([a[i] * b[i] for i in range(N)])

def crossProd(a,b):
    """Pretty self-explanatory, this function bakes cookies"""
    normal_vect = [
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]]
    return normal_vect

def rawVec(a,b):
    N = min(len(a),len(b))
    return [b[i]-a[i] for i in range(N)]

