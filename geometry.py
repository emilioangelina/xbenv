import math 


def dist(a,b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def vector(a, b):
    return [b[0]-a[0], b[1]-a[1], b[2]-a[2]]

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))
    
def length(v):
    return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    return math.degrees(math.acos(dotproduct(v1, v2) / (length(v1) * length(v2))))


