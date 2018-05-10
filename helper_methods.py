# epnurbs.helper_methods module contains a few helper methods for 3d vector calculations

def subtract(u,v):
    return [ u[0]-v[0], u[1]-v[1], u[2]-v[2] ]

def dotproduct(u,v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    
def crossproduct(u, v):
    return [ u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0] ]

def length(u):
    from math import sqrt
    return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2])

def normalize(u):
    l = length(u)
    u[0] /= l
    u[1] /= l
    u[2] /= l
    return u

def distance(a,b):
    return length( (a[0]-b[0], a[1]-b[1], a[2]-b[2]) )

# finds the foot of the perpendicular from X to the plane <X-P, N>=0
# the method assumes that the normal vector n is normalised
# the foot is calculated as F = X - <X-P,N>N
def foot(x, p, n):
    s = (x[0]-p[0])*n[0] + (x[1]-p[1])*n[1] + (x[2]-p[2])*n[2]
    return [ x[0] - s*n[0], x[1] - s*n[1], x[2] - s*n[2] ]
