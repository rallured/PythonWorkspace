import numpy as np
"""
#modified by kov
#x,y,z 3 vectors with coordinates of points (same number of elements).
#return value [A,B,C] of plane Ax + By + C = z
#planesurf is a vector with z of plane points

#compute the average surface, calculate statistical indicator

  # M. Katz 1/26/04
# IDL function to perform a least-squares fit a plane, based on
# Ax + By + C = z
#
# ABC = plane_fit(x, y, z, error=error)
"""

def plane_fit (x, y, z):

    if len(y) != len(x):
        raise ValueError, "wrong number of elements for x and y."
    if len(z) != len(x):
        raise ValueError, "wrong number of elements for z."

    tx2 = (x**2).sum()
    ty2 = (y**2).sum()
    txy = (x*y).sum()
    tx = x.sum()
    ty = y.sum()
    N = len(x)

    a = np.matrix([[tx2, txy, tx],
        [txy, ty2, ty],
        [tx, ty, N ]])

    b = np.array([(z*x).sum(), (z*y).sum(), (z).sum()])
    out = np.array(a.I * b[:,None])
    #plane=(out[0]*x+out[1]*y+out[2])

    return out

if __name__=='__main__':
    #plane parameters
    a=555.5
    b=-3
    c=-15

    #grid
    npx=100
    npy=100

    xp=np.arange(npx)
    yp=np.arange(npy)
    x,y=map(np.ndarray.flatten, np.meshgrid(xp,yp))
    z=x*a+y*b+c

    print "starting value:",a,b,c
    print " fit results: ",plane_fit(x,y,z)

  