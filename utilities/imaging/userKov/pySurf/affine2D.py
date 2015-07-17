import numpy as np
from points import _angpos2
import math
from points import rotate_points,translate_points

## rotate_center -> points.rotate_points
## translate -> points.translate_points

def rot_center_func(theta,center=(0,0)):
    def func(x):
        return rotate_points(x,theta,center)
    return func
    
def translate_func(offset=(0,0)):
    def func(x):
        return translate(x,offset)
    return func
    
def rototrans_func(theta,center=(0,0),offset=(0,0)):
    """return a function that rotate by theta about center THEN translate by offset"""
    def func(x):
        if (x.shape[-1]==3):
            xy=x[:,0:2]
            return np.hstack([translate_points(rotate_points(xy,theta,center),offset),x[:,2,np.newaxis]]) 
        else:  
            return translate_points(rotate_points(x,theta,center),offset)
    return func
    
def find_rototrans(primary, secondary,verbose=False):
    """Return a function that can transform points from the first system to the second by means of a rototranslation. Also (not implemented yet, but kept for interface consistency with find_affine) return the matrix A of the transformation, that can be applied to a vector x with unpad(np.dot(pad(x), A)).
    primary and secondary are sets of points in format [Npoints, 2]. Transformation matrix A is [Ndim+1 x Ndim+1]."""
    
    ang1,r1,b1=_angpos2(primary[:,0:2])   #glass-> 1, mandrel-> 2
    ang2,r2,b2=_angpos2(secondary[:,0:2])
    mrot=(ang2-ang1).mean()
    bartrans=b2-b1
    #define the transformation function. I will transform primary data
    transform=rototrans_func(mrot,b1,bartrans)
    #if gmarkers.shape[0]: transform=a2d.find_affine(mmarkers,gmarkers)
    if verbose:
        print "stdv of markers distance errors from barycenter: %s"%((r1-r2).std())
        print "translation of barycenter (dx,dy,dist): %s"%(bartrans[0],bartrans[1],np.sqrt(np.sum(bartrans**2)))
        print 'rotation angle (degrees): %s +-%s'%(mrot*180./math.pi,(ang1-ang2).std()*180./math.pi)
        print "errors in markers position after rotations:"
        print secondary-transform(primary)
    return transform,[mrot,b1,bartrans] 
    
def find_affine(primary, secondary):
    """Return a function that can transform points from the first system to the second. Also return the matrix A of the transformation, that can be applied to a vector x with unpad(np.dot(pad(x), A)).
    primary and secondary are sets of points in format [Npoints, Ndim]. Transformation matrix A is [Ndim+1 x Ndim+1].
    """
    # Pad the data with ones, so that our transformation can do translations too
    n = primary.shape[0]
    pad = lambda x: np.hstack([x, np.ones((x.shape[0], 1))])
    unpad = lambda x: x[:,:-1]
    X = pad(primary)
    Y = pad(secondary)

    # Solve the least squares problem X * A = Y
    # to find our transformation matrix A
    A, res, rank, s = np.linalg.lstsq(X, Y)
    
    #note: the transformation to and from matrix is needed to work also if x is a single point.
    transform = lambda x: np.array(unpad(np.dot(pad(np.matrix(x)), A)))

    return transform,A   


if __name__=="__main__":
    primary = np.array([[40., 1160., 0.],
                        [40., 40., 0.],
                        [260., 40., 0.],
                        [260., 1160., 0.]])

    secondary = np.array([[610., 560., 0.],
                          [610., -560., 0.],
                          [390., -560., 0.],
                          [390., 560., 0.]])
    
    print "-----------"
    print "Test affine"             
    print "-----------"         
    transform,A= find_affine(primary, secondary)
    print "Starting points:"
    print primary
    print "Target:"
    print secondary
    print "Result:"
    print transform(primary)
    print "Max error:", np.abs(secondary - transform(primary)).max()
    
    
    print "-----------"
    print "Test rototrans"   
    print "-----------"
    secondary2=rotate_points(primary,math.pi/6,[-10,2])
    transform2,A2= find_rototrans(primary, secondary2,verbose=True)
    print "Result:"
    print transform2(primary)
   
